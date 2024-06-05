//! IonSolver supports execution on multiple compute nodes using the mpi protocol if compiled with the `multi-node` feature.
//! 
//! ## MPI Error codes
//! - **100**: Incompatible domain number and number of execution nodes

use mpi::{point_to_point, traits::*};
use ocl_macros::default_device;
use crate::*;

/// Run IonSolver on multiple compute nodes without multi-gpu support.
/// This function starts a single node with control over one `LbmDomain`.
pub fn run_node() {
    let universe = mpi::initialize().unwrap();
    let world: mpi::topology::SimpleCommunicator = universe.world();
    let size = world.size();
    let rank = world.rank();
    if rank == 0 {
        println!("IonSolver - © 2024\n");
    }
    println!("Launched Node {} of {}", rank, size);

    let mut cfg: LbmConfig;
    if rank == 0 { // If root, read and send config to all nodes
        cfg = LbmConfig::new();
        cfg.n_x = 256;
        cfg.n_y = 256;
        cfg.n_z = 256;
        cfg.d_x = 2;
        cfg.nu = cfg.units.si_to_nu(0.1);
        cfg.velocity_set = VelocitySet::D3Q19;
        // Graphics
        cfg.graphics_config.graphics_active = true;
        cfg.graphics_config.streamline_every = 8;
        cfg.graphics_config.vec_vis_mode = graphics::VecVisMode::U;
        //lbm_config.graphics_config.streamline_mode = true;
        cfg.graphics_config.axes_mode = true;
        cfg.graphics_config.q_mode = true;
        cfg.graphics_config.q_min = 0.00001;

        let ser_cfg = bincode::serialize(&cfg).unwrap(); 
        println!("Distributing LbmConfig...");
        for i in 1..size {
            world.process_at_rank(i).send(&ser_cfg[..]);
        }
    } else { // receive config from root node
        let (ser_cfg_bytes, _status) = world.any_process().receive_vec::<u8>();
        cfg = bincode::deserialize(&ser_cfg_bytes).unwrap();
    }

    // Validate domain numbers/world size
    let domain_numbers: u32 = cfg.d_x * cfg.d_y * cfg.d_z;
    if size as u32 != domain_numbers {
        println!("Domain number and node mismatch. {} domains and {} nodes. Aborting", domain_numbers, size);
        world.abort(100);
    }

    let mut domain = node_domain(&mut cfg, rank as u32); // Build domain for node
    println!("Build domain at Node {}", rank);
    world.barrier();
    if rank == 0 {
        println!("Beginning execution");
    }
    domain.node_initialize(&world);
    for _ in 0..1000 {
        domain.node_do_time_step(&world);
    }
}

/// Initialize an `LbmDomain` for a compute node
#[rustfmt::skip]
fn node_domain(cfg: &mut LbmConfig, rank: u32) -> LbmDomain {
    let n_d_x: u32 = (cfg.n_x / cfg.d_x) * cfg.d_x; // Validate and correct simulation size
    let n_d_y: u32 = (cfg.n_y / cfg.d_y) * cfg.d_y;
    let n_d_z: u32 = (cfg.n_z / cfg.d_z) * cfg.d_z;
    if n_d_x != cfg.n_x || n_d_y != cfg.n_y || n_d_z != cfg.n_z {
        println!( "Warning: Resolution {}, {}, {} not divisible by Domains: Overiding resolution.", cfg.n_x, cfg.n_y, cfg.n_z );
    }
    cfg.n_x = n_d_x;
    cfg.n_y = n_d_y;
    cfg.n_z = n_d_z;

    let x = (rank % (cfg.d_x * cfg.d_y)) % cfg.d_x; // Get domain coordinates
    let y = (rank % (cfg.d_x * cfg.d_y)) / cfg.d_x;
    let z =  rank / (cfg.d_x * cfg.d_y);

    LbmDomain::new( &cfg, default_device!(), x, y, z ) // Initialize and return domain
}

// Additional functionality needed for multi-node exectution
impl LbmDomain {
    /// Readies the LBM Simulation to be run.
    /// Executes `kernel_initialize` Kernels for every `LbmDomain` and fills domain transfer buffers.
    pub fn node_initialize(&mut self, world: &mpi::topology::SimpleCommunicator) {
        // the communicate calls at initialization need an odd time step
        self.t = 1;
        self.node_communicate_rho_u_flags(&world);
        self.enqueue_initialize().unwrap();
        self.node_communicate_rho_u_flags(&world);
        self.node_communicate_fi(&world);
        if self.cfg.ext_magneto_hydro && self.cfg.induction_range != 0 {
            self.node_communicate_qi(&world);
            self.enqueue_update_e_b_dyn().unwrap();
        }
        self.queue.finish().unwrap();
        self.t = 0;
    }

    /// Executes one LBM time step.
    /// Executes `kernel_stream_collide` Kernels for every `LbmDomain` and updates domain transfer buffers.
    /// Updates the dynamic E and B fields.
    fn node_do_time_step(&mut self, world: &mpi::topology::SimpleCommunicator) {
        // call kernel stream_collide to perform one LBM time step
        self.enqueue_stream_collide().unwrap();
        if self.cfg.graphics_config.graphics_active {
            self.node_communicate_rho_u_flags(world);
        }
        self.node_communicate_fi(world);

        if self.cfg.ext_magneto_hydro && self.cfg.induction_range != 0 {
            self.node_communicate_qi(world);
            self.enqueue_update_e_b_dyn().unwrap();
        }

        self.queue.finish().unwrap();
        self.t += 1;
    }

    #[rustfmt::skip]
    pub fn node_communicate_field(&mut self, field: TransferField, bytes_per_cell: usize, world: &mpi::topology::SimpleCommunicator) {
        world.barrier(); // Get all nodes ready for transfer
        let d_x = self.cfg.d_x as usize;
        let d_y = self.cfg.d_y as usize;
        let d_z = self.cfg.d_z as usize;
        let d = world.rank() as usize;
        let (x, y, z) = ((d % (d_x * d_y)) % d_x, (d % (d_x * d_y)) / d_x, d / (d_x * d_y)); // Domain x, y and z coord

        if d_x > 1 { // Communicate x-axis
            self.enqueue_transfer_extract_field(field, 0, bytes_per_cell).unwrap(); // Extract into transfer buffers
            let dxp = ((x + 1) % d_x) + (y + z * d_y) * d_x;       // domain index of domain at x+1
            let dxm = ((x + d_x - 1) % d_x) + (y + z * d_y) * d_x; // domain index of domain at x-1
            point_to_point::send_receive_into(&self.transfer_p_host[..], &world.process_at_rank(dxp as i32), &mut self.transfer_temp_host, &world.process_at_rank(dxm as i32));
            point_to_point::send_receive_into(&self.transfer_m_host[..], &world.process_at_rank(dxm as i32), &mut self.transfer_p_host, &world.process_at_rank(dxp as i32));
            unsafe {std::ptr::swap(&mut self.transfer_m_host as *mut _, &mut self.transfer_temp_host as *mut _);} // Swap transfer buffers without copying them
            self.enqueue_transfer_insert_field(field, 0, bytes_per_cell).unwrap(); // Insert from transfer buffers

            
        }
        if d_y > 1 { // Communicate y-axis
            self.enqueue_transfer_extract_field(field, 1, bytes_per_cell).unwrap(); // Extract into transfer buffers
            let dyp = x + (((y + 1) % d_y) + z * d_y) * d_x;       // domain index of domain at y+1
            let dym = x + (((y + d_y - 1) % d_y) + z * d_y) * d_x; // domain index of domain at y-1
            let field_length = self.get_area(1) * bytes_per_cell;
            world.process_at_rank(dyp as i32).send(&self.transfer_p_host[..]); // Communicate in y-positive direction
            self.transfer_m.write(&world.process_at_rank(dym as i32).receive_vec::<u8>().0).len(field_length).enq().unwrap();
            world.process_at_rank(dym as i32).send(&self.transfer_m_host[..]); // Communicate in y-negative direction
            self.transfer_p.write(&world.process_at_rank(dyp as i32).receive_vec::<u8>().0).len(field_length).enq().unwrap();
            self.node_enqueue_transfer_insert_field(field, 1).unwrap(); // Insert from transfer buffers
        }
        if d_z > 1 { // Communicate z-axis
            self.enqueue_transfer_extract_field(field, 2, bytes_per_cell).unwrap(); // Extract into transfer buffers
            let dzp = x + (y + ((z + 1) % d_z) * d_y) * d_x;       // domain index of domain at z+1
            let dzm = x + (y + ((z + d_z - 1) % d_z) * d_y) * d_x; // domain index of domain at z-1
            let field_length = self.get_area(2) * bytes_per_cell;
            world.process_at_rank(dzp as i32).send(&self.transfer_p_host[..]); // Communicate in y-positive direction
            self.transfer_m.write(&world.process_at_rank(dzm as i32).receive_vec::<u8>().0).len(field_length).enq().unwrap();
            world.process_at_rank(dzm as i32).send(&self.transfer_m_host[..]); // Communicate in y-negative direction
            self.transfer_p.write(&world.process_at_rank(dzp as i32).receive_vec::<u8>().0).len(field_length).enq().unwrap();
            self.node_enqueue_transfer_insert_field(field, 2).unwrap(); // Insert from transfer buffers
        }
    }

    /// Communicate Fi across domain boundaries
    #[rustfmt::skip]
    fn node_communicate_fi(&mut self, world: &mpi::topology::SimpleCommunicator) {
        let bytes_per_cell =
            self.cfg.float_type.size_of() * self.cfg.velocity_set.get_transfers(); // FP type size * transfers.
        self.node_communicate_field(TransferField::Fi, bytes_per_cell, world);
    }
    
    /// Communicate rho, u and flags across domain boundaries (needed for graphics)
    #[rustfmt::skip]
    fn node_communicate_rho_u_flags(&mut self, world: &mpi::topology::SimpleCommunicator) {
        self.node_communicate_field(TransferField::RhoUFlags, 17, world);
    }

    fn node_communicate_qi(&mut self, world: &mpi::topology::SimpleCommunicator) {
        let bytes_per_cell = self.cfg.float_type.size_of() * 1; // FP type size * transfers. The fixed D3Q7 lattice has 1 transfer
        self.node_communicate_field(TransferField::Qi, bytes_per_cell, world);
    }

    /// Insert field.
    /// Node insert field function does not need to write transfer buffer content to device
    fn node_enqueue_transfer_insert_field(&mut self,field: TransferField,direction: u32) -> ocl::Result<()> {
        let kernel = self.transfer_kernels[field as usize][1]
            .as_ref()
            .expect("transfer kernel should be initialized"); // [1] is the insertion kernel
        kernel.set_arg("direction", direction)?;
        kernel.set_arg("time_step", self.t)?;
        unsafe {
            kernel
                .cmd()
                .global_work_size(self.get_area(direction))
                .enq()
        }
    }

}
