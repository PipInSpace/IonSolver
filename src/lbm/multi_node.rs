//! # multi-node
//! IonSolver supports execution on multiple compute nodes using the MPI protocol if compiled with the `multi-node` feature.
//! 
//! If compiled with the multi-node feature, multiple programs may need to be started over MPI for IonSover to work properly (depending on the configuration used). 
//! 
//! ## MPI Error codes
//! - **100**: Incompatible domain number and number of execution nodes

use std::{cmp::max, f32::consts::PI};

use image::{ImageBuffer, Rgb};
use mpi::{datatype::PartitionMut, point_to_point, topology::SimpleCommunicator, traits::*, Count};
use ocl_macros::*;

use super::*;


/// Run IonSolver on multiple compute nodes without multi-gpu support.
/// This function starts a single node with control over one `LbmDomain`.
pub fn run_node() {
    let universe = mpi::initialize().unwrap();
    let world: SimpleCommunicator = universe.world();
    let size = world.size();
    let rank = world.rank();
    rprintln("IonSolver - © 2024\n", &world);
    rprintln("IonSolver - © 2024\n", &world);
    println!("Launched Node {} of {}", rank, size);

    let mut cfg: LbmConfig;
    if rank == 0 { // If root, read and send config to all nodes
        cfg = LbmConfig::new();
        cfg.units.set(128.0, 1.0, 1.0, 1.0, 0.1, 1.0, 1.2250, 0.0000000001);
        cfg.n_x = 128;
        cfg.n_y = 128;
        cfg.n_z = 256;
        cfg.d_z = 2;
        cfg.nu = cfg.units.si_to_nu(1.48E-5);
        cfg.velocity_set = VelocitySet::D3Q19;
        cfg.mhd_lod_depth = 4;
        cfg.run_steps = 1000;
        // Extensions
        cfg.ext_volume_force = true;
        cfg.ext_magneto_hydro = true;
        // Graphics
        cfg.graphics_config.graphics_active = true;
        cfg.graphics_config.background_color = 0x1c1b22;
        cfg.graphics_config.camera_width = 1920;
        cfg.graphics_config.camera_height = 1080;
        cfg.graphics_config.streamline_every = 8;
        cfg.graphics_config.vec_vis_mode = graphics::VecVisMode::EDyn;
        cfg.graphics_config.u_max = 100.0;
        cfg.graphics_config.streamline_mode = true;
        cfg.graphics_config.axes_mode = true;
        // Animation
        cfg.graphics_config.render_intervals = true;
        let mut f: Vec<graphics::Keyframe> = vec![];
        for i in 0..100 {
            f.push(graphics::Keyframe {
                time: i * 5,
                repeat: false,
                cam_rot_x: 50.0 + (i * 5) as f32,
                cam_rot_y: -20.0,
                cam_zoom: 2.0,
                ..graphics::Keyframe::default()
            })
        }
        cfg.graphics_config.keyframes = f;

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
    if rank == 0 {
        let cpc = 0.002;
        let charge: Vec<f32> = vec![cpc; 128 * 128 * 32];
        bwrite!(domain.q.as_ref().expect("q"), charge);
    }
    world.barrier();
    rprintln("Beginning execution\n", &world);
    //domain.node_taylor_green(1.0, &world);
    domain.node_setup_velocity_field((0.1, 0.01, 0.0), 1.0);
    domain.node_initialize(&world);

    // Set correct camera parameters
    let mut params = graphics::camera_params_rot(
        50.0 * (PI / 180.0),
        -40.0 * (PI / 180.0),
    );
    params[0] = 2.0;
    bwrite!(domain.graphics.as_ref().expect("graphics").camera_params, params);

    // Run for the configured step amount
    for s in 0..cfg.run_steps {
        // Run simulation time step
        domain.node_do_time_step(&world);
        // Render animation frames if needed
        if cfg.graphics_config.graphics_active {
            domain.node_render_keyframes(s, &world);
        }
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

    LbmDomain::new( &cfg, default_device!(), x, y, z, rank) // Initialize and return domain
}

// Additional functionality needed for multi-node exectution
impl LbmDomain {
    /// Readies the LBM Simulation to be run.
    /// Executes `kernel_initialize` Kernels for the domain and communicates between nodes
    pub fn node_initialize(&mut self, world: &SimpleCommunicator) {
        // the communicate calls at initialization need an odd time step
        self.t = 1;
        self.node_communicate_rho_u_flags(world);
        self.enqueue_initialize().unwrap();
        self.node_communicate_rho_u_flags(world);
        self.node_communicate_fi(world);
        if self.cfg.ext_magneto_hydro {
            self.node_communicate_qi(world);
            self.node_communicate_qu_lods(world);
            self.enqueue_update_e_b_dyn().unwrap();
        }
        self.queue.finish().unwrap();
        self.t = 0;
    }

    /// Executes one LBM time step.
    /// Executes `kernel_stream_collide` Kernels for the node `LbmDomain` and exchanges domain transfer buffers.
    /// Updates the dynamic E and B fields.
    fn node_do_time_step(&mut self, world: &SimpleCommunicator) {
        if self.cfg.ext_magneto_hydro {
            self.enqueue_clear_qu_lod().unwrap(); // Ready LOD Buffer
        }
        // call kernel stream_collide to perform one LBM time step
        self.enqueue_stream_collide().unwrap();
        self.queue.finish().unwrap();
        if self.cfg.graphics_config.graphics_active {
            self.node_communicate_rho_u_flags(world);
        }
        self.node_communicate_fi(world);

        if self.cfg.ext_magneto_hydro {
            self.node_communicate_qi(world);
            self.node_communicate_qu_lods(world);
            self.enqueue_update_e_b_dyn().unwrap();
        }

        self.queue.finish().unwrap();
        self.t += 1;
    }

    // Domain communication
    /// Communicate a field across domain barriers
    #[rustfmt::skip]
    pub fn node_communicate_field(&mut self, field: TransferField, bytes_per_cell: usize, world: &SimpleCommunicator) {
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
            point_to_point::send_receive_into(&self.transfer_p_host[..], &world.process_at_rank(dxp as i32), &mut self.transfer_t_host, &world.process_at_rank(dxm as i32)); // Communicate in x-positive direction
            point_to_point::send_receive_into(&self.transfer_m_host[..], &world.process_at_rank(dxm as i32), &mut self.transfer_p_host, &world.process_at_rank(dxp as i32)); // Communicate in x-negative direction
            unsafe {std::ptr::swap(&mut self.transfer_m_host as *mut _, &mut self.transfer_t_host as *mut _);} // Swap transfer buffers without copying them
            self.enqueue_transfer_insert_field(field, 0, bytes_per_cell).unwrap(); // Insert from transfer buffers
        }
        if d_y > 1 { // Communicate y-axis
            self.enqueue_transfer_extract_field(field, 1, bytes_per_cell).unwrap(); // Extract into transfer buffers
            let dyp = x + (((y + 1) % d_y) + z * d_y) * d_x;       // domain index of domain at y+1
            let dym = x + (((y + d_y - 1) % d_y) + z * d_y) * d_x; // domain index of domain at y-1
            point_to_point::send_receive_into(&self.transfer_p_host[..], &world.process_at_rank(dyp as i32), &mut self.transfer_t_host, &world.process_at_rank(dym as i32)); // Communicate in y-positive direction
            point_to_point::send_receive_into(&self.transfer_m_host[..], &world.process_at_rank(dym as i32), &mut self.transfer_p_host, &world.process_at_rank(dyp as i32)); // Communicate in y-negative direction
            unsafe {std::ptr::swap(&mut self.transfer_m_host as *mut _, &mut self.transfer_t_host as *mut _);} // Swap transfer buffers without copying them
            self.enqueue_transfer_insert_field(field, 1, bytes_per_cell).unwrap(); // Insert from transfer buffers
        }
        if d_z > 1 { // Communicate z-axis
            self.enqueue_transfer_extract_field(field, 2, bytes_per_cell).unwrap(); // Extract into transfer buffers
            let dzp = x + (y + ((z + 1) % d_z) * d_y) * d_x;       // domain index of domain at z+1
            let dzm: usize = x + (y + ((z + d_z - 1) % d_z) * d_y) * d_x; // domain index of domain at z-1
            point_to_point::send_receive_into(&self.transfer_p_host[..], &world.process_at_rank(dzp as i32), &mut self.transfer_t_host, &world.process_at_rank(dzm as i32)); // Communicate in z-positive direction
            point_to_point::send_receive_into(&self.transfer_m_host[..], &world.process_at_rank(dzm as i32), &mut self.transfer_p_host, &world.process_at_rank(dzp as i32)); // Communicate in z-negative direction
            unsafe {std::ptr::swap(&mut self.transfer_m_host as *mut _, &mut self.transfer_t_host as *mut _);} // Swap transfer buffers without copying them
            self.enqueue_transfer_insert_field(field, 2, bytes_per_cell).unwrap(); // Insert from transfer buffers
        }
    }

    /// Communicate Fi across domain boundaries
    #[rustfmt::skip]
    fn node_communicate_fi(&mut self, world: &SimpleCommunicator) {
        let bytes_per_cell =
            self.cfg.float_type.size_of() * self.cfg.velocity_set.get_transfers(); // FP type size * transfers.
        self.node_communicate_field(TransferField::Fi, bytes_per_cell, world);
    }
    
    /// Communicate rho, u and flags across domain boundaries (needed for graphics)
    #[rustfmt::skip]
    fn node_communicate_rho_u_flags(&mut self, world: &SimpleCommunicator) {
        self.node_communicate_field(TransferField::RhoUFlags, 17, world);
    }

    /// Communicate Qi across domain boundaries (needed for magnetohydrodynamics)
    fn node_communicate_qi(&mut self, world: &SimpleCommunicator) {
        let bytes_per_cell = self.cfg.float_type.size_of() * 1; // FP type size * transfers. The fixed D3Q7 lattice has 1 transfer
        self.node_communicate_field(TransferField::Qi, bytes_per_cell, world);
    }

    fn node_communicate_qu_lods(&mut self, world: &SimpleCommunicator) {
        let d_n = world.size(); // Number of nodes/domains
        let dim = self.cfg.velocity_set.get_set_values().0 as u32;
        let d = world.rank();

        fn get_offset(depth: i32, dim: u32) -> usize {
            let mut c = 0;
            for i in 0..=depth {c += ((1<<i) as usize).pow(dim)};
            c
        }

        if d_n > 1 {
            self.read_lods();
            let (x, y, z) = get_coordinates_sl(d as u64, self.cfg.d_x, self.cfg.d_y); // Own domain coordinate

            for dc in 0..d_n { //### Every node loops over every node...
                let gather_node = world.process_at_rank(dc); // Node that gathers data in this cycle.

                if dc == d { //### ...if it reaches itself it gathers data from all other nodes and writes to it's own LOD buffer... 
                    // Construct data partition:
                    let mut counts: Vec<Count> = vec![0; d_n as usize];
                    let mut displs: Vec<Count> = vec![0; d_n as usize];
                    let mut offset = 0;
                    for di in 0..d_n {
                        if di == d {continue;} // gather process sends no data to itself
                        let (dx, dy, dz) = get_coordinates_sl(di as u64, self.cfg.d_x, self.cfg.d_y); // Foreign node coordinate
                        let dist: i32 = max((z as i32 - dz as i32).abs(), max((y as i32 - dy as i32).abs(), (x as i32 - dx as i32).abs()));
                        let depth = max(0, self.cfg.mhd_lod_depth as i32 - dist);
                        let data_length = get_offset(depth, dim) as i32 - get_offset(depth - 1, dim) as i32;
                        counts[di as usize] = data_length * 4;
                        displs[di as usize] = offset;
                        offset += data_length * 4;
                    }
                    // Gather partitioned data from all other nodes
                    let mut partition = PartitionMut::new(&mut self.transfer_lod_host.as_mut().expect("msg")[self.n_lod_own*4..], counts, displs);
                    gather_node.gather_varcount_into_root(&vec![0.0f32; 0], &mut partition);
                } else { //### ...otherwise it simply sends it's own relevant LOD data to the current gathering node.
                    let (dx, dy, dz) = get_coordinates_sl(dc as u64, self.cfg.d_x, self.cfg.d_y); // Gather node coordinate
                    let dist: i32 = max((z as i32 - dz as i32).abs(), max((y as i32 - dy as i32).abs(), (x as i32 - dx as i32).abs()));
                    let depth = max(0, self.cfg.mhd_lod_depth as i32 - dist);
                    // This is the range of relevant LOD data for the current foreign domain
                    let range_s = get_offset(depth - 1, dim); // Range start
                    let range_e = get_offset(depth, dim); // Range end
                    gather_node.gather_varcount_into(&self.transfer_lod_host.as_ref().expect("msg")[range_s*4..range_e*4])
                }
            }
            self.qu_lod.as_ref().expect("msg").write(&self.transfer_lod_host.as_ref().expect("msg")[self.n_lod_own*4..]).offset(self.n_lod_own*4).enq().unwrap();
        }
    }

    // Multi-node graphics
    /// Draw one frame of the simulation.
    /// Each node/domain renders its own frame. Frames are transmitted to the root node, stitched together and saved.
    #[rustfmt::skip]
    fn node_draw_frame(&self, name: String, step: u64, world: &SimpleCommunicator) {
        let width = self.cfg.graphics_config.camera_width;
        let height = self.cfg.graphics_config.camera_height;
        let pixel_n = (width * height) as usize;
        let mut bitmap: Vec<i32> = vec![0; pixel_n]; // Base bitmap
        let mut zbuffer: Vec<i32> = vec![0; pixel_n];
        let d_n = world.size();

        self.enqueue_draw_frame(); // Draw domain frame and read into host buffers
        self.queue.finish().unwrap();
        bread!(self.graphics.as_ref().expect("graphics").bitmap, bitmap);
        bread!(self.graphics.as_ref().expect("graphics").zbuffer, zbuffer);

        if world.rank() != 0 { // Transmit buffers to root
            world.process_at_rank(0).gather_into(&bitmap);
            world.process_at_rank(0).gather_into(&zbuffer);
        } else { // Root: Receive buffers from all nodes and assemble images
            let mut bitmaps = vec![0i32; pixel_n * d_n as usize]; // Master bitmap containing pixels of all domains
            let mut zbuffers = vec![0i32; pixel_n * d_n as usize];
            world.process_at_rank(0).gather_into_root(&bitmap, &mut bitmaps[..]); // Receive buffers
            world.process_at_rank(0).gather_into_root(&zbuffer, &mut zbuffers[..]);

            std::thread::spawn(move || { // Generating images needs own thread for performance reasons
                for d in 1..d_n as usize { // Assemble image from domain pixels
                    let offset = pixel_n * d;
                    for i in 0..pixel_n {
                        let zdi = zbuffers[offset + i];
                        if zdi > zbuffer[i] {
                            bitmap[i] = bitmaps[offset + i];
                            zbuffer[i] = zbuffers[offset + i];
                        }
                    }
                }

                let mut save_buffer: Vec<u8> = Vec::with_capacity(bitmap.len());
                for pixel in &bitmap {
                    let color = pixel & 0xFFFFFF;
                    save_buffer.push(((color >> 16) & 0xFF) as u8);
                    save_buffer.push(((color >> 8) & 0xFF) as u8);
                    save_buffer.push((color & 0xFF) as u8);
                }
                let imgbuffer: ImageBuffer<Rgb<u8>, _> = ImageBuffer::from_raw(width, height, save_buffer).unwrap();
                imgbuffer.save(format!(r"{}/../out/{}_{}.png", std::env::current_exe().unwrap().display(), name, step)).unwrap();
                
            });
        }
    }

    fn node_render_keyframes(&self, step: u64, world: &SimpleCommunicator) {
        if self.cfg.graphics_config.render_intervals {
            for (c, frame) in self.cfg.graphics_config.keyframes.iter().enumerate() {
                if (frame.time % (step+1) == 0 && frame.repeat) || (frame.time == step && !frame.repeat) {
                    let mut params = graphics::camera_params_rot(
                        frame.cam_rot_x * (PI / 180.0),
                        frame.cam_rot_y * (PI / 180.0),
                    );
                    params[0] = frame.cam_zoom;
                    bwrite!(self.graphics.as_ref().expect("graphics").camera_params, params);
                    self.node_draw_frame(format!("frame_{}", c), step, world)
                }
            }
        }
    }

    #[rustfmt::skip]
    #[allow(dead_code)]
    fn node_taylor_green(&mut self, periodicity: f32, world: &SimpleCommunicator) {
        let nx = self.cfg.n_x;
        let ny = self.cfg.n_y;
        let nz = self.cfg.n_z;
        let dx = self.cfg.d_x;
        let dy = self.cfg.d_y;
        let dz = self.cfg.d_z;
        let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64 * 2; // Domain size on each axis
        let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64 * 2; // Needs to account for halo offsets
        let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64 * 2;
        let dst = dsx * dsy * dsz;
        let pif = std::f32::consts::PI;
        let A = 0.25f32;
        let a = nx as f32 / periodicity as f32;
        let b = ny as f32 / periodicity as f32;
        let c = nz as f32 / periodicity as f32;

        let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
        let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
        let d = world.rank() as u32;
        let x = (d % (dx * dy)) % dx; // Current Domain coordinates
        let y = (d % (dx * dy)) / dx;
        let z = d / (dx * dy);
        for zi in 0..dsz { // iterates over every cell in the domain, filling it with the velocity field
            for yi in 0..dsy {
                for xi in 0..dsx {
                    if !self.is_halo(xi as u32, yi as u32, zi as u32){
                        //do not set at halo offsets
                        let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                        let gx = xi - (dx > 1u32) as u64 + x as u64 * (dsx - (dx > 1u32) as u64 * 2); // Global coordinates
                        let gy = yi - (dy > 1u32) as u64 + y as u64 * (dsy - (dy > 1u32) as u64 * 2);
                        let gz = zi - (dz > 1u32) as u64 + z as u64 * (dsz - (dz > 1u32) as u64 * 2);
                        let fx = gx as f32 + 0.5 - 0.5 * nx as f32;
                        let fy = gy as f32 + 0.5 - 0.5 * ny as f32;
                        let fz = gz as f32 + 0.5 - 0.5 * nz as f32;
                        domain_vec_u[(dn) as usize] = A
                            * (2.0 * pif * fx / a).cos()
                            * (2.0 * pif * fy / b).sin()
                            * (2.0 * pif * fz / c).sin(); // x
                        domain_vec_u[(dn + dst) as usize] = -A
                            * (2.0 * pif * fx / a).sin()
                            * (2.0 * pif * fy / b).cos()
                            * (2.0 * pif * fz / c).sin(); // y;
                        domain_vec_u[(dn + dst * 2) as usize] = A
                            * (2.0 * pif * fx / a).sin()
                            * (2.0 * pif * fy / b).sin()
                            * (2.0 * pif * fz / c).cos(); // z
                        domain_vec_rho[(dn) as usize] = 1.0
                            - (A * A) * 3.0 / 4.0 * (4.0 * pif * fx / a).cos()
                            + (4.0 * pif * fy / b).cos();
                    }
                }
            }
        }
        
        // Write to domain buffers
        bwrite!(self.u, domain_vec_u);
        bwrite!(self.rho, domain_vec_rho);
        self.queue.finish().unwrap();
    }

    #[rustfmt::skip]
    #[allow(dead_code)]
    fn node_setup_velocity_field(&mut self, velocity: (f32, f32, f32), density: f32) {
        println!("Setting up velocity field");
        let nx = self.cfg.n_x;
        let ny = self.cfg.n_y;
        let nz = self.cfg.n_z;
        let dx = self.cfg.d_x;
        let dy = self.cfg.d_y;
        let dz = self.cfg.d_z;
        let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64 * 2; // Domain size on each axis
        let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64 * 2; // Needs to account for halo offsets
        let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64 * 2;
        let dtotal = dsx * dsy * dsz;
        let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
        let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
        for zi in 0..dsz {
            for yi in 0..dsy {
                for xi in 0..dsx {
                    if !(((xi == 0 || xi == dsx - 1) && dx > 1)
                        || ((yi == 0 || yi == dsy - 1) && dy > 1)
                        || ((zi == 0 || zi == dsz - 1) && dz > 1))
                    {
                        //do not set at halo offsets
                        let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                        domain_vec_u[(dn) as usize] = velocity.0; // x
                        domain_vec_u[(dn + dtotal) as usize] = velocity.1; // y;
                        domain_vec_u[(dn + dtotal * 2) as usize] = velocity.2; // z
                        domain_vec_rho[(dn) as usize] = density;
                    }
                }
            }
        }

        // Write to domain buffers
        bwrite!(self.u, domain_vec_u);
        bwrite!(self.rho, domain_vec_rho);
        self.queue.finish().unwrap();
    }
}

/// Print only if root process. Needs access to comm-world 
fn rprintln(str: &str, world: &SimpleCommunicator) {
    if world.rank() == 0 {
        println!("{str}");
    }
}