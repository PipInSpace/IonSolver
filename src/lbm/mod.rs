//! # lbm
//! 
//! Contains methods for creating and running magnetohydrodynamic LBM simulations
//! 
//! Structured into:
//! - `lbm`: LBM configuration structs, LBM simulation master struct, LbmDomain structs for single devices
//! - `mod graphics`: Graphics functionality for visualization
//! - `mod multi-node`: Extensions for simulation execution on multiple compute nodes (optional)
//! - `mod precompute`: Precomputation of electric and magnetic fields for simulation start up
//! - `mod units`: Units struct for unit conversion between simulation and SI units
//! 
//! This is the core funtionality needed for IonSolver.

use ocl::{Buffer, Context, Device, Kernel, Platform, Program, Queue};
use std::cmp;

pub mod graphics;
#[cfg(feature = "multi-node")]
pub mod multi_node;
pub mod precompute;
pub mod units;
mod types;

use crate::*;
use crate::lbm::{graphics::GraphicsConfig, units::Units};
pub use crate::lbm::types::*;

// Helper Functions
/// Get `x, y, z` coordinates from 1D index `n` and side lengths `n_x` and `n_y`.
pub fn get_coordinates_sl(n: u64, n_x: u32, n_y: u32) -> (u32, u32, u32) {
    let t: u64 = n % (n_x as u64 * n_y as u64);
    //x, y, z
    (
        (t % n_x as u64) as u32,
        (t / n_x as u64) as u32,
        (n / (n_x as u64 * n_y as u64)) as u32,
    )
}

/// Struct used to bundle arguments for LBM simulation setup.
/// 
/// Use this struct to configure a simulation and then instantiate the simulation with `Lbm::new(LbmConfig)`.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct LbmConfig {
    // Holds init information about the simulation
    /// Velocity discretization mode
    pub velocity_set: VelocitySet,
    /// Simulation relaxation time type
    pub relaxation_time: RelaxationTime,
    /// Simulation float type for memory compression of simulation data (DDFs)
    pub float_type: FloatType,
    /// Struct for unit conversion
    pub units: Units,
    /// Size of simulation on each axis
    pub n_x: u32,
    pub n_y: u32,
    pub n_z: u32,

    /// Number of domains on each axis
    pub d_x: u32,
    pub d_y: u32,
    pub d_z: u32,

    /// Kinematic viscosity
    pub nu: f32,
    /// Force on each axis. Needs ext_volume_force to work
    pub f_x: f32,
    pub f_y: f32,
    pub f_z: f32,

    //Extensions
    /// Enable equilibrium boudaries extension (for inlets/outlets)
    pub ext_equilibrium_boudaries: bool,
    /// Enable volume force extension (for handling of complex forces)
    pub ext_volume_force: bool,
    /// Enable force field extension (i. E. for gravity). Needs ext_volume_force to work
    pub ext_force_field: bool,
    /// Enable magnetohydrodynamics extension. Needs ext_volume_force to work 
    pub ext_magneto_hydro: bool,

    /// LOD option for dynamic fields.
    /// 
    /// Dynamic field quality improves with higher values: 1 is very coarse, 4 is higher quality (Performance does not peek at lowest values).
    /// Set this to 0 to disable LODs (VERY SLOW). 
    pub mhd_lod_depth: u8,

    /// Configuration struct for the built-in graphics engine
    pub graphics_config: GraphicsConfig,

    /// Run the simulation for x steps
    pub run_steps: u64,
}

impl LbmConfig {
    /// Returns `LbmConfig` with default values
    pub fn new() -> LbmConfig {
        LbmConfig {
            velocity_set: VelocitySet::D2Q9,
            relaxation_time: RelaxationTime::Srt,
            float_type: FloatType::FP16S,
            units: units::Units::new(),
            n_x: 1,
            n_y: 1,
            n_z: 1,
            d_x: 1,
            d_y: 1,
            d_z: 1,
            nu: 1.0f32 / 6.0f32,
            f_x: 0.0f32,
            f_y: 0.0f32,
            f_z: 0.0f32,

            ext_equilibrium_boudaries: false,
            ext_volume_force: false,
            ext_magneto_hydro: false,
            ext_force_field: false,

            mhd_lod_depth: 4, // Dynamic field LODs

            graphics_config: graphics::GraphicsConfig::new(),
            run_steps: 0,
        }
    }
}

/// To start a LBM simulation, initialise an Lbm struct:
/// ```
/// Lbm::new(lbm_config: LbmConfig)
/// ```
/// The `new()` function takes in another struct, the [`LbmConfig`], which contains all necessary arguments.
/// [`LbmConfig`] needs to be configured beforehand.
///
/// The Lbm struct contains one or multiple [`LbmDomain`] structs.
/// The Simulation is actually run on these [`LbmDomain`] structs, each Domain corresponding to an OpenCL device,
/// enabling multi-device parallelization.
/// [`LbmDomain`] initialisation is handled automatically in `Lbm::new()` using `LbmDomain::new()`
///
/// [`LbmConfig`]: crate::lbm::LbmConfig
/// [`LbmDomain`]: crate::lbm::LbmDomain
pub struct Lbm {
    /// Vector of [`LbmDomain`]s that are part of this simulation.
    pub domains: Vec<LbmDomain>,
    /// A copy of the [`LbmConfig`] used to initialize the simulation.
    pub config: LbmConfig,
    /// A vector of Charges positioned in the simulation via a 1D index. Used in static field computation. 
    pub charges: Option<Vec<(u64, f32)>>,
    /// A vector of Magnets positioned in the simulation via a 1D index. Used in static field computation. 
    pub magnets: Option<Vec<(u64, [f32; 3])>>,
    initialized: bool,
}

impl Lbm {
    // TODO: add read for charges and magnets
    /// Returns new `Lbm` struct from pre-configured `LbmConfig` struct. `LbmDomain` setup is handled automatically.
    /// Configures Domains
    pub fn new(mut lbm_config: LbmConfig) -> Lbm {
        let n_d_x: u32 = (lbm_config.n_x / lbm_config.d_x) * lbm_config.d_x;
        let n_d_y: u32 = (lbm_config.n_y / lbm_config.d_y) * lbm_config.d_y;
        let n_d_z: u32 = (lbm_config.n_z / lbm_config.d_z) * lbm_config.d_z;
        if n_d_x != lbm_config.n_x || n_d_y != lbm_config.n_y || n_d_z != lbm_config.n_z {
            println!(
                "Warning: Resolution {}, {}, {} not divisible by Domains: Overiding resolution.",
                lbm_config.n_x, lbm_config.n_y, lbm_config.n_z
            )
        }

        lbm_config.n_x = n_d_x;
        lbm_config.n_y = n_d_y;
        lbm_config.n_z = n_d_z;

        let domain_numbers: u32 = lbm_config.d_x * lbm_config.d_y * lbm_config.d_z;

        let device_infos = opencl::device_selection(domain_numbers);
        //TODO: sanity check

        let mut lbm_domains: Vec<LbmDomain> = Vec::new();
        for d in 0..domain_numbers {
            println!("Initializing domain {}/{}", d + 1, domain_numbers);
            let x = (d % (lbm_config.d_x * lbm_config.d_y)) % lbm_config.d_x;
            let y = (d % (lbm_config.d_x * lbm_config.d_y)) / lbm_config.d_x;
            let z = d / (lbm_config.d_x * lbm_config.d_y);
            println!(
                "    Using \"{}\" for domain {}",
                device_infos[d as usize].name().unwrap(),
                d + 1
            );
            lbm_domains.push(LbmDomain::new(
                &lbm_config,
                device_infos[d as usize],
                x,
                y,
                z,
                d,
            ))
        }
        println!("All domains initialized.\n");

        Lbm {
            domains: lbm_domains,
            config: lbm_config,
            charges: None,
            magnets: None,
            initialized: false,
        }
    }

    /// Readies the LBM Simulation to be run.
    /// Executes `kernel_initialize` Kernels for every `LbmDomain` and fills domain transfer buffers.
    pub fn initialize(&mut self) {
        // the communicate calls at initialization need an odd time step
        self.increment_timestep(1);
        self.communicate_rho_u_flags();
        self.kernel_initialize();
        self.communicate_rho_u_flags();
        self.communicate_fi();
        if self.config.ext_magneto_hydro {
            self.communicate_qi();
            self.communicate_qu_lods();
            self.update_e_b_dynamic();
        }

        self.finish_queues();
        self.reset_timestep();
        self.initialized = true;
    }

    /// Runs Simulations for `steps`
    #[allow(unused)]
    pub fn run(&mut self, steps: u64) {
        //Initialize, then run simulation for steps
        if !self.initialized {
            //Run initialization Kernel
            self.initialize();
        }
        for i in 0..steps {
            //println!("Step {}", i);
            self.do_time_step();
        }
    }

    /// Executes one LBM time step.
    /// Executes `kernel_stream_collide` Kernels for every `LbmDomain` and updates domain transfer buffers.
    /// Updates the dynamic E and B fields.
    pub fn do_time_step(&mut self) {
        if self.config.ext_magneto_hydro {
            self.clear_qu_lod(); // Ready LOD Buffer
        }
        // call kernel stream_collide to perform one LBM time step
        self.stream_collide();
        if self.config.graphics_config.graphics_active {
            self.communicate_rho_u_flags();
        }
        self.communicate_fi();
        if self.config.ext_magneto_hydro {
            self.communicate_qi();
            //let lod = bget!(self.domains[0].qu_lod.as_ref().expect("msg"));
            //println!("1b: {:?}", lod);
            //let lod = bget!(self.domains[1].qu_lod.as_ref().expect("msg"));
            //println!("2b: {:?}", lod);
            self.communicate_qu_lods();
            //let lod = bget!(self.domains[0].qu_lod.as_ref().expect("msg"));
            //println!("1a: {:?}", lod);
            //let lod = bget!(self.domains[1].qu_lod.as_ref().expect("msg"));
            //println!("2a: {:?}", lod);
            self.update_e_b_dynamic();
        }
        if self.get_d_n() == 1 || self.config.ext_magneto_hydro {
            // Additional synchronization only needed in single-GPU or after E&B update
            self.finish_queues();
        }
        self.increment_timestep(1);
    }

    /// Blocks execution until all `LbmDomain` OpenCL queues have finished.
    pub fn finish_queues(&self) {
        for d in 0..self.domains.len() {
            self.domains[d].queue.finish().unwrap();
        }
    }

    /// Executes initialize kernel
    fn kernel_initialize(&self) {
        for d in 0..self.get_d_n() {
            self.domains[d].enqueue_initialize().unwrap();
        }
    }

    /// Execute stream_collide kernel on all domains
    fn stream_collide(&self) {
        for d in 0..self.get_d_n() {
            self.domains[d].enqueue_stream_collide().unwrap();
        }
    }

    /// Execute E and B dynamic update kernel on all domains
    fn update_e_b_dynamic(&self) {
        for d in 0..self.get_d_n() {
            self.domains[d].enqueue_update_e_b_dyn().unwrap();
        }
    }

    /// Reset charge and velocity LOD
    fn clear_qu_lod(&self) {
        for d in 0..self.get_d_n() {
            self.domains[d].enqueue_clear_qu_lod().unwrap();
        }
    }

    // Domain communication
    /// Communicate a field across domain barriers
    #[rustfmt::skip]
    fn communicate_field(&mut self, field: TransferField, bytes_per_cell: usize) {
        let d_x = self.config.d_x as usize;
        let d_y = self.config.d_y as usize;
        let d_z = self.config.d_z as usize;
        let d_n = self.get_d_n();

        if d_x > 1 { // Communicate x-axis
            for d in 0..d_n {self.domains[d].enqueue_transfer_extract_field(field, 0, bytes_per_cell).unwrap();} // Extract into transfer buffers
            self.finish_queues(); // Synchronize domains
            for d in 0..d_n { // Swap transfer buffers at domain boundaries
                let (x, y, z) = ((d % (d_x * d_y)) % d_x, (d % (d_x * d_y)) / d_x, d / (d_x * d_y)); // Domain x, y and z coord
                let dxp = ((x + 1) % d_x) + (y + z * d_y) * d_x; // domain index of domain at x+1
                unsafe {std::ptr::swap(&mut self.domains[d].transfer_p_host as *mut _, &mut self.domains[dxp].transfer_m_host as *mut _);} // Swap transfer buffers without copying them
            }
            for d in 0..d_n {self.domains[d].enqueue_transfer_insert_field(field, 0, bytes_per_cell).unwrap();} // Insert from transfer buffers
        }
        if d_y > 1 { // Communicate y-axis
            for d in 0..d_n {self.domains[d].enqueue_transfer_extract_field(field, 1, bytes_per_cell).unwrap();} // Extract into transfer buffers
            self.finish_queues(); // Synchronize domains
            for d in 0..d_n { // Swap transfer buffers at domain boundaries
                let (x, y, z) = ((d % (d_x * d_y)) % d_x, (d % (d_x * d_y)) / d_x, d / (d_x * d_y)); // Domain x, y and z coord
                let dyp = x + (((y + 1) % d_y) + z * d_y) * d_x; // domain index of domain at y+1
                unsafe {std::ptr::swap(&mut self.domains[d].transfer_p_host as *mut _, &mut self.domains[dyp].transfer_m_host as *mut _);} // Swap transfer buffers without copying them
            }
            for d in 0..d_n {self.domains[d].enqueue_transfer_insert_field(field, 1, bytes_per_cell).unwrap();} // Insert from transfer buffers
        }
        if d_z > 1 { // Communicate z-axis
            for d in 0..d_n {self.domains[d].enqueue_transfer_extract_field(field, 2, bytes_per_cell).unwrap();} // Extract into transfer buffers
            self.finish_queues(); // Synchronize domains
            for d in 0..d_n { // Swap transfer buffers at domain boundaries
                let (x, y, z) = ((d % (d_x * d_y)) % d_x, (d % (d_x * d_y)) / d_x, d / (d_x * d_y)); // Domain x, y and z coord
                let dzp = x + (y + ((z + 1) % d_z) * d_y) * d_x; // domain index of domain at z+1
                unsafe {std::ptr::swap(&mut self.domains[d].transfer_p_host as *mut _, &mut self.domains[dzp].transfer_m_host as *mut _);} // Swap transfer buffers without copying them
            }
            for d in 0..d_n {self.domains[d].enqueue_transfer_insert_field(field, 2, bytes_per_cell).unwrap();} // Insert from transfer buffers
        }
    }

    /// Communicate Fi across domain boundaries
    fn communicate_fi(&mut self) {
        let bytes_per_cell =
            self.config.float_type.size_of() * self.config.velocity_set.get_transfers(); // FP type size * transfers.
        self.communicate_field(TransferField::Fi, bytes_per_cell);
    }
    /// Communicate rho, u and flags across domain boundaries (needed for graphics)
    fn communicate_rho_u_flags(&mut self) {
        self.communicate_field(TransferField::RhoUFlags, 17);
    }

    /// Communicate Qi across domain boundaries
    fn communicate_qi(&mut self) {
        let bytes_per_cell = self.config.float_type.size_of() * 1; // FP type size * transfers. The fixed D3Q7 lattice has 1 transfer
        self.communicate_field(TransferField::Qi, bytes_per_cell);
    }

    /// Communicate charge and velocity LODs across domains
    //#[allow(unused)]
    #[rustfmt::skip]
    fn communicate_qu_lods(&mut self) {
        let d_n = self.get_d_n();
        let dim = self.config.velocity_set.get_set_values().0 as u32;

        fn get_offset(depth: i32, dim: u32) -> usize {
            let mut c = 0;
            for i in 0..=depth {c += ((1<<i) as usize).pow(dim)};
            c
        }

        if d_n > 1 {
            for d in 0..d_n { self.domains[d].read_lods(); }
            for d in 0..d_n { // Every domain...
                let (x, y, z) = get_coordinates_sl(d as u64, self.config.d_x, self.config.d_y); // Own domain coordinate
                let mut offset = self.domains[d].n_lod_own; // buffer write offset
                for dc in 0..d_n { // ...loops over every other domain and writes the extracted LOD to its own LOD buffer
                    if d != dc {
                        let (dx, dy, dz) = get_coordinates_sl(dc as u64, self.config.d_x, self.config.d_y);
                        let dist: i32 = cmp::max((z as i32 - dz as i32).abs(), cmp::max((y as i32 - dy as i32).abs(), (x as i32 - dx as i32).abs()));
                        let depth = cmp::max(0, self.config.mhd_lod_depth as i32 - dist);
                        // This is the range of relevant LOD data for the current foreign domain
                        let range_s = get_offset(depth - 1, dim); // Range start
                        let range_e = get_offset(depth, dim); // Range end
                        // Write to device
                        self.domains[d].qu_lod.as_ref().expect("msg")
                            .write(&self.domains[dc].transfer_lod_host.as_ref().expect("msg")[range_s*4..range_e*4])
                            .offset(offset*4).enq().unwrap();
                        offset += range_e - range_s;
                    }
                }
            }
        }
    }


    // Helper functions
    /// Increments time steps variable for every `LbmDomain`
    fn increment_timestep(&mut self, steps: u32) {
        for d in 0..self.domains.len() {
            self.domains[d].t += steps as u64;
        }
    }

    /// Resets timestep variable for every `LbmDomain`
    fn reset_timestep(&mut self) {
        for d in 0..self.domains.len() {
            self.domains[d].t = 0;
        }
    }

    /// Returns the amount of domains for this lbm
    pub fn get_d_n(&self) -> usize {
        self.domains.len()
    }

    /// Returns the current simulation timestep
    pub fn get_time_step(&self) -> u64 {
        self.domains[0].t
    }
}

/// The `LbmDomain` struct holds all information to run a LBM-Simulation on one OpenCL Device.
/// `LbmDomain` should not be initialized on it's own, but automatically through the `Lbm::new()` function when initializing a new [`Lbm`].
/// This ensures all arguments are correctly set.
/// 
/// [`Lbm`]: crate::lbm::Lbm
pub struct LbmDomain {
    // FluidX3D creates contexts/queues/programs for each device automatically through another class
    pub queue: Queue,

    kernel_initialize: Kernel, // Basic Kernels
    kernel_stream_collide: Kernel,
    kernel_update_fields: Kernel,
    pub transfer_kernels: [[Option<Kernel>; 2]; 3],

    kernel_update_e_b_dyn: Option<Kernel>, // Optional Kernels
    kernel_clear_qu_lod: Option<Kernel>,

    pub n_x: u32, // Domain size
    pub n_y: u32,
    pub n_z: u32,

    fx: f32, // Volume force
    fy: f32,
    fz: f32,

    pub fi: VariableFloatBuffer, // Buffers
    pub rho: Buffer<f32>,
    pub u: Buffer<f32>,
    pub flags: Buffer<u8>,

    pub transfer_p: Buffer<u8>, // Transfer buffers positive/negative
    pub transfer_m: Buffer<u8>, // Size is maximum of (17 bytes (rho, u and flags) or bytes needed for fi)
    pub transfer_p_host: Vec<u8>,
    pub transfer_m_host: Vec<u8>,
    #[cfg(feature = "multi-node")]
    pub transfer_t_host: Vec<u8>,
    pub transfer_lod_host: Option<Vec<f32>>,
    pub n_lod_own: usize, // The number of LOD data chunks that belong to this domain

    pub e: Option<Buffer<f32>>, // Optional Buffers
    pub b: Option<Buffer<f32>>,
    pub e_dyn: Option<Buffer<f32>>,
    pub b_dyn: Option<Buffer<f32>>,
    pub qi: Option<VariableFloatBuffer>,
    pub q: Option<Buffer<f32>>,
    pub qu_lod: Option<Buffer<f32>>,
    pub f: Option<Buffer<f32>>,
    pub t: u64, // Timestep

    pub graphics: Option<graphics::Graphics>, // Graphics struct, handles rendering
    pub cfg: LbmConfig,
}

impl LbmDomain {
    /// Returns new `LbmDomain` from provided arguments.
    ///
    /// `LbmDomain` should not be initialized on it's own, but automatically through the `Lbm::new()` function.
    /// This ensures all arguments are correctly set.
    #[rustfmt::skip]
    pub fn new(lbm_config: &LbmConfig, device: Device, x: u32, y: u32, z: u32, i: u32) -> LbmDomain {
        let n_x = lbm_config.n_x / lbm_config.d_x + 2u32 * (lbm_config.d_x > 1u32) as u32; // Size + Halo offsets
        let n_y = lbm_config.n_y / lbm_config.d_y + 2u32 * (lbm_config.d_y > 1u32) as u32; // When multiple domains on axis -> add 2 cells of padding
        let n_z = lbm_config.n_z / lbm_config.d_z + 2u32 * (lbm_config.d_z > 1u32) as u32;
        let n = Self::get_n(n_x, n_y, n_z);

        let d_x = lbm_config.d_x;
        let d_y = lbm_config.d_y;
        let d_z = lbm_config.d_z;
        let d_n = d_x * d_y * d_z;

        let o_x = (x * lbm_config.n_x / lbm_config.d_x) as i32 - (lbm_config.d_x > 1u32) as i32;
        let o_y = (y * lbm_config.n_y / lbm_config.d_y) as i32 - (lbm_config.d_y > 1u32) as i32;
        let o_z = (z * lbm_config.n_z / lbm_config.d_z) as i32 - (lbm_config.d_z > 1u32) as i32;

        let (dimensions, velocity_set, transfers) = lbm_config.velocity_set.get_set_values();

        let t = 0;

        // MAGNETO_HYDRO LOD buffer size
        let n_lod_own; // Amount of lod blocks of this domain
        let n_lod = {  // Amount of lod blocks of this domain + lod blocks from neighbouring domains
            let mut c = 1;
            // Own domain LODs
            for i in 0..lbm_config.mhd_lod_depth {c += ((1<<(i+1)) as usize).pow(dimensions as u32)}
            n_lod_own = c;
            // Other domain LODs
            for d in 0..d_n {
                let (dx, dy, dz) = get_coordinates_sl(d as u64, d_x, d_y);
                // single-axis manhattan distance to other domain, controls lod amount
                let dist = cmp::max((z as i32 - dz as i32).abs(), cmp::max((y as i32 - dy as i32).abs(), (x as i32 - dx as i32).abs()));
                if dist != 0 {
                    c += ((1<<cmp::max(lbm_config.mhd_lod_depth as i32 - dist, 0)) as usize).pow(dimensions as u32)
                }
            }
            c
        };

        let ocl_code = Self::get_device_defines(n_x, n_y, n_z, d_x, d_y, d_z, i, o_x, o_y, o_z, dimensions, velocity_set, transfers, lbm_config.nu, n_lod, n_lod_own, lbm_config)
            + &if lbm_config.graphics_config.graphics_active { graphics::get_graphics_defines(&lbm_config.graphics_config) } else { "".to_string() }
            + &opencl::get_opencl_code(); // Only appends graphics defines if needed

        // OCL variables are directly exposed, due to no other device struct.
        let platform = Platform::default();
        let context = Context::builder().platform(platform).devices(device).build().unwrap();
        let queue = Queue::new(&context, device, None).unwrap();
        println!("    Compiling Program...");
        let program = Program::builder().devices(device).src(&ocl_code).build(&context).unwrap();

        // Allocate Buffers
        let fi: VariableFloatBuffer = match lbm_config.float_type {
            FloatType::FP32 => { VariableFloatBuffer::F32(buffer!(&queue,[n * velocity_set as u64],0.0f32)) } // Float Type F32
            _               => { VariableFloatBuffer::U16(buffer!(&queue,[n * velocity_set as u64],0u16)) }   // Float Type F16S/F16C
        };
        //let rho = opencl::create_buffer(&queue, [n], 1.0f32);
        let rho =  buffer!(&queue, [n], 1.0f32);
        let u =    buffer!(&queue, [n * 3], 0f32);
        let flags = buffer!(&queue, [n], 0u8);

        // VOLUME_FORCE extension
        // Force field buffer as 3D Vectors
        let f: Option<Buffer<f32>> = if lbm_config.ext_force_field { Some(buffer!(&queue, [n * 3], 0f32)) } else { None };
        // MAGNETO_HYDRO extension
        // Electric field buffer as 3D Vectors
        let e: Option<Buffer<f32>> =     if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n * 3], 0f32)) } else { None };
        let e_dyn: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n * 3], 0f32)) } else { None };
        // Magnetic field buffers as 3D Vectors
        let b: Option<Buffer<f32>> =     if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n * 3], 0f32)) } else { None };
        let b_dyn: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n * 3], 0f32)) } else { None };
        // Charge ddfs and charge field buffers
        let qi: Option<VariableFloatBuffer> = if lbm_config.ext_magneto_hydro {
            Some(match lbm_config.float_type {
                FloatType::FP32 => { VariableFloatBuffer::F32(buffer!(&queue, [n * 7], 0.0f32)) } // Float Type F32
                _               => { VariableFloatBuffer::U16(buffer!(&queue, [n * 7], 0u16)) }   // Float Type F16S/F16C
            })
        } else { None };
        let q: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n], 0f32)) } else { None };
        // Level of detail charge and velocity for field updates
        let qu_lod: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n_lod * 4], 0f32)) } else { None };
        #[cfg(feature = "multi-node")]
        let transfer_lod_host: Option<Vec<f32>> = if lbm_config.ext_magneto_hydro { Some(vec![0.0f32; n_lod * 4]) } else { None }; // Transfer buffer needs to be bigger in multi-node mode
        #[cfg(not(feature = "multi-node"))]
        let transfer_lod_host: Option<Vec<f32>> = if lbm_config.ext_magneto_hydro { Some(vec![0.0f32; n_lod_own * 4]) } else { None };

        // Initialize Kernels
        let mut initialize_builder = kernel_builder!(program, queue, "initialize", [n]);
        let mut stream_collide_builder = kernel_builder!(program, queue, "stream_collide", [n]);
        let mut update_fields_builder = kernel_builder!(program, queue, "update_fields", [n]);
        let mut kernel_update_e_b_dyn: Option<Kernel> = None;
        let mut kernel_clear_qu_lod: Option<Kernel> = None;
        match &fi {
            //Initialize kernels. Different Float types need different arguments (Fi-Buffer specifically)
            VariableFloatBuffer::F32(fif32) => { // Float Type F32
                kernel_args!(initialize_builder, ("fi", fif32));
                kernel_args!(stream_collide_builder, ("fi", fif32));
                kernel_args!(update_fields_builder, ("fi", fif32));
            }
            VariableFloatBuffer::U16(fiu16) => { // Float Type F16S/F16C
                kernel_args!(initialize_builder, ("fi", fiu16));
                kernel_args!(stream_collide_builder, ("fi", fiu16));
                kernel_args!(update_fields_builder, ("fi", fiu16));
            }
        }
        kernel_args!(initialize_builder,     ("rho", &rho), ("u", &u), ("flags", &flags));
        kernel_args!(stream_collide_builder, ("rho", &rho), ("u", &u), ("flags", &flags), ("t", t), ("fx", lbm_config.f_x), ("fy", lbm_config.f_y), ("fz", lbm_config.f_z));
        kernel_args!(update_fields_builder,  ("rho", &rho), ("u", &u), ("flags", &flags), ("t", t), ("fx", lbm_config.f_x), ("fy", lbm_config.f_y), ("fz", lbm_config.f_z));

        // Conditional arguments. Place at end of kernel functions
        if lbm_config.ext_force_field { kernel_args!(stream_collide_builder, ("F", f.as_ref().expect("F"))); }
        if lbm_config.ext_magneto_hydro {
            kernel_args!(stream_collide_builder, ("E", e.as_ref().expect("e")), ("B", b.as_ref().expect("b")), ("E_dyn", e_dyn.as_ref().expect("e_dyn")), ("B_dyn", b_dyn.as_ref().expect("b_dyn")));
            kernel_args!(initialize_builder,     ("E", e.as_ref().expect("e")), ("B", b.as_ref().expect("b")), ("E_dyn", e_dyn.as_ref().expect("e_dyn")), ("B_dyn", b_dyn.as_ref().expect("b_dyn")));
            
            match qi.as_ref().expect("qi should be initialized") {
                VariableFloatBuffer::U16(qi_u16) => {kernel_args!(stream_collide_builder, ("qi", qi_u16)); kernel_args!(initialize_builder, ("qi", qi_u16));},
                VariableFloatBuffer::F32(qi_f32) => {kernel_args!(stream_collide_builder, ("qi", qi_f32)); kernel_args!(initialize_builder, ("qi", qi_f32));},
            };

            kernel_args!(initialize_builder,     ("Q", q.as_ref().expect("Q")));
            kernel_args!(stream_collide_builder, ("Q", q.as_ref().expect("Q")), ("QU_lod", qu_lod.as_ref().expect("QU_lod")));

            // Dynamic E/B kernel
            kernel_update_e_b_dyn = Some(
                kernel!(program, queue, "update_e_b_dynamic", [n], ("E_dyn", e_dyn.as_ref().expect("e_dyn")), ("B_dyn", b_dyn.as_ref().expect("b_dyn")), ("Q", q.as_ref().expect("q")), ("u", &u), ("QU_lod", qu_lod.as_ref().expect("QU_lod")), ("flags", &flags))
            );
            // Clear LOD
            kernel_clear_qu_lod = Some(
                kernel!(program, queue, "clear_qu_lod", [n_lod], ("QU_lod", qu_lod.as_ref().expect("QU_lod")))
            );
        }
        
        let kernel_stream_collide: Kernel = stream_collide_builder.build().unwrap();
        let kernel_initialize: Kernel = initialize_builder.build().unwrap();
        let kernel_update_fields: Kernel = update_fields_builder.build().unwrap();


        // Multi-Domain-Transfers:
        // Transfer buffer initializaton:
        let mut a_max: usize = 0;
        if lbm_config.d_x > 1 { a_max = cmp::max(a_max, n_y as usize * n_z as usize); } // Ax
        if lbm_config.d_y > 1 { a_max = cmp::max(a_max, n_x as usize * n_z as usize); } // Ay
        if lbm_config.d_z > 1 { a_max = cmp::max(a_max, n_x as usize * n_y as usize); } // Az
        let transfer_buffer_size = a_max * cmp::max(17, transfers as usize * lbm_config.float_type.size_of());

        let transfer_p = buffer!(&queue, transfer_buffer_size, 0_u8); // Size is maxiumum 17 bytes or bytes needed for fi
        let transfer_m = buffer!(&queue, transfer_buffer_size, 0_u8); // Size is maxiumum 17 bytes or bytes needed for fi
        let transfer_p_host: Vec<u8> = vec![0_u8; transfer_buffer_size];
        let transfer_m_host: Vec<u8> = vec![0_u8; transfer_buffer_size];
        #[cfg(feature = "multi-node")]
        let transfer_t_host: Vec<u8> = vec![0_u8; transfer_buffer_size]; // Temporary buffer for node communication

        // Transfer kernels
        let transfer_kernels = [
            [
                Some(match &fi {
                    VariableFloatBuffer::U16(fi_u16) => kernel!(program, queue, "transfer_extract_fi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fi", fi_u16)),
                    VariableFloatBuffer::F32(fi_f32) => kernel!(program, queue, "transfer_extract_fi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fi", fi_f32)),
                }), // Extract Fi
                Some(match &fi {
                    VariableFloatBuffer::U16(fi_u16) => kernel!(program, queue, "transfer__insert_fi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fi", fi_u16)),
                    VariableFloatBuffer::F32(fi_f32) => kernel!(program, queue, "transfer__insert_fi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fi", fi_f32)),
                }), // Insert Fi
            ], // Fi
            [
                Some(kernel!(program, queue, "transfer_extract_rho_u_flags", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("rho", &rho), ("u", &u), ("flags", &flags))), // Extract rho, u, flags
                Some(kernel!(program, queue, "transfer__insert_rho_u_flags", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("rho", &rho), ("u", &u), ("flags", &flags))), // Insert  rho, u, flags
            ], // Rho, u and flags (needed for graphics)
            [
                if lbm_config.ext_magneto_hydro {
                    Some(match &qi.as_ref().expect("qi should be defined") {
                        VariableFloatBuffer::U16(qi_u16) => kernel!(program, queue, "transfer_extract_qi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("qi", qi_u16)),
                        VariableFloatBuffer::F32(qi_f32) => kernel!(program, queue, "transfer_extract_qi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("qi", qi_f32)),
                    })
                } else { None } , // Extract Qi
                if lbm_config.ext_magneto_hydro {
                    Some(match &qi.as_ref().expect("qi should be defined") {
                        VariableFloatBuffer::U16(qi_u16) => kernel!(program, queue, "transfer__insert_qi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("qi", qi_u16)),
                        VariableFloatBuffer::F32(qi_f32) => kernel!(program, queue, "transfer__insert_qi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("qi", qi_f32)),
                    })
                } else { None }, // Insert Qi
            ], // Qi charge ddfs
        ];
        println!("    Kernels for domain compiled.");


        let graphics: Option<graphics::Graphics> = if lbm_config.graphics_config.graphics_active { Some(graphics::Graphics::new(lbm_config, &program, &queue, &flags, &u, (n_x, n_y, n_z))) } else { None };

        LbmDomain {
            queue,

            kernel_initialize,
            kernel_stream_collide,
            kernel_update_fields,
            transfer_kernels,

            kernel_update_e_b_dyn,
            kernel_clear_qu_lod,

            n_x, n_y, n_z,
            fx: lbm_config.f_x, fy: lbm_config.f_y, fz: lbm_config.f_z,
            fi, rho, u, flags,

            transfer_p, transfer_m,
            transfer_p_host, transfer_m_host,
            #[cfg(feature = "multi-node")]
            transfer_t_host,
            transfer_lod_host,
            n_lod_own,

            e, b, e_dyn, b_dyn, qi, q, qu_lod, f, t, graphics, cfg: lbm_config.clone()
        } //Returns initialised domain
    }

    /// Returns a string of OpenCL C `#define`s from the provided arguments that are appended to the base OpenCl code at runtime.
    /// These are unique for every domain.
    fn get_device_defines(
        n_x: u32,
        n_y: u32,
        n_z: u32,
        d_x: u32,
        d_y: u32,
        d_z: u32,
        d_i: u32,
        o_x: i32,
        o_y: i32,
        o_z: i32,
        dimensions: u8,
        velocity_set: u8,
        transfers: u8,
        nu: f32,
        n_lod: usize,
        n_lod_own: usize,
        lbm_config: &LbmConfig,
    ) -> String {
        //Conditional Defines:
        //Velocity set types
        let d2q9 = "\n	#define DEF_W0 (1.0f/2.25f)".to_owned() // center (0)
        +"\n	#define DEF_WS (1.0f/9.0f)" // straight (1-4)
        +"\n	#define DEF_WE (1.0f/36.0f)";
        let d3q15 = "\n	#define DEF_W0 (1.0f/4.5f)".to_owned() // center (0)
        +"\n	#define DEF_WS (1.0f/9.0f)" // straight (1-6)
        +"\n	#define DEF_WC (1.0f/72.0f)";
        let d3q19 = "\n	#define DEF_W0 (1.0f/3.0f)".to_owned() // center (0)
        +"\n	#define DEF_WS (1.0f/18.0f)" // straight (1-6)
        +"\n	#define DEF_WE (1.0f/36.0f)";
        let d3q27 = "\n	#define DEF_W0 (1.0f/3.0f)".to_owned() // center (0)
        +"\n	#define DEF_WS (1.0f/18.0f)" // straight (1-6)
        +"\n	#define DEF_WE (1.0f/36.0f)";
        //Relaxation time types
        let srt = "\n	#define SRT";
        let trt = "\n	#define TRT";
        //Float types
        let fp16s = "\n	#define fpxx half".to_owned() // switchable data type (scaled IEEE-754 16-bit floating-point format: 1-5-10, exp-30, +-1.99902344, +-1.86446416E-9, +-1.81898936E-12, 3.311 digits)
        +"\n	#define fpxx_copy ushort" // switchable data type for direct copying (scaled IEEE-754 16-bit floating-point format: 1-5-10, exp-30, +-1.99902344, +-1.86446416E-9, +-1.81898936E-12, 3.311 digits)
        +"\n	#define load(p,o) vload_half(o,p)*3.0517578E-5f" // special function for loading half
        +"\n	#define store(p,o,x) vstore_half_rte((x)*32768.0f,o,p)"; // special function for storing half
        let fp16c = "\n	#define fpxx ushort".to_owned() // switchable data type (custom 16-bit floating-point format: 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits), 12.5% slower than IEEE-754 16-bit
        +"\n	#define fpxx_copy ushort" // switchable data type for direct copying (custom 16-bit floating-point format: 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits), 12.5% slower than IEEE-754 16-bit
        +"\n	#define load(p,o) half_to_float_custom(p[o])" // special function for loading half
        +"\n	#define store(p,o,x) p[o]=float_to_half_custom(x)"; // special function for storing half
        let fp32 = "\n	#define fpxx float".to_owned() // switchable data type (regular 32-bit float)
        +"\n	#define fpxx_copy float" // switchable data type for direct copying (regular 32-bit float)
        +"\n	#define load(p,o) p[o]" // regular float read
        +"\n	#define store(p,o,x) p[o]=x"; // regular float write

        // Return String
        "\n    #define DEF_NX ".to_owned() + &n_x.to_string()+"u"
        +"\n	#define DEF_NY "+ &n_y.to_string()+"u"
        +"\n	#define DEF_NZ "+ &n_z.to_string()+"u"
        +"\n	#define DEF_N  "+ &Self::get_n(n_x, n_y, n_z).to_string()+"ul"

        +"\n	#define DEF_DX "+ &d_x.to_string()+"u"
        +"\n	#define DEF_DY "+ &d_y.to_string()+"u"
        +"\n	#define DEF_DZ "+ &d_z.to_string()+"u"
        +"\n	#define DEF_DI "+ &d_i.to_string()+"u"

        +"\n	#define DEF_OX "+ &o_x.to_string()+"" // offsets are signed integer!
        +"\n	#define DEF_OY "+ &o_y.to_string()+""
        +"\n	#define DEF_OZ "+ &o_z.to_string()+""

        +"\n	#define DEF_AX "+ &(n_y * n_z).to_string()+"u"
        +"\n	#define DEF_AY "+ &(n_z * n_x).to_string()+"u"
        +"\n	#define DEF_AZ "+ &(n_x * n_y).to_string()+"u"

        +"\n	#define DEF_DOMAIN_OFFSET_X "+ &format!("{:?}", (o_x as f32+((d_x>1) as i32 as f32) - 0.5*(d_x as f32 - 1.0) * (n_x - 2u32 * ((d_x>1u32) as i32 as u32)) as f32))+"f"
        +"\n	#define DEF_DOMAIN_OFFSET_Y "+ &format!("{:?}", (o_y as f32+((d_y>1) as i32 as f32) - 0.5*(d_y as f32 - 1.0) * (n_y - 2u32 * ((d_y>1u32) as i32 as u32)) as f32))+"f"
        +"\n	#define DEF_DOMAIN_OFFSET_Z "+ &format!("{:?}", (o_z as f32+((d_z>1) as i32 as f32) - 0.5*(d_z as f32 - 1.0) * (n_z - 2u32 * ((d_z>1u32) as i32 as u32)) as f32))+"f"

        +"\n	#define D"+ &dimensions.to_string()+"Q"+ &velocity_set.to_string()+"" // D2Q9/D3Q15/D3Q19/D3Q27
        +"\n	#define DEF_VELOCITY_SET "+ &velocity_set.to_string()+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
        +"\n	#define DEF_DIMENSIONS "  + &dimensions.to_string()+"u"   // number spatial dimensions (2D or 3D)
        +"\n	#define DEF_TRANSFERS "   + &transfers.to_string()+"u"    // number of DDFs that are transferred between multiple domains

        +"\n	#define DEF_C 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
        +"\n	#define DEF_W " + &format!("{:?}", 1.0f32/(3.0f32*nu+0.5f32))+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
        + match lbm_config.velocity_set {
            VelocitySet::D2Q9 => &d2q9,
            VelocitySet::D3Q15 => &d3q15,
            VelocitySet::D3Q19 => &d3q19,
            VelocitySet::D3Q27 => &d3q27
        }
        + match lbm_config.relaxation_time {
            RelaxationTime::Srt => srt,
            RelaxationTime::Trt => trt
        }
        +"\n	#define TYPE_S 0x01" // 0b00000001 // (stationary or moving) solid boundary
        +"\n	#define TYPE_E 0x02" // 0b00000010 // equilibrium boundary (inflow/outflow)
        +"\n	#define TYPE_T 0x04" // 0b00000100 // temperature boundary
        +"\n	#define TYPE_F 0x08" // 0b00001000 // fluid
        +"\n	#define TYPE_I 0x10" // 0b00010000 // interface
        +"\n	#define TYPE_G 0x20" // 0b00100000 // gas
        +"\n	#define TYPE_X 0x40" // 0b01000000 // reserved type X
        +"\n	#define TYPE_Y 0x80" // 0b10000000 // reserved type Y

        +"\n	#define TYPE_MS 0x03" // 0b00000011 // cell next to moving solid boundary
        +"\n	#define TYPE_BO 0x03" // 0b00000011 // any flag bit used for boundaries (temperature excluded)
        +"\n	#define TYPE_IF 0x18" // 0b00011000 // change from interface to fluid
        +"\n	#define TYPE_IG 0x30" // 0b00110000 // change from interface to gas
        +"\n	#define TYPE_GI 0x38" // 0b00111000 // change from gas to interface
        +"\n	#define TYPE_SU 0x38" // 0b00111000 // any flag bit used for SURFACE
        + match lbm_config.float_type { //Floatingpoint types
            FloatType::FP16S => &fp16s,
            FloatType::FP16C => &fp16c,
            FloatType::FP32 => &fp32
        }
        + if lbm_config.ext_equilibrium_boudaries {"\n	#define EQUILIBRIUM_BOUNDARIES"} else {""}
        + if lbm_config.ext_volume_force {         "\n	#define VOLUME_FORCE"} else {""}
        + &if lbm_config.ext_magneto_hydro {
         "\n	#define MAGNETO_HYDRO".to_owned()
        +"\n	#define DEF_KE "    + &format!("{:?}f", lbm_config.units.si_to_ke()) // coulomb constant scaled by distance per lattice cell
        +"\n	#define DEF_KMU "   + &format!("{:?}f", lbm_config.units.si_to_mu_0() / (4.0 * PI))
        +"\n	#define DEF_LOD_DEPTH " + &format!("{}u", lbm_config.mhd_lod_depth)
        +"\n    #define DEF_NUM_LOD " + &format!("{}u", n_lod)
        +"\n    #define DEF_NUM_LOD_OWN " + &format!("{}u", n_lod_own)
        +"\n	#define DEF_WQ  "  + &format!("{:?}f", 1.0/(2.0*lbm_config.units.si_to_k_charge_expansion()+0.5))
        } else {"".to_string()}
        + if lbm_config.ext_force_field {                "\n	#define FORCE_FIELD"} else {""}
        + if lbm_config.graphics_config.graphics_active {"\n	#define UPDATE_FIELDS"} else {""}
        //Extensions
    }

    /// Enqueues the `kernel_intialize` Kernel. Needs to be called after first setup to ready the simulation.
    pub fn enqueue_initialize(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, arguments are already set
        unsafe { self.kernel_initialize.cmd().enq()? }
        self.queue.finish()
    }

    /// Enqueues the `kernel_stream_collide` Kernel. This is the main simulation kernel and is executed every time step.
    pub fn enqueue_stream_collide(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, some arguments are already set
        unsafe {
            self.kernel_stream_collide.set_arg("t", self.t)?;
            self.kernel_stream_collide.set_arg("fx", self.fx)?;
            self.kernel_stream_collide.set_arg("fy", self.fy)?;
            self.kernel_stream_collide.set_arg("fz", self.fz)?;
            self.kernel_stream_collide.cmd().enq()
        }
    }

    /// Enqueues the `kernel_update_fields` Kernel. This functionality is automatically handled in the stream_collide Kernel if graphics are enabled.
    #[allow(unused)]
    pub fn enqueue_update_fields(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, arguments are already set
        unsafe {
            self.kernel_update_fields.set_arg("t", self.t)?;
            self.kernel_update_fields.set_arg("fx", self.fx)?;
            self.kernel_update_fields.set_arg("fy", self.fy)?;
            self.kernel_update_fields.set_arg("fz", self.fz)?;
            self.kernel_update_fields.cmd().enq()
        }
    }

    pub fn enqueue_update_e_b_dyn(&self) -> ocl::Result<()> {
        unsafe {
            self.kernel_update_e_b_dyn
                .as_ref()
                .expect("kernel should be initialized")
                .cmd()
                .enq()
        }
    }

    pub fn enqueue_clear_qu_lod(&self) -> ocl::Result<()> {
        unsafe {
            self.kernel_clear_qu_lod
                .as_ref()
                .expect("kernel should be initialized")
                .cmd()
                .enq()
        }
    }

    // Transfer fields
    pub fn get_area(&self, direction: u32) -> usize {
        let a = [
            self.n_y as usize * self.n_z as usize,
            self.n_x as usize * self.n_z as usize,
            self.n_x as usize * self.n_y as usize,
        ];
        a[direction as usize]
    }
    /// Extract field
    pub fn enqueue_transfer_extract_field(
        &mut self,
        field: TransferField,
        direction: u32,
        bytes_per_cell: usize,
    ) -> ocl::Result<()> {
        // Enqueue extraction kernel
        let kernel = self.transfer_kernels[field as usize][0]
            .as_ref()
            .expect("transfer kernel should be initialized"); // [0] is the extraction kernel
        kernel.set_arg("direction", direction)?;
        kernel.set_arg("time_step", self.t)?;
        unsafe {
            kernel
                .cmd()
                .global_work_size(self.get_area(direction))
                .enq()?;
        }
        // Read result into host transfer buffers
        let field_length = self.get_area(direction) * bytes_per_cell;
        self.transfer_p
            .read(&mut self.transfer_p_host)
            .len(field_length)
            .enq()
            .unwrap();
        self.transfer_m
            .read(&mut self.transfer_m_host)
            .len(field_length)
            .enq()
    }

    /// Insert field
    pub fn enqueue_transfer_insert_field(
        &mut self,
        field: TransferField,
        direction: u32,
        bytes_per_cell: usize,
    ) -> ocl::Result<()> {
        // Write transfer data to device buffers
        let field_length = self.get_area(direction) * bytes_per_cell;
        self.transfer_p
            .write(&self.transfer_p_host)
            .len(field_length)
            .enq()?;
        self.transfer_m
            .write(&self.transfer_m_host)
            .len(field_length)
            .enq()?;
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

    // Transfer LODs
    /// Read own LODs into host buffer
    pub fn read_lods(&mut self) {
        self.qu_lod.as_ref().expect("msg").read(self.transfer_lod_host.as_mut().expect("msg")).len(self.n_lod_own*4).enq().unwrap();
    }

    #[allow(unused)]
    /// Debug function. Prints out all cell information in Lattice and SI units.
    ///
    /// This is really slow.
    pub fn dump_cell(&self, c: usize, cfg: &LbmConfig) {
        let n = self.n_x * self.n_y * self.n_z;

        let c_x = c;
        let c_y = c + n as usize;
        let c_z = c + n as usize * 2;
        let (x, y, z) = self.get_coordinates(c as u64);

        let mut rho: Vec<f32> = vec![0.0; n as usize];
        let mut u: Vec<f32> = vec![0.0; n as usize * 3];
        let mut flags: Vec<u8> = vec![0; n as usize];
        let mut q: Vec<f32> = vec![0.0; n as usize];
        let mut e: Vec<f32> = vec![0.0; n as usize * 3];
        let mut b: Vec<f32> = vec![0.0; n as usize * 3];
        let mut e_dyn: Vec<f32> = vec![0.0; n as usize * 3];
        let mut b_dyn: Vec<f32> = vec![0.0; n as usize * 3];

        self.rho.read(&mut rho).enq().unwrap();
        self.u.read(&mut u).enq().unwrap();
        self.flags.read(&mut flags).enq().unwrap();
        match &self.q {
            Some(q_s) => {
                q_s.read(&mut q).enq().unwrap();
            }
            None => {}
        }
        match &self.e {
            Some(e_s) => {
                e_s.read(&mut e).enq().unwrap();
            }
            None => {}
        }
        match &self.b {
            Some(b_s) => {
                b_s.read(&mut b).enq().unwrap();
            }
            None => {}
        }
        match &self.e_dyn {
            Some(e_s) => {
                e_s.read(&mut e_dyn).enq().unwrap();
            }
            None => {}
        }
        match &self.b_dyn {
            Some(b_s) => {
                b_s.read(&mut b_dyn).enq().unwrap();
            }
            None => {}
        }

        let rho_si = cfg.units.dens_to_si(rho[c]);
        let u_x_si = cfg.units.speed_to_si(u[c_x]);
        let u_y_si = cfg.units.speed_to_si(u[c_y]);
        let u_z_si = cfg.units.speed_to_si(u[c_z]);
        let b_x_si = cfg.units.mag_flux_to_si(b[c_x]);
        let b_y_si = cfg.units.mag_flux_to_si(b[c_y]);
        let b_z_si = cfg.units.mag_flux_to_si(b[c_z]);
        let b_dyn_x_si = cfg.units.mag_flux_to_si(b_dyn[c_x]);
        let b_dyn_y_si = cfg.units.mag_flux_to_si(b_dyn[c_y]);
        let b_dyn_z_si = cfg.units.mag_flux_to_si(b_dyn[c_z]);
        let e_x_si = cfg.units.e_field_to_si(e[c_x]);
        let e_y_si = cfg.units.e_field_to_si(e[c_y]);
        let e_z_si = cfg.units.e_field_to_si(e[c_z]);
        let e_dyn_x_si = cfg.units.e_field_to_si(e_dyn[c_x]);
        let e_dyn_y_si = cfg.units.e_field_to_si(e_dyn[c_y]);
        let e_dyn_z_si = cfg.units.e_field_to_si(e_dyn[c_z]);
        let charge = cfg.units.charge_to_si(q[c]);

        println!(
            "Dumping cell {}:
    x: {}, y: {}, z: {}
    rho:    {} / {} kg/m³
    u:      {}, {}, {} / {}, {}, {} m/s
    flags:  {}
    e:      {}, {}, {} / {}, {}, {} V/m
    b:      {}, {}, {} / {}, {}, {} T
    e_dyn:  {}, {}, {} / {}, {}, {} V/m
    b_dyn:  {}, {}, {} / {}, {}, {} T
    charge: {} / {} As",
            c,
            x,
            y,
            z,
            rho[c],
            rho_si,
            u[c_x],
            u[c_y],
            u[c_z],
            u_x_si,
            u_y_si,
            u_z_si,
            flags[c],
            e[c_x],
            e[c_y],
            e[c_z],
            e_x_si,
            e_y_si,
            e_z_si,
            b[c_x],
            b[c_y],
            b[c_z],
            b_x_si,
            b_y_si,
            b_z_si,
            e_dyn[c_x],
            e_dyn[c_y],
            e_dyn[c_z],
            e_dyn_x_si,
            e_dyn_y_si,
            e_dyn_z_si,
            b_dyn[c_x],
            b_dyn[c_y],
            b_dyn[c_z],
            b_dyn_x_si,
            b_dyn_y_si,
            b_dyn_z_si,
            q[c],
            charge
        )
    }

    /// Get total simulation size.
    fn get_n(n_x: u32, n_y: u32, n_z: u32) -> u64 {
        n_x as u64 * n_y as u64 * n_z as u64
    }

    #[allow(dead_code)]
    pub fn is_halo(&mut self, x: u32, y: u32, z: u32) -> bool {
	    return ((self.cfg.d_x>1)&(x==0||x>=self.cfg.n_x-1))||((self.cfg.d_y>1)&(y==0||y>=self.cfg.n_y-1))||((self.cfg.d_z>1)&(z==0||z>=self.cfg.n_z-1));
    }

    /// Get `x, y, z` coordinates from 1D index `n`.
    pub fn get_coordinates(&self, n: u64) -> (u32, u32, u32) {
        get_coordinates_sl(n, self.n_x, self.n_y)
    }
}
