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

use std::cmp;

mod domain;
pub mod graphics;
#[cfg(feature = "multi-node")]
pub mod multi_node;
pub mod precompute;
mod units;
mod types;

use crate::opencl;
use {domain::LbmDomain, graphics::GraphicsConfig};
pub use types::*;
pub use units::Units;

// Helper Functions
/// Get `x, y, z` coordinates from 1D index `n` and side lengths `n_x` and `n_y`.
fn get_coordinates_sl(n: u64, n_x: u32, n_y: u32) -> (u32, u32, u32) {
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
            self.communicate_fqi();
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
            self.communicate_fqi();
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
    fn communicate_fqi(&mut self) {
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
    +"\n	#define DEF_N  "+ &(n_x as u64 * n_y as u64 * n_z as u64).to_string()+"ul"
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
    +"\n	#define DEF_KMU "   + &format!("{:?}f", lbm_config.units.si_to_mu_0() / (4.0 * std::f32::consts::PI))
    +"\n	#define DEF_LOD_DEPTH " + &format!("{}u", lbm_config.mhd_lod_depth)
    +"\n    #define DEF_NUM_LOD " + &format!("{}u", n_lod)
    +"\n    #define DEF_NUM_LOD_OWN " + &format!("{}u", n_lod_own)
    +"\n	#define DEF_WQ  "  + &format!("{:?}f", 1.0/(2.0*lbm_config.units.si_to_k_charge_expansion()+0.5))
    } else {"".to_string()}
    + if lbm_config.ext_force_field {                "\n	#define FORCE_FIELD"} else {""}
    + if lbm_config.graphics_config.graphics_active {"\n	#define UPDATE_FIELDS"} else {""}
    //Extensions
}