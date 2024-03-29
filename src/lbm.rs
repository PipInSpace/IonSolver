use std::cmp;

use crate::*;
use crate::{graphics::GraphicsConfig, units::Units};
use ocl::{Buffer, Context, Device, Kernel, Platform, Program, Queue};

/// Velocity discretizations in 2D and 3D.
/// ```
/// D2Q9:  2D
/// D3Q15: 3D low precision
/// D3Q19: 3D recommended
/// D3Q27: 3D highest precision
/// ```
#[allow(dead_code)]
#[derive(Clone, Copy, Default)]
pub enum VelocitySet {
    #[default]
    /// 2D
    D2Q9,
    /// 3D low precision
    D3Q15,
    /// 3D recommended
    D3Q19,
    /// 3D highest precision
    D3Q27,
}

impl VelocitySet {
    fn get_transfers(&self) -> usize {
        match self {
            VelocitySet::D2Q9 => 3_usize,
            VelocitySet::D3Q15 => 5_usize,
            VelocitySet::D3Q19 => 5_usize,
            VelocitySet::D3Q27 => 9_usize,
        }
    }
    fn get_set_values(&self) -> (u8, u8, u8) {
        match self {
            //Set dimensions/velocitys/transfers from Enum
            VelocitySet::D2Q9 => (2, 9, 3),
            VelocitySet::D3Q15 => (3, 15, 5),
            VelocitySet::D3Q19 => (3, 19, 5),
            VelocitySet::D3Q27 => (3, 27, 9),
        }
    }
}

/// LBM relaxation time type.
/// ```
/// Srt: Single relaxation time type, more efficient
/// Trt: Two-relaxation time type, more precise
/// ```
#[allow(dead_code)]
#[derive(Clone, Copy, Default)]
pub enum RelaxationTime {
    #[default]
    /// Single relaxation time, more efficient
    Srt,
    /// Two-relaxation time, more precise
    Trt,
}

/// Enum for different floating-point number types used in the simulation.
///
/// Types: `FP32`,`FP16S`,`FP16C`
///
/// `FP32` represents the normal floating-point number type `f32`. It takes the most memory.
///
/// `FP16S` and `FP16C` are custom floating-point number types, represented as `u16`. They take less memory.
/// `FP16S` is recommended for best precision/memory footprint.
///
/// [Learn more at this paper about custom float types.](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats)
#[allow(dead_code)]
#[derive(Clone, Copy, Default)]
pub enum FloatType {
    #[default]
    /// Custom float type represented as a u16, recommended
    FP16S,
    /// Custom float type represented as a u16
    FP16C,
    /// Default float type
    FP32,
}

impl FloatType {
    fn size_of(&self) -> usize {
        match self {
            FloatType::FP16S => 2_usize,
            FloatType::FP16C => 2_usize,
            FloatType::FP32 => 4_usize,
        }
    }
}

/// Enum representing a buffer of variable [`FloatType`]
///
/// [`FloatType`]: crate::lbm::FloatType
#[derive(Clone)]
pub enum VariableFloatBuffer {
    U16(Buffer<u16>), // Buffers for variable float types
    F32(Buffer<f32>),
}

/// Enum to identify a field transfered between domain boundaries
#[derive(Clone, Copy)]
enum TransferField {
    Fi,
    RhoUFlags,
}

/// Struct used to bundle arguments for LBM simulation setup.
/// ```
/// velocity_set: VelocitySet, // Velocity discretization mode
/// relaxation_time: RelaxationTime, // Sim relaxation time type
/// float_type: FloatType, // Sim float type for memory compression of ddf's
/// n_x: u32, // Sim size on each axis
/// n_y: u32,
/// n_z: u32,
///
/// d_x: u32, // Sim domains on each axis
/// d_y: u32, // Multi-domain support is planned, but only partially implemented currently
/// d_z: u32,
///
/// nu: f32, // Kinematic shear viscosity
/// fx: f32, // Volume force on each axis
/// fy: f32,
/// fz: f32,
///
/// ext_equilibrium_boudaries: bool, // Extensions, provide additional functionality
/// ext_volume_force: bool,
/// ext_force_field: bool,    // Needs volume_force to work
/// ext_electro_hydro: bool, // Needs force_field to work
///
/// graphics_config: GraphicsConfig, // Grapics config struct
/// ```
#[derive(Clone)]
pub struct LbmConfig {
    // Holds init information about the simulation
    pub velocity_set: VelocitySet,
    pub relaxation_time: RelaxationTime,
    pub float_type: FloatType,
    pub units: Units,
    pub n_x: u32, //Size
    pub n_y: u32,
    pub n_z: u32,

    pub d_x: u32, //Domains
    pub d_y: u32,
    pub d_z: u32,

    pub nu: f32,
    pub fx: f32,
    pub fy: f32,
    pub fz: f32,

    pub ext_equilibrium_boudaries: bool, //Extensions
    pub ext_volume_force: bool,
    pub ext_force_field: bool,   // Needs volume_force to work
    pub ext_magneto_hydro: bool, // Needs volume_force to work

    /// Cell range of each cells induction (Keep this small)
    pub induction_range: u8,
    /// Charge per density (unit: (A/s)/(kg/m^3))
    pub charge_per_dens: f32,

    pub graphics_config: GraphicsConfig,
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
            fx: 0.0f32,
            fy: 0.0f32,
            fz: 0.0f32,

            ext_equilibrium_boudaries: false,
            ext_volume_force: false,
            ext_magneto_hydro: false,
            ext_force_field: false,

            induction_range: 5, // Range of the cells induction (Keep this small)
            charge_per_dens: 0.0,

            graphics_config: graphics::GraphicsConfig::new(),
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
    pub domains: Vec<LbmDomain>,
    pub config: LbmConfig,
    initialized: bool,
}

impl Lbm {
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
            ))
        }
        println!("All domains initialized.\n");

        Lbm {
            domains: lbm_domains,
            config: lbm_config,
            initialized: false,
        }
    }

    /// Readies the LBM Simulation to be run.
    /// Executes `kernel_initialize` Kernels for every `LbmDomain` and fills domain transfer buffers.
    pub fn initialize(&mut self) {
        // the communicate calls at initialization need an odd time step
        self.increment_timestep(1);
        //communicate_rho_u_flags
        for d in 0..self.get_domain_numbers() {
            self.domains[d].enqueue_initialize().unwrap();
        }
        //communicate_rho_u_flags
        //communicate_fi
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
        // call kernel stream_collide to perform one LBM time step
        self.stream_collide();
        if self.config.graphics_config.graphics_active {
            self.communicate_rho_u_flags();
        }
        self.communicate_fi();
        
        if self.config.ext_magneto_hydro && self.config.induction_range != 0 {
            self.update_e_b_dynamic();
        }

        self.finish_queues();
        self.increment_timestep(1);
    }

    /// Blocks execution until all `LbmDomain` OpenCL queues have finished.
    pub fn finish_queues(&self) {
        for d in 0..self.domains.len() {
            self.domains[d].queue.finish().unwrap();
        }
    }

    /// Execute stream_collide kernel on all domains
    fn stream_collide(&self) {
        for d in 0..self.get_domain_numbers() {
            self.domains[d].enqueue_stream_collide().unwrap();
        }
    }

    /// Execute E and B dynamic update kernel on all domains
    fn update_e_b_dynamic(&self) {
        for d in 0..self.get_domain_numbers() {
            self.domains[d].enqueue_update_e_b_dyn().unwrap();
        }
    }

    // Domain communication
    /// Communicate a field across domain barriers
    #[rustfmt::skip]
    fn communicate_field(&mut self, field: TransferField, bytes_per_cell: usize) {
        let d_x = self.config.d_x as usize;
        let d_y = self.config.d_y as usize;
        let d_z = self.config.d_z as usize;
        let d_n = self.get_domain_numbers();

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
            self.config.float_type.size_of() * self.config.velocity_set.get_transfers();
        self.communicate_field(TransferField::Fi, bytes_per_cell);
    }
    /// Communicate rho, u and flags across domain boundaries (needed for graphics)
    fn communicate_rho_u_flags(&mut self) {
        self.communicate_field(TransferField::RhoUFlags, 17);
    }

    // Helper functions

    /// Increments time steps variable for every `LbmDomain`
    fn increment_timestep(&mut self, steps: u32) {
        for d in 0..self.domains.len() {
            self.domains[d].t += steps as u64;
        }
    }

    /// Resets tims steps variable for every `LbmDomain`
    fn reset_timestep(&mut self) {
        for d in 0..self.domains.len() {
            self.domains[d].t = 0;
        }
    }

    pub fn get_domain_numbers(&self) -> usize {
        self.domains.len()
    }

    /// Get `x, y, z` coordinates from 1D index `n`.
    #[allow(unused)]
    pub fn get_coordinates(&self, n: u64) -> (u32, u32, u32) {
        // disassembles 1D index into 3D index
        let t: u64 = n % (self.config.n_x as u64 * self.config.n_y as u64);
        //x, y, z
        (
            (t % self.config.n_x as u64) as u32,
            (t / self.config.n_x as u64) as u32,
            (n / (self.config.n_x as u64 * self.config.n_y as u64)) as u32,
        )
    }
}

/// The `LbmDomain` struct holds all information to run a LBM-Simulation on one OpenCL Device.
/// It is initialized with:
/// ```
/// LbmDomain::new(lbm_config: LbmConfig, device: Device, x: u32, y: u32, z: u32)
/// ```
/// `LbmDomain` should not be initialized on it's own, but automatically through the `Lbm::new()` function.
/// This ensures all arguments are correctly set.
pub struct LbmDomain {
    // FluidX3D creates contexts/queues/programs for each device automatically through another class
    pub queue: Queue,

    kernel_initialize: Kernel, // Basic Kernels
    kernel_stream_collide: Kernel,
    kernel_update_fields: Kernel,
    transfer_kernels: [[Kernel; 2]; 2],

    kernel_update_e_b_dyn: Option<Kernel>, // Optional Kernels

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

    pub e: Option<Buffer<f32>>, // Optional Buffers
    pub b: Option<Buffer<f32>>,
    pub e_dyn: Option<Buffer<f32>>,
    pub b_dyn: Option<Buffer<f32>>,
    pub f: Option<Buffer<f32>>,
    pub t: u64, // Timestep

    pub graphics: Option<graphics::Graphics>, // Graphics struct, handles rendering
}

impl LbmDomain {
    /// Returns new `LbmDomain` from provided arguments.
    ///
    /// `LbmDomain` should not be initialized on it's own, but automatically through the `Lbm::new()` function.
    /// This ensures all arguments are correctly set.
    pub fn new(lbm_config: &LbmConfig, device: Device, x: u32, y: u32, z: u32) -> LbmDomain {
        let n_x = lbm_config.n_x / lbm_config.d_x + 2u32 * (lbm_config.d_x > 1u32) as u32; // Size + Halo offsets
        let n_y = lbm_config.n_y / lbm_config.d_y + 2u32 * (lbm_config.d_y > 1u32) as u32; // When multiple domains on axis -> add 2 cells of padding
        let n_z = lbm_config.n_z / lbm_config.d_z + 2u32 * (lbm_config.d_z > 1u32) as u32;
        let n = Self::get_n(n_x, n_y, n_z);

        let d_x = lbm_config.d_x;
        let d_y = lbm_config.d_y;
        let d_z = lbm_config.d_z;

        let o_x = (x * n_x / d_x) as i32 - (lbm_config.d_x > 1u32) as i32;
        let o_y = (y * n_y / d_y) as i32 - (lbm_config.d_y > 1u32) as i32;
        let o_z = (z * n_z / d_z) as i32 - (lbm_config.d_z > 1u32) as i32;

        let (dimensions, velocity_set, transfers) = lbm_config.velocity_set.get_set_values();

        let t = 0;

        let ocl_code = Self::get_device_defines(
            n_x,
            n_y,
            n_z,
            d_x,
            d_y,
            d_z,
            o_x,
            o_y,
            o_z,
            dimensions,
            velocity_set,
            transfers,
            lbm_config.nu,
            lbm_config,
        ) + &if lbm_config.graphics_config.graphics_active {
            graphics::get_graphics_defines(lbm_config.graphics_config)
        } else {
            "".to_string()
        } + &opencl::get_opencl_code(); // Only appends graphics defines if needed

        // OCL variables are directly exposed, due to no other device struct.
        let platform = Platform::default();
        let context = Context::builder()
            .platform(platform)
            .devices(device)
            .build()
            .unwrap();
        let queue = Queue::new(&context, device, None).unwrap();
        println!("    Compiling Program...");
        let program = Program::builder()
            .devices(device)
            .src(&ocl_code)
            .build(&context)
            .unwrap();

        // Initialize Buffers
        let fi: VariableFloatBuffer = match lbm_config.float_type {
            FloatType::FP32 => {
                // Float Type F32
                VariableFloatBuffer::F32(opencl::create_buffer(
                    &queue,
                    [n * velocity_set as u64],
                    0.0f32,
                ))
            }
            _ => {
                // Float Type F16S/F16C
                VariableFloatBuffer::U16(opencl::create_buffer(
                    &queue,
                    [n * velocity_set as u64],
                    0u16,
                ))
            }
        };
        let rho = opencl::create_buffer(&queue, [n], 1.0f32);
        let u = opencl::create_buffer(&queue, [n * 3], 0f32);
        let flags = opencl::create_buffer(&queue, [n], 0u8);

        // Force field buffer as 3D Vectors
        let f: Option<Buffer<f32>> = if lbm_config.ext_force_field {
            Some(opencl::create_buffer(&queue, [n * 3], 0f32))
        } else {
            None
        };
        // Electric field buffer as 3D Vectors
        let e: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro {
            Some(opencl::create_buffer(&queue, [n * 3], 0f32))
        } else {
            None
        };
        let e_dyn: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro {
            Some(opencl::create_buffer(&queue, [n * 3], 0f32))
        } else {
            None
        };
        // Magnetic field buffers as 3D Vectors
        let b: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro {
            Some(opencl::create_buffer(&queue, [n * 3], 0f32))
        } else {
            None
        };
        let b_dyn: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro {
            Some(opencl::create_buffer(&queue, [n * 3], 0f32))
        } else {
            None
        };

        // Initialize Kernels
        let mut kernel_initialize_builder = Kernel::builder();
        let mut kernel_stream_collide_builder = Kernel::builder();
        let mut kernel_update_fields_builder = Kernel::builder();
        kernel_initialize_builder
            .program(&program)
            .name("initialize")
            .queue(queue.clone())
            .global_work_size([n]);
        kernel_stream_collide_builder
            .program(&program)
            .name("stream_collide")
            .queue(queue.clone())
            .global_work_size([n]);
        kernel_update_fields_builder
            .program(&program)
            .name("update_fields")
            .queue(queue.clone())
            .global_work_size([n]);
        match &fi {
            //Initialize kernels. Different Float types need different arguments (Fi-Buffer specifically)
            VariableFloatBuffer::F32(fif32) => {
                // Float Type F32
                kernel_initialize_builder.arg_named("fi", fif32);
                kernel_stream_collide_builder.arg_named("fi", fif32);
                kernel_update_fields_builder.arg_named("fi", fif32);
            }
            VariableFloatBuffer::U16(fiu16) => {
                // Float Type F16S/F16C
                kernel_initialize_builder.arg_named("fi", fiu16);
                kernel_stream_collide_builder.arg_named("fi", fiu16);
                kernel_update_fields_builder.arg_named("fi", fiu16);
            }
        }
        kernel_initialize_builder
            .arg_named("rho", &rho)
            .arg_named("u", &u)
            .arg_named("flags", &flags);
        kernel_stream_collide_builder
            .arg_named("rho", &rho)
            .arg_named("u", &u)
            .arg_named("flags", &flags)
            .arg_named("t", t)
            .arg_named("fx", lbm_config.fx)
            .arg_named("fy", lbm_config.fy)
            .arg_named("fz", lbm_config.fz);
        kernel_update_fields_builder
            .arg_named("rho", &rho)
            .arg_named("u", &u)
            .arg_named("flags", &flags)
            .arg_named("t", t)
            .arg_named("fx", lbm_config.fx)
            .arg_named("fy", lbm_config.fy)
            .arg_named("fz", lbm_config.fz);
        // Conditional arguments. Place at end of kernel functions
        if lbm_config.ext_force_field {
            kernel_stream_collide_builder
                .arg_named("F", f.as_ref().expect("f buffer used but not initialized"));
        }

        let mut kernel_update_e_b_dyn: Option<Kernel> = None;

        if lbm_config.ext_magneto_hydro {
            kernel_stream_collide_builder
                .arg_named("E", e.as_ref().expect("e buffer used but not initialized"))
                .arg_named("B", b.as_ref().expect("b buffer used but not initialized"))
                .arg_named(
                    "E_dyn",
                    e_dyn
                        .as_ref()
                        .expect("e_dyn buffer used but not initialized"),
                )
                .arg_named(
                    "B_dyn",
                    b_dyn
                        .as_ref()
                        .expect("b_dyn buffer used but not initialized"),
                );
            kernel_initialize_builder
                .arg_named("E", e.as_ref().expect("e buffer used but not initialized"))
                .arg_named("B", b.as_ref().expect("b buffer used but not initialized"))
                .arg_named(
                    "E_dyn",
                    e_dyn
                        .as_ref()
                        .expect("e_dyn buffer used but not initialized"),
                )
                .arg_named(
                    "B_dyn",
                    b_dyn
                        .as_ref()
                        .expect("b_dyn buffer used but not initialized"),
                );

            // Dynamic B kernel
            kernel_update_e_b_dyn = Some(
                Kernel::builder()
                    .program(&program)
                    .name("update_e_b_dynamic")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named(
                        "E_dyn",
                        e_dyn
                            .as_ref()
                            .expect("e_dyn buffer used but not initialized"),
                    )
                    .arg_named(
                        "B_dyn",
                        b_dyn
                            .as_ref()
                            .expect("b_dyn buffer used but not initialized"),
                    )
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("flags", &flags)
                    .build()
                    .unwrap(),
            )
        }

        let kernel_initialize: Kernel = kernel_initialize_builder.build().unwrap();
        let kernel_stream_collide: Kernel = kernel_stream_collide_builder.build().unwrap();
        let kernel_update_fields: Kernel = kernel_update_fields_builder.build().unwrap();

        // Multi-Domain-Transfers:
        // Transfer buffer initializaton:
        let mut a_max: usize = 0;
        if lbm_config.d_x > 1 {
            a_max = cmp::max(a_max, n_y as usize * n_z as usize);
        } // Ax
        if lbm_config.d_y > 1 {
            a_max = cmp::max(a_max, n_x as usize * n_z as usize);
        } // Ay
        if lbm_config.d_z > 1 {
            a_max = cmp::max(a_max, n_x as usize * n_y as usize);
        } // Az
        let transfer_buffer_size =
            a_max * cmp::max(17, transfers as usize * lbm_config.float_type.size_of());

        let transfer_p = opencl::create_buffer(&queue, transfer_buffer_size, 0_u8); // Size is maxiumum 17 bytes or bytes needed for fi
        let transfer_m = opencl::create_buffer(&queue, transfer_buffer_size, 0_u8); // Size is maxiumum 17 bytes or bytes needed for fi
        let transfer_p_host: Vec<u8> = vec![0_u8; transfer_buffer_size];
        let transfer_m_host: Vec<u8> = vec![0_u8; transfer_buffer_size];

        // Transfer kernels
        let transfer_kernels = [
            [
                match &fi {
                    VariableFloatBuffer::U16(fi_u16) => Kernel::builder()
                        .program(&program)
                        .name("transfer_extract_fi")
                        .queue(queue.clone())
                        .global_work_size(1)
                        .arg_named("direction", 0_u32)
                        .arg_named("time_step", 0_u64)
                        .arg_named("transfer_buffer_p", &transfer_p)
                        .arg_named("transfer_buffer_m", &transfer_m)
                        .arg_named("fi", fi_u16)
                        .build()
                        .unwrap(),
                    VariableFloatBuffer::F32(fi_f32) => Kernel::builder()
                        .program(&program)
                        .name("transfer_extract_fi")
                        .queue(queue.clone())
                        .global_work_size(1)
                        .arg_named("direction", 0_u32)
                        .arg_named("time_step", 0_u64)
                        .arg_named("transfer_buffer_p", &transfer_p)
                        .arg_named("transfer_buffer_m", &transfer_m)
                        .arg_named("fi", fi_f32)
                        .build()
                        .unwrap(),
                }, // Extract Fi
                match &fi {
                    VariableFloatBuffer::U16(fi_u16) => Kernel::builder()
                        .program(&program)
                        .name("transfer__insert_fi")
                        .queue(queue.clone())
                        .global_work_size(1)
                        .arg_named("direction", 0_u32)
                        .arg_named("time_step", 0_u64)
                        .arg_named("transfer_buffer_p", &transfer_p)
                        .arg_named("transfer_buffer_m", &transfer_m)
                        .arg_named("fi", fi_u16)
                        .build()
                        .unwrap(),
                    VariableFloatBuffer::F32(fi_f32) => Kernel::builder()
                        .program(&program)
                        .name("transfer__insert_fi")
                        .queue(queue.clone())
                        .global_work_size(1)
                        .arg_named("direction", 0_u32)
                        .arg_named("time_step", 0_u64)
                        .arg_named("transfer_buffer_p", &transfer_p)
                        .arg_named("transfer_buffer_m", &transfer_m)
                        .arg_named("fi", fi_f32)
                        .build()
                        .unwrap(),
                }, // Insert Fi
            ], // Fi
            [
                Kernel::builder()
                    .program(&program)
                    .name("transfer_extract_rho_u_flags")
                    .queue(queue.clone())
                    .global_work_size(1)
                    .arg_named("direction", 0_u32)
                    .arg_named("time_step", 0_u64)
                    .arg_named("transfer_buffer_p", &transfer_p)
                    .arg_named("transfer_buffer_m", &transfer_m)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("flags", &flags)
                    .build()
                    .unwrap(), // Extract rho, u, flags
                Kernel::builder()
                    .program(&program)
                    .name("transfer__insert_rho_u_flags")
                    .queue(queue.clone())
                    .global_work_size(1)
                    .arg_named("direction", 0_u32)
                    .arg_named("time_step", 0_u64)
                    .arg_named("transfer_buffer_p", &transfer_p)
                    .arg_named("transfer_buffer_m", &transfer_m)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("flags", &flags)
                    .build()
                    .unwrap(), // Extract rho, u, flags
            ], // Rho, u and flags (needed for graphics)
        ];
        println!("    Kernels for domain compiled.");
        //TODO: allocate transfer buffers

        let graphics: Option<graphics::Graphics> = if lbm_config.graphics_config.graphics_active {
            Some(graphics::Graphics::new(
                lbm_config, &program, &queue, &flags, &u,
            ))
        } else {
            None
        };

        LbmDomain {
            queue,

            kernel_initialize,
            kernel_stream_collide,
            kernel_update_fields,
            transfer_kernels,

            kernel_update_e_b_dyn,

            n_x,
            n_y,
            n_z,

            fx: lbm_config.fx,
            fy: lbm_config.fy,
            fz: lbm_config.fz,

            fi,
            rho,
            u,
            flags,

            transfer_p,
            transfer_m,
            transfer_p_host,
            transfer_m_host,

            e,
            b,
            e_dyn,
            b_dyn,
            f,
            t,

            graphics,
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
        o_x: i32,
        o_y: i32,
        o_z: i32,
        dimensions: u8,
        velocity_set: u8,
        transfers: u8,
        nu: f32,
        lbm_config: &LbmConfig,
    ) -> String {
        //Conditional Defines:
        //Velocity set types
        let d2q9 = "\n	#define def_w0 (1.0f/2.25f)".to_owned() // center (0)
        +"\n	#define def_ws (1.0f/9.0f)" // straight (1-4)
        +"\n	#define def_we (1.0f/36.0f)";
        let d3q15 = "\n	#define def_w0 (1.0f/4.5f)".to_owned() // center (0)
        +"\n	#define def_ws (1.0f/9.0f)" // straight (1-6)
        +"\n	#define def_wc (1.0f/72.0f)";
        let d3q19 = "\n	#define def_w0 (1.0f/3.0f)".to_owned() // center (0)
        +"\n	#define def_ws (1.0f/18.0f)" // straight (1-6)
        +"\n	#define def_we (1.0f/36.0f)";
        let d3q27 = "\n	#define def_w0 (1.0f/3.0f)".to_owned() // center (0)
        +"\n	#define def_ws (1.0f/18.0f)" // straight (1-6)
        +"\n	#define def_we (1.0f/36.0f)";
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
        "\n    #define def_Nx ".to_owned() + &n_x.to_string()+"u"
        +"\n	#define def_Ny "+ &n_y.to_string()+"u"
        +"\n	#define def_Nz "+ &n_z.to_string()+"u"
        +"\n	#define def_N  "+ &Self::get_n(n_x, n_y, n_z).to_string()+"ul"

        +"\n	#define def_Dx "+ &d_x.to_string()+"u"
        +"\n	#define def_Dy "+ &d_y.to_string()+"u"
        +"\n	#define def_Dz "+ &d_z.to_string()+"u"

        +"\n	#define def_Ox "+ &o_x.to_string()+"" // offsets are signed integer!
        +"\n	#define def_Oy "+ &o_y.to_string()+""
        +"\n	#define def_Oz "+ &o_z.to_string()+""

        +"\n	#define def_Ax "+ &(n_y * n_z).to_string()+"u"
        +"\n	#define def_Ay "+ &(n_z * n_x).to_string()+"u"
        +"\n	#define def_Az "+ &(n_x * n_y).to_string()+"u"

        +"\n	#define def_domain_offset_x "+ &format!("{:?}", (o_x as f32+(d_x>1) as i32 as f32 - 0.5*(d_x as f32 - 1.0) * (n_x - 2u32 * (d_x>1u32) as i32 as u32) as f32))+"f"
        +"\n	#define def_domain_offset_y "+ &format!("{:?}", (o_y as f32+(d_y>1) as i32 as f32 - 0.5*(d_y as f32 - 1.0) * (n_y - 2u32 * (d_y>1u32) as i32 as u32) as f32))+"f"
        +"\n	#define def_domain_offset_z "+ &format!("{:?}", (o_z as f32+(d_z>1) as i32 as f32 - 0.5*(d_z as f32 - 1.0) * (n_z - 2u32 * (d_z>1u32) as i32 as u32) as f32))+"f"

        +"\n	#define D"+ &dimensions.to_string()+"Q"+ &velocity_set.to_string()+"" // D2Q9/D3Q15/D3Q19/D3Q27
        +"\n	#define def_velocity_set "+ &velocity_set.to_string()+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
        +"\n	#define def_dimensions "+ &dimensions.to_string()+"u" // number spatial dimensions (2D or 3D)
        +"\n	#define def_transfers "+ &transfers.to_string()+"u" // number of DDFs that are transferred between multiple domains

        +"\n	#define def_c 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
        +"\n	#define def_w " + &format!("{:?}", 1.0f32/(3.0f32*nu+0.5f32))+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
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
        + if lbm_config.ext_volume_force {"\n	        #define VOLUME_FORCE"} else {""}
        + &if lbm_config.ext_magneto_hydro {"\n	        #define MAGNETO_HYDRO".to_owned()
        +"\n	#define def_ke "+ &format!("{:?}f", lbm_config.units.si_to_ke()) // coulomb constant scaled by distance per lattice cell
        +"\n	#define def_kmu "+ &format!("{:?}f", lbm_config.units.si_to_mu_0() / (4.0 * PI))
        +"\n	#define def_charge "+ &format!("{:?}f", lbm_config.units.si_to_charge_per_dens(lbm_config.charge_per_dens)) // charge held per density unit
        +"\n	#define def_ind_r "+ &lbm_config.induction_range.to_string()
        } else {"".to_string()}
        + if lbm_config.ext_force_field {"\n	        #define FORCE_FIELD"} else {""}
        + if lbm_config.graphics_config.graphics_active {"\n	#define UPDATE_FIELDS"} else {""}
        //Extensions
    }

    /// Enqueues the `kernel_intialize` Kernel. Needs to be called after first setup to ready the simulation.
    fn enqueue_initialize(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, arguments are already set
        unsafe { self.kernel_initialize.cmd().enq()? }
        self.queue.finish()
    }

    /// Enqueues the `kernel_stream_collide` Kernel. This is the main simulation kernel and is executed every time step.
    fn enqueue_stream_collide(&self) -> ocl::Result<()> {
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
    fn enqueue_update_fields(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, arguments are already set
        unsafe {
            self.kernel_update_fields.set_arg("t", self.t)?;
            self.kernel_update_fields.set_arg("fx", self.fx)?;
            self.kernel_update_fields.set_arg("fy", self.fy)?;
            self.kernel_update_fields.set_arg("fz", self.fz)?;
            self.kernel_update_fields.cmd().enq()
        }
    }

    fn enqueue_update_e_b_dyn(&self) -> ocl::Result<()> {
        unsafe {
            self.kernel_update_e_b_dyn
                .as_ref()
                .expect("kernel should be initialized")
                .cmd()
                .enq()
        }
    }

    // Transfer fields
    fn get_area(&self, direction: u32) -> usize {
        let a = [
            self.n_y as usize * self.n_z as usize,
            self.n_x as usize * self.n_z as usize,
            self.n_x as usize * self.n_y as usize,
        ];
        a[direction as usize]
    }
    /// Extract field
    fn enqueue_transfer_extract_field(
        &mut self,
        field: TransferField,
        direction: u32,
        bytes_per_cell: usize,
    ) -> ocl::Result<()> {
        // Enqueue extraction kernel
        let kernel = &self.transfer_kernels[field as usize][0]; // [0] is the extraction kernel
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
        println!("Field_length: {}, Len: {}, Dir: {}", field_length, self.transfer_p_host.len(), direction);
        self.transfer_p
            //.create_sub_buffer(None, [0], field_length).unwrap()
            .read(&mut self.transfer_p_host)
            .enq().unwrap();
        self.transfer_m
            //.create_sub_buffer(None, [0], field_length)?
            .read(&mut self.transfer_m_host)
            .enq()
    }

    /// Insert field
    fn enqueue_transfer_insert_field(
        &mut self,
        field: TransferField,
        direction: u32,
        bytes_per_cell: usize,
    ) -> ocl::Result<()> {
        let field_length = self.get_area(direction) * bytes_per_cell;
        self.transfer_p
            //.create_sub_buffer(None, [0], field_length)?
            .write(&self.transfer_p_host)
            .enq()?;
        self.transfer_m
            //.create_sub_buffer(None, [0], field_length)?
            .write(&self.transfer_m_host)
            .enq()?;
        let kernel = &self.transfer_kernels[field as usize][1]; // [1] is the insertion kernel
        kernel.set_arg("direction", direction)?;
        kernel.set_arg("time_step", self.t)?;
        unsafe {
            kernel
                .cmd()
                .global_work_size(self.get_area(direction))
                .enq()
        }
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
        let mut e: Vec<f32> = vec![0.0; n as usize * 3];
        let mut b: Vec<f32> = vec![0.0; n as usize * 3];
        let mut e_dyn: Vec<f32> = vec![0.0; n as usize * 3];
        let mut b_dyn: Vec<f32> = vec![0.0; n as usize * 3];

        self.rho.read(&mut rho).enq().unwrap();
        self.u.read(&mut u).enq().unwrap();
        self.flags.read(&mut flags).enq().unwrap();
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

        let rho_si = cfg.units.mass_to_si(rho[c]);
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
        let charge = cfg.units.charge_to_si(rho[c] * 0.0000000005);

        println!(
            "Dumping cell {}:
    x: {}, y: {}, z: {}
    rho:    {} / {} kg
    u:      {}x, {}y {}z / {}x, {}y {}z m/s
    flags:  {}
    e:      {}x, {}y {}z / {}x, {}y {}z V/m
    b:      {}x, {}y {}z / {}x, {}y {}z T
    e_dyn:  {}x, {}y {}z / {}x, {}y {}z V/m
    b_dyn:  {}x, {}y {}z / {}x, {}y {}z T
    charge: {} / {} A/s",
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
            rho[c] * 0.0000000005,
            charge
        )
    }

    /// Get total simulation size.
    fn get_n(n_x: u32, n_y: u32, n_z: u32) -> u64 {
        n_x as u64 * n_y as u64 * n_z as u64
    }

    /// Get `x, y, z` coordinates from 1D index `n`.
    #[allow(unused)]
    pub fn get_coordinates(&self, n: u64) -> (u32, u32, u32) {
        let t: u64 = n % (self.n_x as u64 * self.n_y as u64);
        //x, y, z
        (
            (t % self.n_x as u64) as u32,
            (t / self.n_x as u64) as u32,
            (n / (self.n_x as u64 * self.n_y as u64)) as u32,
        )
    }
}
