use crate::{graphics::Graphics, graphics::GraphicsConfig, *};
use ocl::{flags, Buffer, Context, Device, Kernel, Platform, Program, Queue};

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum VelocitySet {
    D2Q9,  // 2D
    D3Q15, // 3D low precision
    D3Q19, // 3D recommended
    D3Q27, // 3D highest precision
}

impl Default for VelocitySet {
    fn default() -> Self {
        VelocitySet::D2Q9
    }
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum RelaxationTime {
    SRT, // Single relaxation time, more efficient
    TRT, // Two-relaxation time, more precise
}

impl Default for RelaxationTime {
    fn default() -> Self {
        RelaxationTime::SRT
    }
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum FloatType {
    FP16S, // Custom float type represented as a u16, recommended
    FP16C, // Custom float type represented as a u16
    FP32,  // Default float type
}

impl Default for FloatType {
    fn default() -> Self {
        FloatType::FP16S
    }
}

//#[allow(dead_code)]
#[derive(Clone)]
pub enum VariableFloatBuffer {
    U16(Buffer<u16>), // Buffers for variable float types
    F32(Buffer<f32>),
}

// ############################################# Code Structure #############################################
// To start a LBM simulation, one needs to initialise an Lbm struct using the Lbm::init() function.
// The init() function takes in another struct, the LbmConfig, which contains all necessary arguments.
// LbmConfig needs to be configured beforehand, then the init() function takes care of Lbm initialisation.
// The Lbm struct contains one or multiple LbmDomain structs.
// The Simulation is actually run on these LbmDomain structs, each Domain corresponding to an OpenCL device.
// The simulation can be split up into multiple domains. This enables multi-device parallelization.
// LbmDomain initialisation is handled automatically in Lbm::init() using LbmDomain::init()

#[derive(Clone, Copy, Default)]
pub struct LbmConfig {
    // Holds init information about the simulation
    pub velocity_set: VelocitySet,
    pub relaxation_time: RelaxationTime,
    pub float_type: FloatType,
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
    pub sigma: f32,
    pub alpha: f32,
    pub beta: f32,

    pub particles_n: u32,
    pub particles_rho: f32,

    pub ext_equilibrium_boudaries: bool, //Extensions
    pub ext_volume_force: bool,
    pub ext_force_field: bool, // Needs volume_force to work
    pub ext_electric_force: bool, // Needs force_field to work

    pub graphics_config: GraphicsConfig,
}

impl LbmConfig {
    pub fn new() -> LbmConfig {
        // provides a LbmConfig struct with default values
        LbmConfig {
            velocity_set: VelocitySet::D2Q9,
            relaxation_time: RelaxationTime::SRT,
            float_type: FloatType::FP16S,
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
            sigma: 0.0f32,
            alpha: 0.0f32,
            beta: 0.0f32,

            particles_n: 0,
            particles_rho: 0.0f32,

            ext_equilibrium_boudaries: false,
            ext_volume_force: false,
            ext_electric_force: false,
            ext_force_field: false,

            graphics_config: GraphicsConfig::new(),
        }
    }
}

pub struct Lbm {
    // main simulation struct. Holds one or multiple subdomains
    pub domains: Vec<LbmDomain>,
    pub config: LbmConfig,
    initialized: bool,
}

impl Lbm {
    pub fn init(mut lbm_config: LbmConfig) -> Lbm {
        //Returns new Lbm from config struct. Domain setup handled automatically.
        //Configures Domains
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

        let h_x = (lbm_config.d_x > 1u32) as u32; // Halo offsets
        let h_y = (lbm_config.d_y > 1u32) as u32;
        let h_z = (lbm_config.d_z > 1u32) as u32;

        let device_infos = opencl::device_selection(domain_numbers);
        //TODO: sanity check

        let mut lbm_domains: Vec<LbmDomain> = Vec::new();
        for d in 0..domain_numbers {
            println!("Initializing domain {}/{}", d + 1, domain_numbers);
            let x = (d % (lbm_config.d_x * lbm_config.d_y)) % lbm_config.d_x;
            let y = (d % (lbm_config.d_x * lbm_config.d_y)) / lbm_config.d_x;
            let z = d / (lbm_config.d_x * lbm_config.d_y);
            println!(
                "Using device {} for domain {}",
                device_infos[d as usize].name().unwrap(),
                d + 1
            );
            lbm_domains.push(LbmDomain::init(
                lbm_config,
                device_infos[d as usize],
                h_x,
                h_y,
                h_z,
                x,
                y,
                z,
            ))
        }
        println!("All domains initialized.");

        let lbm = Lbm {
            domains: lbm_domains,
            config: lbm_config,
            initialized: false,
        };
        lbm
    }

    pub fn initialize(&mut self) {
        for d in 0..self.get_domain_numbers() {
            // the communicate calls at initialization need an odd time step
            self.domains[d].increment_timestep(1);
        }
        //communicate_rho_u_flags
        for d in 0..self.get_domain_numbers() {
            self.domains[d].enqueue_initialize().unwrap();
        }
        //communicate_rho_u_flags
        //communicate_fi
        for d in 0..self.get_domain_numbers() {
            self.domains[d].queue.finish().unwrap();
        }
        for d in 0..self.get_domain_numbers() {
            self.domains[d].reset_timestep();
        }
        self.initialized = true;
    }

    #[allow(unused)]
    pub fn run(&mut self, steps: u64) {
        //Initialize, then run simulation for steps
        //TODO: Display info in command line
        if !self.initialized {
            //Run initialization Kernel
            self.initialize();
            //TODO: Display initialization info in comman line
        }
        for i in 0..steps {
            println!("Step {}", i);
            self.do_time_step();
        }
    }

    pub fn do_time_step(&mut self) {
        // call kernel stream_collide to perform one LBM time step
        for d in 0..self.get_domain_numbers() {
            self.domains[d].enqueue_stream_collide().unwrap();
        }
        //communicate_fi
        for d in 0..self.get_domain_numbers() {
            self.domains[d].queue.finish().unwrap();
        }
        //for d in 0..self.get_domain_numbers() {
        //    self.domains[d].enqueue_update_fields().unwrap();
        //}
        //for d in 0..self.get_domain_numbers() {
        //    self.domains[d].queue.finish().unwrap();
        //}
        for d in 0..self.get_domain_numbers() {
            self.domains[d].increment_timestep(1);
        }
    }

    pub fn get_domain_numbers(&self) -> usize {
        self.domains.len()
    }

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

pub struct LbmDomain {
    device: Device, // FluidX3D creates contexts/queues/programs for each device automatically through another class
    context: Context, // IonSolver creates and saves this information directly.
    program: Program,
    pub queue: Queue,
    lbm_config: LbmConfig,
    ocl_code: String,

    kernel_initialize: Kernel, // Basic Kernels
    kernel_stream_collide: Kernel,
    kernel_update_fields: Kernel,

    pub n_x: u32, // Domain size
    pub n_y: u32,
    pub n_z: u32,

    d_x: u32, // Domain
    d_y: u32,
    d_z: u32,

    o_x: i32, // Offset
    o_y: i32,
    o_z: i32,

    nu: f32, // Kinematic shear viscosity
    fx: f32, // Volume force
    fy: f32,
    fz: f32,
    sigma: f32, // surface tension coefficient
    alpha: f32,
    beta: f32,

    particles_n: u32,   // Number of particles
    particles_rho: f32, // Particle mass

    pub fi32: Buffer<f32>, // Buffers
    pub fi16: Buffer<u16>,
    pub rho: Buffer<f32>,
    pub u: Buffer<f32>,
    pub q: Buffer<f32>,
    pub flags: Buffer<u8>,
    pub t: u64, // Timestep

    pub graphics: Graphics, // Graphics struct, handles rendering
}

impl LbmDomain {
    pub fn init(
        lbm_config: LbmConfig,
        device: Device,
        h_x: u32,
        h_y: u32,
        h_z: u32,
        x: u32,
        y: u32,
        z: u32,
    ) -> LbmDomain {
        let n_x = lbm_config.n_x / lbm_config.d_x + 2u32 * h_x; // Size + Halo offsets
        let n_y = lbm_config.n_y / lbm_config.d_y + 2u32 * h_y; // When multiple domains on axis -> add 2 cells of padding
        let n_z = lbm_config.n_z / lbm_config.d_z + 2u32 * h_z;
        let n = Self::get_n(n_x, n_y, n_z);

        let d_x = lbm_config.d_x;
        let d_y = lbm_config.d_y;
        let d_z = lbm_config.d_z;

        let o_x = (x * n_x / d_x) as i32 - h_x as i32;
        let o_y = (y * n_y / d_y) as i32 - h_y as i32;
        let o_z = (z * n_z / d_z) as i32 - h_z as i32;

        let dimensions;
        let velocity_set;
        let transfers;

        match lbm_config.velocity_set {
            //Set dimensions/velocitys/transfers from Enum
            VelocitySet::D2Q9 => {
                dimensions = 2;
                velocity_set = 9;
                transfers = 3;
            }
            VelocitySet::D3Q15 => {
                dimensions = 3;
                velocity_set = 15;
                transfers = 5;
            }
            VelocitySet::D3Q19 => {
                dimensions = 3;
                velocity_set = 19;
                transfers = 5;
            }
            VelocitySet::D3Q27 => {
                dimensions = 3;
                velocity_set = 27;
                transfers = 9;
            }
        }

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
        ) + &if lbm_config.graphics_config.graphics {
            graphics::get_graphics_defines(lbm_config.graphics_config)
        } else {
            "".to_string()
        } + &opencl::get_opencl_code(); // Only appends graphics defines if needed

        // OCL variables are directly exposed, due to no other device struct.
        let platform = Platform::default();
        let context = Context::builder()
            .platform(platform)
            .devices(device.clone())
            .build()
            .unwrap();
        let queue = Queue::new(&context, device, None).unwrap();
        println!("Compiling Program...");
        let program = Program::builder()
            .devices(device)
            .src(&ocl_code)
            .build(&context)
            .unwrap();

        // Initialize Buffers
        let fi32: Buffer<f32>;
        let fi16: Buffer<u16>;
        match lbm_config.float_type {
            // Initialise two fiBuffer variants for rust. Only one will be used.
            FloatType::FP32 => {
                // Float Type F32
                fi32 = Buffer::<f32>::builder()
                    .queue(queue.clone())
                    .len([n * velocity_set as u64])
                    .fill_val(0.0f32)
                    .build()
                    .unwrap();
                //Evil memory hack. fi16 won't be used, but needs to be allocated for rust :)
            }
            _ => {
                // Float Type F16S/F16C
                fi16 = Buffer::<u16>::builder()
                    .queue(queue.clone())
                    .len([n * velocity_set as u64])
                    .fill_val(0u16)
                    .build()
                    .unwrap();
                //Evil memory hack. fi32 won't be used, but needs to be allocated for rust :)
            }
        };
        let rho = Buffer::<f32>::builder()
            .queue(queue.clone())
            .len([n])
            .fill_val(1.0f32)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
        let u = Buffer::<f32>::builder()
            .queue(queue.clone())
            .len([n * 3])
            .fill_val(0.0f32)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
        let flags = Buffer::<u8>::builder()
            .queue(queue.clone())
            .len([n])
            .fill_val(0u8)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();

        // Force field buffer
        let f: Buffer<f32>;
        if lbm_config.ext_force_field {
            f = Buffer::<f32>::builder()
                .queue(queue.clone())
                .len([n * 3])
                .fill_val(0.0f32)
                .flags(flags::MEM_READ_WRITE)
                .build()
                .unwrap();
        }
        // Electric charge buffer.
        let q: Buffer<f32>;
        if lbm_config.ext_electric_force {
            q = Buffer::<f32>::builder()
                .queue(queue.clone())
                .len([n])
                .fill_val(0.0f32)
                .flags(flags::MEM_READ_WRITE)
                .build()
                .unwrap();
        }

        // Initialize Kernels
        let kernel_initialize: Kernel;
        let kernel_stream_collide: Kernel;
        let kernel_update_fields: Kernel;
        match lbm_config.float_type {
            //Initialize kernels. Different Float types need different arguments (Fi-Buffer specifically)
            FloatType::FP32 => {
                // Float Type F32
                kernel_initialize = Kernel::builder()
                    .program(&program)
                    .name("initialize")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi32)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("q", &q)
                    .arg_named("flags", &flags)
                    .build()
                    .unwrap();
                kernel_stream_collide = Kernel::builder()
                    .program(&program)
                    .name("stream_collide")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi32)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("q", &q)
                    .arg_named("flags", &flags)
                    .arg_named("t", &t)
                    .arg_named("fx", &lbm_config.fx)
                    .arg_named("fy", &lbm_config.fy)
                    .arg_named("fz", &lbm_config.fz)
                    .build()
                    .unwrap();
                kernel_update_fields = Kernel::builder()
                    .program(&program)
                    .name("update_fields")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi32)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("flags", &flags)
                    .arg_named("t", &t)
                    .arg_named("fx", &lbm_config.fx)
                    .arg_named("fy", &lbm_config.fy)
                    .arg_named("fz", &lbm_config.fz)
                    .build()
                    .unwrap();
            }
            _ => {
                // Float Type F16S/F16C
                kernel_initialize = Kernel::builder()
                    .program(&program)
                    .name("initialize")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi16)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("q", &q)
                    .arg_named("flags", &flags)
                    .build()
                    .unwrap();
                kernel_stream_collide = Kernel::builder()
                    .program(&program)
                    .name("stream_collide")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi16)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("q", &q)
                    .arg_named("flags", &flags)
                    .arg_named("t", &t)
                    .arg_named("fx", &lbm_config.fx)
                    .arg_named("fy", &lbm_config.fy)
                    .arg_named("fz", &lbm_config.fz)
                    .build()
                    .unwrap();
                kernel_update_fields = Kernel::builder()
                    .program(&program)
                    .name("update_fields")
                    .queue(queue.clone())
                    .global_work_size([n])
                    .arg_named("fi", &fi16)
                    .arg_named("rho", &rho)
                    .arg_named("u", &u)
                    .arg_named("flags", &flags)
                    .arg_named("t", &t)
                    .arg_named("fx", &lbm_config.fx)
                    .arg_named("fy", &lbm_config.fy)
                    .arg_named("fz", &lbm_config.fz)
                    .build()
                    .unwrap();
            }
        }
        println!("Kernels for domain compiled.");
        //TODO: allocate transfer buffers

        let graphics = Graphics::new(lbm_config.clone(), &program, &queue, &flags, &u);

        let domain = LbmDomain {
            device,
            context,
            queue,
            program,
            ocl_code,
            lbm_config,

            kernel_initialize,
            kernel_stream_collide,
            kernel_update_fields,

            n_x,
            n_y,
            n_z,

            d_x,
            d_y,
            d_z,

            o_x,
            o_y,
            o_z,

            nu: lbm_config.nu,
            fx: lbm_config.fx,
            fy: lbm_config.fy,
            fz: lbm_config.fz,
            sigma: lbm_config.sigma,
            alpha: lbm_config.alpha,
            beta: lbm_config.beta,

            particles_n: lbm_config.particles_n,
            particles_rho: lbm_config.particles_rho,

            fi32,
            fi16,
            rho,
            u,
            q,
            flags,
            t,

            graphics,
        };
        domain //Returns initialised domain
    }

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
        lbm_config: LbmConfig,
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

        return "\n    #define def_Nx ".to_owned() + &n_x.to_string()+"u"
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

        +"\n	#define def_domain_offset_x "+ &format!("{:.5}", (o_x as f32+(d_x>1) as i32 as f32 - 0.5*(d_x as f32 - 1.0) * (n_x - 2u32 * (d_x>1u32) as i32 as u32) as f32))+"f"
        +"\n	#define def_domain_offset_y "+ &format!("{:.5}", (o_y as f32+(d_y>1) as i32 as f32 - 0.5*(d_y as f32 - 1.0) * (n_y - 2u32 * (d_y>1u32) as i32 as u32) as f32))+"f"
        +"\n	#define def_domain_offset_z "+ &format!("{:.5}", (o_z as f32+(d_z>1) as i32 as f32 - 0.5*(d_z as f32 - 1.0) * (n_z - 2u32 * (d_z>1u32) as i32 as u32) as f32))+"f"

        +"\n	#define D"+ &dimensions.to_string()+"Q"+ &velocity_set.to_string()+"" // D2Q9/D3Q15/D3Q19/D3Q27
        +"\n	#define def_velocity_set "+ &velocity_set.to_string()+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
        +"\n	#define def_dimensions "+ &dimensions.to_string()+"u" // number spatial dimensions (2D or 3D)
        +"\n	#define def_transfers "+ &transfers.to_string()+"u" // number of DDFs that are transferred between multiple domains

        +"\n	#define def_c 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
        +"\n	#define def_w " + &format!("{:.5}", 1.0f32/(3.0f32*nu+0.5f32))+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
        + match lbm_config.velocity_set {
            VelocitySet::D2Q9 => &d2q9,
            VelocitySet::D3Q15 => &d3q15,
            VelocitySet::D3Q19 => &d3q19,
            VelocitySet::D3Q27 => &d3q27
        }
        + match lbm_config.relaxation_time {
            RelaxationTime::SRT => srt,
            RelaxationTime::TRT => trt
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
        + if lbm_config.ext_electric_force {"\n	        #define ELECTRIC_FORCE"} else {""}
        + if lbm_config.ext_force_field {"\n	        #define FORCE_FIELD"} else {""}
        + if lbm_config.graphics_config.graphics {"\n	#define UPDATE_FIELDS"} else {""};
        //Extensions
    }

    fn enqueue_initialize(&self) -> ocl::Result<()> {
        //Enqueues Initialization kernel, arguments are already set
        unsafe { self.kernel_initialize.cmd().enq()? }
        self.queue.finish()
    }

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

    // Manually update fields. Is autmatically handled in stream-collide if graphics are enabled
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

    fn increment_timestep(&mut self, steps: u32) {
        self.t += steps as u64;
    }

    fn reset_timestep(&mut self) {
        self.t = 0;
    }

    fn get_n(n_x: u32, n_y: u32, n_z: u32) -> u64 {
        n_x as u64 * n_y as u64 * n_z as u64
    }

    fn self_get_n(&self) -> u64 {
        self.lbm_config.n_x as u64 * self.lbm_config.n_y as u64 * self.lbm_config.n_z as u64
    }

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
