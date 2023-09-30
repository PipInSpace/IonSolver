use crate::*;
use ocl::Device;

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum VelocitySet {
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27,
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum RelaxationTime {
    SRT,
    TRT,
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum FloatType {
    FP16S,
    FP16C,
    FP32,
}

// ############################################# Code Structure #############################################
// To start a LBM simulation, one needs to initialise an Lbm struct using the Lbm::init() function.
// The init() function takes in another struct, the LbmConfig, which contains all necessary arguments.
// LbmConfig needs to be configured beforehand, then the init() function takes care of Lbm initialisation.
// The Lbm struct contains one or multiple LbmDomain structs.
// The Simulation is actually run on these LbmDomain structs, each Domain corresponding to an OpenCL device.
// The simulation can be split up into multiple domains. This enables multi-device parallelization.
// LbmDomain initialisation is handled automatically in Lbm::init() using LbmDomain::init()

#[derive(Clone, Copy)]
pub struct LbmConfig {
    velocity_set: VelocitySet,
    relaxation_time: RelaxationTime,
    float_type: FloatType,
    n_x: u32, //Size
    n_y: u32,
    n_z: u32,

    d_x: u32, //Domains
    d_y: u32,
    d_z: u32,

    nu: f32,
    fx: f32,
    fy: f32,
    fz: f32,
    sigma: f32,
    alpha: f32,
    beta: f32,

    particles_n: u32,
    particles_rho: f32,

    ext_equilibrium_boudaries: bool, //Extensions
}

impl LbmConfig {
    pub fn new() -> LbmConfig {
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
        }
    }
}

#[derive(Clone)]
pub struct Lbm {
    domains: Vec<LbmDomain>,
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
            let x = (d % (lbm_config.d_x * lbm_config.d_y)) % lbm_config.d_x;
            let y = (d % (lbm_config.d_x * lbm_config.d_y)) / lbm_config.d_x;
            let z = d / (lbm_config.d_x * lbm_config.d_y);
            //Get devices
            lbm_domains.push(LbmDomain::init(lbm_config, device_infos[d as usize], h_x, h_y, h_z, x, y, z))
        }

        let lbm = Lbm {
            domains: vec![LbmDomain::from_velocity_set(lbm_config.velocity_set)],
            config: lbm_config,
            initialized: false,
        };
        lbm
    }

    fn initialize(self) {
        //for d in 0..self.get_domain_numbers() {self.domains[d]}
    }

    pub fn run(self, steps: u64) {
        //Run simulation for steps
        //TODO: Display info in command line
        if !self.initialized {
            //Run initialization Kernel
            self.initialize();
            //TODO: Display initialization info in comman line
        }
    }

    //Helper functions:
}

#[derive(Clone)]
pub struct LbmDomain {
    device_info: Device, //TODO: Compiles OpenCL C code
    lbm_config: LbmConfig,
    ocl_code: String,

    n_x: u32, //Size
    n_y: u32,
    n_z: u32,

    d_x: u32, //Domain
    d_y: u32,
    d_z: u32,

    o_x: i32, //Offset
    o_y: i32,
    o_z: i32,

    nu: f32,
    fx: f32,
    fy: f32,
    fz: f32,
    sigma: f32,
    alpha: f32,
    beta: f32,

    particles_n: u32,
    particles_rho: f32,

    velocity_set_enum: VelocitySet,
    dimensions: u8,
    velocity_set: u8,
    transfers: u8,
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
        let lbm_domain_velocity_set = LbmDomain::from_velocity_set(lbm_config.velocity_set);

        let mut domain = LbmDomain {
            device_info: device,
            ocl_code: String::new(),
            lbm_config,

            n_x: lbm_config.n_x / lbm_config.d_x + 2u32 * h_x,
            n_y: lbm_config.n_y / lbm_config.d_y + 2u32 * h_y,
            n_z: lbm_config.n_z / lbm_config.d_z + 2u32 * h_z,

            d_x: lbm_config.d_x,
            d_y: lbm_config.d_y,
            d_z: lbm_config.d_z,

            o_x: (x * lbm_config.n_x / lbm_config.d_x) as i32 - h_x as i32,
            o_y: (y * lbm_config.n_y / lbm_config.d_y) as i32 - h_y as i32,
            o_z: (z * lbm_config.n_z / lbm_config.d_z) as i32 - h_z as i32,

            nu: lbm_config.nu,
            fx: lbm_config.fx,
            fy: lbm_config.fy,
            fz: lbm_config.fz,
            sigma: lbm_config.sigma,
            alpha: lbm_config.alpha,
            beta: lbm_config.beta,

            particles_n: lbm_config.particles_n,
            particles_rho: lbm_config.particles_rho,

            velocity_set_enum: lbm_config.velocity_set,
            dimensions: lbm_domain_velocity_set.dimensions,
            velocity_set: lbm_domain_velocity_set.velocity_set,
            transfers: lbm_domain_velocity_set.transfers, //Completes ONLY velocity set values, check for other completions
        };
        let ocl_code = domain.clone().get_device_defines() + &LbmDomain::get_opencl_code();
        domain.ocl_code = ocl_code;
        domain //Returns initialised domain
    }

    fn get_default_domain() -> LbmDomain {
        LbmDomain {
            device_info: opencl::get_devices()[0],
            ocl_code: String::new(),
            lbm_config: LbmConfig::new(),
            n_x: 1,
            n_y: 1,
            n_z: 1,
            d_x: 1,
            d_y: 1,
            d_z: 1,
            o_x: 0,
            o_y: 0,
            o_z: 0,
            fx: 0.0f32,
            fy: 0.0f32,
            fz: 0.0f32,
            nu: 0.0f32,
            sigma: 0.0f32,
            alpha: 0.0f32,
            beta: 0.0f32,

            particles_n: 0,
            particles_rho: 0.0f32,

            velocity_set_enum: VelocitySet::D2Q9,
            dimensions: 0,
            velocity_set: 0,
            transfers: 0,
        }
    }

    fn from_velocity_set(vel_set_cfg: VelocitySet) -> LbmDomain {
        let default_domain = LbmDomain::get_default_domain();

        match vel_set_cfg {
            VelocitySet::D2Q9 => LbmDomain {
                dimensions: 2,
                velocity_set: 9,
                transfers: 3,
                ..default_domain
            },
            VelocitySet::D3Q15 => LbmDomain {
                dimensions: 3,
                velocity_set: 15,
                transfers: 5,
                ..default_domain
            },
            VelocitySet::D3Q19 => LbmDomain {
                dimensions: 3,
                velocity_set: 19,
                transfers: 5,
                ..default_domain
            },
            VelocitySet::D3Q27 => LbmDomain {
                dimensions: 3,
                velocity_set: 27,
                transfers: 9,
                ..default_domain
            },
        }
    }

    fn get_opencl_code() -> String {
        return include_str!("kernels.cl").to_string(); //TODO Kernel needs to be preproccessed first
    }

    fn get_device_defines(self) -> String {
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

        return "\n    #define def_Nx ".to_owned() + &self.n_x.to_string()+"u"
        +"\n	#define def_Ny "+ &self.n_y.to_string()+"u"
        +"\n	#define def_Nz "+ &self.n_z.to_string()+"u"
        +"\n	#define def_N  "+ &self.clone().get_n().to_string()+"ul"
    
        +"\n	#define def_Dx "+ &self.d_x.to_string()+"u"
        +"\n	#define def_Dy "+ &self.d_y.to_string()+"u"
        +"\n	#define def_Dz "+ &self.d_z.to_string()+"u"
    
        +"\n	#define def_Ox "+ &self.o_x.to_string()+"" // offsets are signed integer!
        +"\n	#define def_Oy "+ &self.o_y.to_string()+""
        +"\n	#define def_Oz "+ &self.o_z.to_string()+""
    
        +"\n	#define def_Ax "+ &(self.n_y * self.n_z).to_string()+"u"
        +"\n	#define def_Ay "+ &(self.n_z * self.n_x).to_string()+"u"
        +"\n	#define def_Az "+ &(self.n_x * self.n_y).to_string()+"u"
    
        +"\n	#define D"+ &self.dimensions.to_string()+"Q"+ &self.velocity_set.to_string()+"" // D2Q9/D3Q15/D3Q19/D3Q27
        +"\n	#define def_velocity_set "+ &self.velocity_set.to_string()+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
        +"\n	#define def_dimensions "+ &self.dimensions.to_string()+"u" // number spatial dimensions (2D or 3D)
        +"\n	#define def_transfers "+ &self.transfers.to_string()+"u" // number of DDFs that are transferred between multiple domains
    
        +"\n	#define def_c 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
        +"\n	#define def_w " + &(1.0f32/(3.0f32*self.nu+0.5f32)).to_string()+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
        + match self.velocity_set_enum {
            VelocitySet::D2Q9 => &d2q9,
            VelocitySet::D3Q15 => &d3q15,
            VelocitySet::D3Q19 => &d3q19,
            VelocitySet::D3Q27 => &d3q27
        }
        + match self.lbm_config.relaxation_time {
            RelaxationTime::SRT => srt,
            RelaxationTime::TRT => trt
        }
        + match self.lbm_config.float_type { //Floatingpoint types
            FloatType::FP16S => &fp16s,
            FloatType::FP16C => &fp16c,
            FloatType::FP32 => &fp32
        }
        + if self.lbm_config.ext_equilibrium_boudaries {"\n	#define EQUILIBRIUM_BOUNDARIES"} else {""};
        //Extensions
    }

    fn get_n(self) -> u64 {
        self.n_x as u64 * self.n_y as u64 * self.n_z as u64
    }
}

fn get_tau() -> f32 {
    3.0f32 * get_nu() + 0.5f32
}
fn get_nu() -> f32 {
    1.0f32 / 6.0f32
}
