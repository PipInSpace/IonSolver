//! # domain
//! 
//! Defines the `LbmDomain` struct that holds all information to run a LBM-Simulation on one OpenCL Device.
//! 
//! `LbmDomain` should not be initialized on it's own, but automatically through the `Lbm::new()` function when initializing a new [`Lbm`].
//! This ensures all arguments are correctly set.

use ocl::{Buffer, Context, Device, Kernel, Platform, Program, Queue};
use std::cmp;
use crate::lbm::*;
use ocl_macros::*;

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
    pub fqi: Option<VariableFloatBuffer>,
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
        let n = n_x as u64 * n_y as u64 * n_z as u64;

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

        let ocl_code = get_device_defines(n_x, n_y, n_z, d_x, d_y, d_z, i, o_x, o_y, o_z, dimensions, velocity_set, transfers, lbm_config.nu, n_lod, n_lod_own, lbm_config)
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
        let fqi: Option<VariableFloatBuffer> = if lbm_config.ext_magneto_hydro {
            Some(match lbm_config.float_type {
                FloatType::FP32 => { VariableFloatBuffer::F32(buffer!(&queue, [n * 7], 0.0f32)) } // Float Type F32
                _               => { VariableFloatBuffer::U16(buffer!(&queue, [n * 7], 0u16)) }   // Float Type F16S/F16C
            })
        } else { None };
        // Electron gas ddfs
        let ei: Option<VariableFloatBuffer> = if lbm_config.ext_magneto_hydro {
            Some(match lbm_config.float_type {
                FloatType::FP32 => { VariableFloatBuffer::F32(buffer!(&queue, [n * velocity_set as u64], 0.0f32)) } // Float Type F32
                _               => { VariableFloatBuffer::U16(buffer!(&queue, [n * velocity_set as u64], 0u16)) }   // Float Type F16S/F16C
            })
        } else { None };
        let q: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n], 0f32)) } else { None };
        // Level of detail charge and velocity for field updates
        let qu_lod: Option<Buffer<f32>> = if lbm_config.ext_magneto_hydro { Some(buffer!(&queue, [n_lod * 4], 0f32)) } else { None };
        // Transfer buffer for LODs
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
            
            match fqi.as_ref().expect("fqi should be initialized") {
                VariableFloatBuffer::U16(fqi_u16) => {kernel_args!(stream_collide_builder, ("fqi", fqi_u16)); kernel_args!(initialize_builder, ("fqi", fqi_u16));},
                VariableFloatBuffer::F32(fqi_f32) => {kernel_args!(stream_collide_builder, ("fqi", fqi_f32)); kernel_args!(initialize_builder, ("fqi", fqi_f32));},
            };
            match ei.as_ref().expect("ei should be initialized") {
                VariableFloatBuffer::U16(ei_u16) => {kernel_args!(stream_collide_builder, ("ei", ei_u16)); kernel_args!(initialize_builder, ("ei", ei_u16));},
                VariableFloatBuffer::F32(ei_f32) => {kernel_args!(stream_collide_builder, ("ei", ei_f32)); kernel_args!(initialize_builder, ("ei", ei_f32));},
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
                    Some(match &fqi.as_ref().expect("fqi should be defined") {
                        VariableFloatBuffer::U16(fqi_u16) => kernel!(program, queue, "transfer_extract_fqi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fqi", fqi_u16)),
                        VariableFloatBuffer::F32(fqi_f32) => kernel!(program, queue, "transfer_extract_fqi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fqi", fqi_f32)),
                    })
                } else { None } , // Extract Qi
                if lbm_config.ext_magneto_hydro {
                    Some(match &fqi.as_ref().expect("fqi should be defined") {
                        VariableFloatBuffer::U16(fqi_u16) => kernel!(program, queue, "transfer__insert_fqi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fqi", fqi_u16)),
                        VariableFloatBuffer::F32(fqi_f32) => kernel!(program, queue, "transfer__insert_fqi", 1, ("direction", 0_u32), ("time_step", 0_u64), ("transfer_buffer_p", &transfer_p), ("transfer_buffer_m", &transfer_m), ("fqi", fqi_f32)),
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

            e, b, e_dyn, b_dyn, fqi, q, qu_lod, f, t, graphics, cfg: lbm_config.clone()
        } //Returns initialised domain
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

    #[allow(dead_code)]
    pub fn is_halo(&mut self, x: u32, y: u32, z: u32) -> bool {
	    return ((self.cfg.d_x>1)&(x==0||x>=self.cfg.n_x-1))||((self.cfg.d_y>1)&(y==0||y>=self.cfg.n_y-1))||((self.cfg.d_z>1)&(z==0||z>=self.cfg.n_z-1));
    }

    /// Get `x, y, z` coordinates from 1D index `n`.
    pub fn get_coordinates(&self, n: u64) -> (u32, u32, u32) {
        get_coordinates_sl(n, self.n_x, self.n_y)
    }
}