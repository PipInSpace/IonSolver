extern crate ocl;
use crate::lbm::*;
use crate::*;
use ocl::ProQue;
use std::f32::consts::PI;
use std::fs;
use std::sync::mpsc;

#[allow(unused)]
pub fn simloop(
    sim: SimState,
    sim_tx: mpsc::Sender<SimState>,
    ctrl_rx: mpsc::Receiver<SimControlTx>,
) {
    //TODO: Simulation here
    //sim_tx.send(sim) sends data to the main window loop
    let mut state = SimControlTx {
        paused: true,
        save: false,
        clear_images: true,
        frame_spacing: 10,
        active: true,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    };
    let mut i = 0;

    let mut lbm_config = LbmConfig::new();
    lbm_config.n_x = 256;
    lbm_config.n_y = 256;
    lbm_config.n_z = 128;
    lbm_config.d_x = 2;
    lbm_config.nu = 0.01;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    let mut test_lbm = Lbm::init(lbm_config);
    test_lbm.setup();
    test_lbm.domains[0].graphics.streamline_mode = true;
    test_lbm.domains[0].graphics.q_mode = false;
    test_lbm.initialize();

    // get initial config from ui
    let recieve_result = ctrl_rx.try_recv();
    if let Ok(recieve) = recieve_result {
        state = recieve;
    }

    // Clearing out folder if requested
    if state.save && state.clear_images {
        match fs::remove_dir_all("out") {
            Ok(_) => (),
            Err(_) => println!("Did not find out folder. Creating it."),
        }
        fs::create_dir("out").unwrap();
    }

    let mut has_commenced = false; //has the simulation started
    let mut cached_rot: Vec<f32> = vec![0.0; 2];
    let mut cached_zoom = 0.0;
    loop {
        //This is the master loop, cannot be paused
        if !state.paused {
            has_commenced = true;
            loop {
                //This is the loop of the simulation. Can be paused by receiving a control message
                let recieve_result = ctrl_rx.try_iter().last();
                if let Some(recieve) = recieve_result {
                    state = recieve;
                }
                if state.paused || !state.active {
                    break;
                }

                //Simulation commences here
                test_lbm.do_time_step();
                if cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom
                {
                    //only update params when camera updated
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &test_lbm.domains {
                        d.graphics.camera_params.write(&params).enq().unwrap();
                    }
                }

                if i % state.frame_spacing == 0 {
                    println!("Step {}", i);
                    if lbm_config.graphics_config.graphics {
                        test_lbm.draw_frame(state.save, state.frame_spacing, sim_tx.clone(), i);
                    }
                }
                i += 1;
            }
        }
        if state.paused && state.active {
            let mut cached_rot: Vec<f32> = vec![0.0; 2];
            let mut cached_zoom = 0.0;
            loop {
                //This is the loop of the simulation if it is paused but active. useful for displaying the simulation
                let recieve_result = ctrl_rx.try_iter().last();
                if let Some(recieve) = recieve_result {
                    state = recieve;
                }
                if !state.paused || !state.active {
                    break;
                }

                if cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom
                    || !has_commenced
                {
                    //only draw when camera updated
                    has_commenced = true;
                    cached_rot = state.camera_rotation.clone();
                    cached_zoom = state.camera_zoom;
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &test_lbm.domains {
                        d.graphics.camera_params.write(&params).enq().unwrap();
                    }
                    if lbm_config.graphics_config.graphics {
                        test_lbm.draw_frame(state.save, state.frame_spacing, sim_tx.clone(), i);
                    }
                    thread::sleep(Duration::from_millis(33)) // about 30 FPS
                }
            }
        }
        if !state.active {
            println!("Exiting Simulation Loop");
            break;
        }
        let recieve_result = ctrl_rx.try_recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        }
    }
}

impl Lbm {
    pub fn setup(&mut self) {
        // 3D Taylor-Green vortices
        println!("Setting up Taylor-Green vorticies");
        let nx = self.config.n_x;
        let ny = self.config.n_y;
        let nz = self.config.n_z;
        let ntotal = nx as u64 * ny as u64 * nz as u64;
        let pif = std::f32::consts::PI;
        #[allow(non_snake_case)]
        let A = 0.25f32;
        let periodicity = 1u32;
        let a = nx as f32 / periodicity as f32;
        let b = ny as f32 / periodicity as f32;
        let c = nz as f32 / periodicity as f32;
        let mut vec_u: Vec<f32> = vec![0.0; (nx as u64 * ny as u64 * nz as u64 * 3) as usize];
        let mut vec_rho: Vec<f32> = vec![0.0; (nx as u64 * ny as u64 * nz as u64) as usize];
        for n in 0..(nx as u64 * ny as u64 * nz as u64) {
            let (x, y, z) = self.get_coordinates(n);
            let fx = x as f32 + 0.5 - 0.5 * nx as f32;
            let fy = y as f32 + 0.5 - 0.5 * ny as f32;
            let fz = z as f32 + 0.5 - 0.5 * nz as f32;
            vec_u[(n) as usize] = A
                * (2.0 * pif * fx / a).cos()
                * (2.0 * pif * fy / b).sin()
                * (2.0 * pif * fz / c).sin(); // x
            vec_u[(n + ntotal) as usize] = -A
                * (2.0 * pif * fx / a).sin()
                * (2.0 * pif * fy / b).cos()
                * (2.0 * pif * fz / c).sin(); // y
            vec_u[(n + (ntotal * 2)) as usize] = A
                * (2.0 * pif * fx / a).sin()
                * (2.0 * pif * fy / b).sin()
                * (2.0 * pif * fz / c).cos(); // z
            vec_rho[n as usize] =
                1.0 - (A * A) * 3.0 / 4.0 * (4.0 * pif * fx / a).cos() + (4.0 * pif * fy / b).cos();
        }
        let domain_numbers: u32 = self.config.d_x * self.config.d_y * self.config.d_z;
        for d in 0..domain_numbers {
            let dx = self.config.d_x;
            let dy = self.config.d_y;
            let dz = self.config.d_z;
            println!("Initializing domain {}/{}", d + 1, domain_numbers);
            let x = (d % (dx * dy)) % dx; // Domain coordinates
            let y = (d % (dx * dy)) / dx;
            let z = d / (dx * dy);
            let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64*2; // Domain size on each axis
            let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64*2; // Needs to account for halo offsets
            let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64*2;
            let dtotal = dsx * dsy * dsz;

            let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
            let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
            for zi in 0..dsz as u64 {
                // iterates over every cell in the domain, loading the information from the precomputed all-domain-vector
                for yi in 0..dsy as u64 {
                    for xi in 0..dsx as u64 {
                        if !(((xi==0||xi==dsx-1) && dx > 1)||((yi==0||yi==dsy-1) && dy > 1)||((zi==0||zi==dsz-1) && dz > 1)){//do not set at halo offsets
                            let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                            let gx = xi - (dx > 1u32) as u64 + x as u64 * (dsx - (dx > 1u32) as u64*2);
                            let gy = yi - (dy > 1u32) as u64 + y as u64 * (dsy - (dy > 1u32) as u64*2);
                            let gz = zi - (dz > 1u32) as u64 + z as u64 * (dsz - (dz > 1u32) as u64*2);
                            let gn = gx + (gy * nx as u64) + (gz * nx as u64 * ny as u64);
                            domain_vec_u[(dn) as usize] = vec_u[(gn) as usize];
                            domain_vec_u[(dn + dtotal) as usize] = vec_u[(gn + ntotal) as usize];
                            domain_vec_u[(dn + dtotal * 2) as usize] = vec_u[(gn + ntotal * 2) as usize];
                            domain_vec_rho[(dn) as usize] = vec_rho[(gn) as usize];
                        }
                    }
                }
            }

            self.domains[d as usize]
                .u
                .write(&domain_vec_u)
                .enq()
                .unwrap();
            self.domains[d as usize]
                .rho
                .write(&domain_vec_rho)
                .enq()
                .unwrap();
            self.domains[d as usize].queue.finish().unwrap();
        }
        //self.u.write(&vec_u).enq().unwrap();
        //self.rho.write(&vec_rho).enq().unwrap();
        //self.queue.finish().unwrap();
        println!("Finished setting up Taylor-Green vorticies");
    }
}

#[allow(unused)]
fn test_function() -> ocl::Result<()> {
    let src = include_str!("kernels.cl");

    let pro_que = ProQue::builder().src(src).dims([1 << 6, 1 << 6]).build()?;

    let buffer_a = pro_que.create_buffer::<f32>()?;
    let buffer_b = pro_que.create_buffer::<f32>()?;

    let setup_kernel = pro_que
        .kernel_builder("add")
        .arg(&buffer_a)
        .arg(1.0f32)
        .build()?;

    let gauss_seidel_step_kernel = pro_que
        .kernel_builder("gauss_seidel_step")
        .arg(&buffer_a)
        .arg(&buffer_b)
        .arg(4.096f32)
        .arg(17.384f32)
        .build()?;

    unsafe {
        setup_kernel.enq()?;
        setup_kernel.set_arg(0, &buffer_b)?;
        setup_kernel.enq()?;
        for _i in 1..20 {
            gauss_seidel_step_kernel.enq()?;
        }
    }

    let mut vec = vec![0.0f32; buffer_a.len()];
    buffer_a.read(&mut vec).enq()?;

    println!("Gauss: The value at index [{}] is now '{}'!", 131, vec[131]);
    Ok(())
}
