extern crate ocl;
use crate::lbm::*;
use crate::*;
use ocl::ProQue;
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
        save: true,
        clear_images: true,
        frame_spacing: 10,
        active: true,
    };
    let mut i = 0;

    //OpenCL Test function
    //test_function().unwrap();

    let mut lbm_config = LbmConfig::new();
    lbm_config.n_x = 256;
    lbm_config.n_y = 256;
    lbm_config.n_z = 256;
    lbm_config.nu = 0.01;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    let mut test_lbm = Lbm::init(lbm_config);
    test_lbm.domains[0].setup();
    test_lbm.initialize();
    let mut test_vec:Vec<f32> = vec![0.0; 256*256*256];
    test_lbm.domains[0].rho.read(&mut test_vec).enq().unwrap();
    println!("rho at index 5000: {}", test_vec[5000]);

    // get initial config from ui
    let recieve_result = ctrl_rx.try_recv();
    if let Ok(recieve) = recieve_result {
        state = recieve;
    }

    // Clearing out folder if requested
    if state.save && state.clear_images {
        fs::remove_dir_all("out").unwrap(); // TODO: make not bad
        fs::create_dir("out").unwrap();
    }

    loop {
        //This is the master loop, cannot be paused
        if !state.paused {
            loop {
                //This is the loop of the simulation. Can be paused by receiving a control message
                let recieve_result = ctrl_rx.try_recv();
                if let Ok(recieve) = recieve_result {
                    state = recieve;
                }
                if state.paused || !state.active {
                    break;
                }

                //Simulation commences here
                test_lbm.do_time_step();
                test_lbm.domains[0].rho.read(&mut test_vec).enq().unwrap();
                println!("rho at index 5000: {}", test_vec[5000]);

                if i % state.frame_spacing == 0 {
                    println!("Step {}", i);
                    if lbm_config.graphics_config.graphics {
                        test_lbm.draw_frame(state.save, state.frame_spacing, sim_tx.clone(), i);
                    }
                }
                i += 1;
            }
        }
        if !state.active {
            println!("Exiting Simulation Loop");
            break;
        }
        let recieve_result = ctrl_rx.recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        }
    }
}

impl LbmDomain {
    pub fn setup(&mut self) {
        // 3D Taylor-Green vortices
        println!("Setting up Taylor-Green vorticies");
        let nx = self.n_x;
        let ny = self.n_y;
        let nz = self.n_z;
        let pif = std::f32::consts::PI;
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
            vec_u[(n * 3) as usize] = A
                * (2.0 * pif * fx / a).cos()
                * (2.0 * pif * fy / b).sin()
                * (2.0 * pif * fz / c).sin(); // x
            vec_u[(n * 3 + 1) as usize] = A
                * (2.0 * pif * fx / a).sin()
                * (2.0 * pif * fy / b).cos()
                * (2.0 * pif * fz / c).sin(); // y
            vec_u[(n * 3 + 2) as usize] = A
                * (2.0 * pif * fx / a).sin()
                * (2.0 * pif * fy / b).sin()
                * (2.0 * pif * fz / c).cos(); // z
            vec_rho[n as usize] =
                1.0 - (A * A) * 3.0 / 4.0 * (4.0 * pif * fx / a).cos() + (4.0 * pif * fy / b).cos();
        }
        self.u.write(&vec_u).enq().unwrap();
        self.rho.write(&vec_rho).enq().unwrap();
        self.queue.finish().unwrap();
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
