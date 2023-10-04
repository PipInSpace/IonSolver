extern crate ocl;
use crate::lbm::*;
use crate::*;
use image::{ImageBuffer, Rgb};
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
    lbm_config.n_x = 512;
    lbm_config.n_y = 512;
    lbm_config.n_z = 256;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    lbm_config.ext_equilibrium_boudaries = true;
    let mut test_lbm = Lbm::init(lbm_config);
    test_lbm.initialize();

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
