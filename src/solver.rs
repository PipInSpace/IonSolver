extern crate ocl;
use crate::lbm::*;
use crate::*;
use lbm::VelocitySet;
use ocl::ProQue;
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
        paused: false,
        save: false,
        active: true,
    };
    let mut i = 0;

    //OpenCL Test function
    test_function().unwrap();

    let lbm_config = LbmConfig::new();
    let test_lbm = Lbm::init(lbm_config);

    //OpenCL setup
    let src = include_str!("kernels.cl");

    let pro_que = ProQue::builder()
        .src(src)
        .dims([1 << 6, 1 << 6])
        .build()
        .unwrap();

    let buffer_a = pro_que.create_buffer::<f32>().unwrap();
    let buffer_b = pro_que.create_buffer::<f32>().unwrap();

    let setup_kernel = pro_que
        .kernel_builder("add")
        .arg(&buffer_a)
        .arg(1.0f32)
        .build()
        .unwrap();

    let gauss_seidel_step_kernel = pro_que
        .kernel_builder("gauss_seidel_step")
        .arg(&buffer_a)
        .arg(&buffer_b)
        .arg(4.096f32)
        .arg(17.384f32)
        .build()
        .unwrap();

    //Execute setup kernels
    unsafe {
        setup_kernel.enq().unwrap();
        setup_kernel.set_arg(0, &buffer_b).unwrap();
        setup_kernel.enq().unwrap();
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
                unsafe {
                    for _i in 1..20 {
                        gauss_seidel_step_kernel.enq().unwrap();
                    }
                }

                if i % 1000 == 0 {
                    let mut vec = vec![0.0f32; buffer_a.len()];
                    buffer_a.read(&mut vec).enq().unwrap();

                    println!(
                        "Step {}: The value at index [{}] is now '{}'!",
                        i, 130, vec[130]
                    );
                }

                i += 1;
            }
        }
        if !state.active {
            break;
        }
        let recieve_result = ctrl_rx.recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        }
    }
}

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
