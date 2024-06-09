//! # IonSolver
//! 
//! IonSolver is a magnetohydrodynamic simulation software utilizing an extended Lattice-Boltzmann-Model.

#![allow(non_snake_case)]
extern crate ocl;
use ocl_macros::*;
use std::time::{Duration, Instant};
use std::{f32::consts::PI, fs, sync::mpsc, thread};

mod file;
mod info;
mod lbm;
mod opencl;
mod setup;
#[cfg(feature = "gui")]
mod ui;

use lbm::*;

pub struct SimState {
    pub step: u32,
    #[cfg(feature = "gui")]
    pub img: egui::ColorImage,
}

impl SimState {
    pub fn new() -> SimState {
        Self {
            step: 0,
            #[cfg(feature = "gui")]
            img: egui::ColorImage::default(),
        }
    }
}

/// SimControlTx sends information/commands to the simulation
pub struct SimControlTx {
    paused: bool,
    save_img: bool,
    save_file: bool,
    load_file: bool,
    clear_images: bool,
    frame_spacing: u32,
    file_name: String,
    active: bool,
    camera_rotation: Vec<f32>,
    camera_zoom: f32,
}

fn main() {
    #[cfg(feature = "multi-node")]
    multi_node::run_node();
    #[cfg(not(feature = "multi-node"))]
    run_gpu();
}

#[allow(dead_code)]
#[allow(unused_variables)]
#[cfg(not(feature = "multi-node"))]
/// Run IonSolver on a single compute node with multi-gpu support.
fn run_gpu() {
    //println!("{}", info::LOGO_COLOR);
    println!("IonSolver - © 2024\n");

    // Create channels
    let (sim_tx, sim_rx) = mpsc::channel();
    let (ctrl_tx, ctrl_rx) = mpsc::channel();
    let exit_ctrl_tx = ctrl_tx.clone();

    // Start Simulation loop
    let handle = thread::spawn(move || {
        simloop(sim_tx, ctrl_rx);
    });

    #[cfg(feature = "gui")]
    ui::run_ui(sim_rx, ctrl_tx);

    // Send exit control message
    _ = exit_ctrl_tx.send(SimControlTx {
        paused: true,
        save_img: false,
        save_file: false,
        load_file: false,
        clear_images: false,
        frame_spacing: 10,
        file_name: "".to_string(),
        active: false,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    });

    handle.join().unwrap();
}

/// Runs the simulation and it's control logic in a different thread, used in single-node multi-gpu mode
#[allow(dead_code)]
fn simloop(sim_tx: mpsc::Sender<SimState>, ctrl_rx: mpsc::Receiver<SimControlTx>) {
    //sim_tx.send(state) sends data to the main window loop
    let mut state = SimControlTx {
        paused: true,
        save_img: false,
        save_file: false,
        load_file: false,
        clear_images: true,
        frame_spacing: 10,
        file_name: "".to_string(),
        active: true,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    };
    let mut step_c: u32 = 0;

    let mut lbm = setup::setup();
    lbm.initialize();

    // Clearing out folder if requested
    if state.save_img && state.clear_images {
        match fs::remove_dir_all("out") {
            Ok(_) => (),
            Err(_) => println!("Did not find out folder. Creating it."),
        }
        fs::create_dir("out").unwrap();
    }

    // Create graphics thread with evil pointer hacks
    if lbm.config.graphics_config.graphics_active {
        let lbm_ptr = &lbm as *const _ as usize;
        let state_ptr = &state as *const _ as usize;
        let step_ptr = &step_c as *const _ as usize;
        let sim_tx_g = sim_tx.clone();
        thread::spawn(move || {
            let lbm = unsafe { &*(lbm_ptr as *const Lbm) };
            let state = unsafe { &*(state_ptr as *const SimControlTx) };
            let step_count = unsafe { &*(step_ptr as *const u32) };
            let mut drawn_step: u32 = 0;
            let mut cached_rot: Vec<f32> = vec![0.0; 2];
            let mut cached_zoom = 0.0;

            loop {
                let now = Instant::now();
                if !state.active {
                    println!("Exiting Graphics Loop");
                    break;
                }

                // If sim not paused, only draw frame when sim time step updated
                let frame_changed = if *step_count > drawn_step {
                    drawn_step = *step_count;
                    true
                } else {
                    false
                };

                let camera_changed = cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom;
                if camera_changed {
                    //only update camera params when camera updated
                    cached_rot = state.camera_rotation.clone();
                    cached_zoom = state.camera_zoom;
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &lbm.domains {
                        bwrite!(d.graphics.as_ref().expect("graphics").camera_params, params);
                    }
                }
                if (camera_changed || (!state.paused && frame_changed))
                    && lbm.config.graphics_config.graphics_active
                {
                    // Only draws frames, never saves them
                    lbm.draw_frame(false, "frame".to_string(), sim_tx_g.clone(), step_count);
                }
                thread::sleep(Duration::from_millis(std::cmp::max(
                    0,
                    33 - (now.elapsed().as_millis() as i64),
                ) as u64));
            }
        });
    }

    let mn = (lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64) / 1000000; // N / 10E6
    let mut t_p_s = 0; // Time per step
    loop {
        //This is the master loop, cannot be paused
        if !state.paused {
            let mut step_count_time = 1;
            let loop_time = std::time::Instant::now();
            loop {
                //This is the loop of the simulation. Can be paused by receiving a control message
                pull_state(&ctrl_rx, &mut state);

                if state.paused || !state.active {
                    info::sim_speed(step_c, t_p_s, mn);
                    break;
                }

                lbm.do_time_step();

                t_p_s = (loop_time.elapsed().as_micros() / step_count_time) as u32;
                //Calculate simulation speed
                if step_c % 200 == 0 {
                    info::sim_speed(step_c, t_p_s, mn)
                }

                // Saves frames if needed
                if step_c % state.frame_spacing == 0
                    && state.save_img
                    && lbm.config.graphics_config.graphics_active
                {
                    lbm.draw_frame(true, "frame".to_string(), sim_tx.clone(), &step_c);
                }
                // Render predefined keyframes if needed
                if lbm.config.graphics_config.graphics_active {
                    lbm.render_keyframes(sim_tx.clone(), &step_c);
                }
                step_c += 1;
                step_count_time += 1;
            }
        } else {
            if state.save_file {
                match file::write(&lbm, state.file_name.as_str()) {
                    Ok(_) => {},
                    Err(msg) => {println!("ERROR: Could not write file \"{}\"!\nReason: {}", state.file_name.as_str(), msg)},
                }
                state.save_file = false;
            }
            if state.load_file {
                // copy config to keep graphics settings
                let mut new_lbm_config: LbmConfig = LbmConfig::new();
                new_lbm_config.graphics_config.graphics_active = true;
                new_lbm_config.graphics_config.background_color = 0x1c1b22;
                new_lbm_config.graphics_config.camera_width = 1920;
                new_lbm_config.graphics_config.camera_height = 1080;
                new_lbm_config.graphics_config.streamline_every = 8;
                new_lbm_config.graphics_config.vec_vis_mode = graphics::VecVisMode::U;
                new_lbm_config.graphics_config.streamline_mode = true;
                new_lbm_config.graphics_config.axes_mode = true;
                new_lbm_config.graphics_config.q_mode = false;
                new_lbm_config.graphics_config.flags_surface_mode = true;
                new_lbm_config.graphics_config.flags_mode = true;

                match file::read(state.file_name.as_str(), &mut new_lbm_config) {
                    Ok(new_lbm) => { 
                        lbm = new_lbm;
                        lbm.initialize();
                    },
                    Err(msg) => { println!("ERROR: Could not read file \"{}\"!\nReason: {}", state.file_name.as_str(), msg)},
                }
                state.load_file = false;
            }
        }
        if !state.active {
            println!("\nExiting Simulation Loop");
            break;
        }
        let recieve_result = ctrl_rx.try_recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        } else {
            // Prevent spam polling while allowing smooth control of graphics
            thread::sleep(Duration::from_millis(10))
        }
    }
}

fn pull_state(rx: &mpsc::Receiver<SimControlTx>, state: &mut SimControlTx) {
    let recieve_result = rx.try_iter().last();
    if let Some(recieve) = recieve_result {
        *state = recieve;
    }
}
