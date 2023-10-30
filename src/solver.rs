extern crate ocl;
use crate::lbm::*;
use crate::*;
use std::f32::consts::PI;
use std::io::Write;
use std::sync::mpsc;
use std::{fs, io};

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
    lbm_config.n_z = 256;
    lbm_config.d_x = 1;
    lbm_config.nu = 0.1;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    let mut test_lbm = Lbm::new(lbm_config);
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
                    print!("\rStep {}", i);
                    io::stdout().flush().unwrap();
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
            println!("\nExiting Simulation Loop");
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

        let pif = std::f32::consts::PI;
        #[allow(non_snake_case)]
        let A = 0.25f32;
        let periodicity = 2u32;
        let a = nx as f32 / periodicity as f32;
        let b = ny as f32 / periodicity as f32;
        let c = nz as f32 / periodicity as f32;

        let domain_numbers: u32 = self.config.d_x * self.config.d_y * self.config.d_z;
        let dx = self.config.d_x;
        let dy = self.config.d_y;
        let dz = self.config.d_z;
        let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64 * 2; // Domain size on each axis
        let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64 * 2; // Needs to account for halo offsets
        let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64 * 2;
        let dtotal = dsx * dsy * dsz;

        print!("");
        for d in 0..domain_numbers {
            let x = (d % (dx * dy)) % dx; // Current Domain coordinates
            let y = (d % (dx * dy)) / dx;
            let z = d / (dx * dy);

            let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
            #[allow(unused_mut)]
            let mut domain_vec_q: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize]; // temporary q (all 0)
            let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
            for zi in 0..dsz {
                // iterates over every cell in the domain, filling it with  Taylor-Green-vortex
                print!(
                    "\r{}",
                    info::progressbar(
                        ((zi as f32 + 1.0) / dsz as f32) * (1.0 / domain_numbers as f32)
                            + ((d as f32) / domain_numbers as f32)
                    )
                );
                for yi in 0..dsy {
                    for xi in 0..dsx {
                        if !(((xi == 0 || xi == dsx - 1) && dx > 1)
                            || ((yi == 0 || yi == dsy - 1) && dy > 1)
                            || ((zi == 0 || zi == dsz - 1) && dz > 1))
                        {
                            //do not set at halo offsets
                            let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                            let gx =
                                xi - (dx > 1u32) as u64 + x as u64 * (dsx - (dx > 1u32) as u64 * 2); // Global coordinates
                            let gy =
                                yi - (dy > 1u32) as u64 + y as u64 * (dsy - (dy > 1u32) as u64 * 2);
                            let gz =
                                zi - (dz > 1u32) as u64 + z as u64 * (dsz - (dz > 1u32) as u64 * 2);

                            let fx = gx as f32 + 0.5 - 0.5 * nx as f32;
                            let fy = gy as f32 + 0.5 - 0.5 * ny as f32;
                            let fz = gz as f32 + 0.5 - 0.5 * nz as f32;

                            domain_vec_u[(dn) as usize] = A
                                * (2.0 * pif * fx / a).cos()
                                * (2.0 * pif * fy / b).sin()
                                * (2.0 * pif * fz / c).sin(); // x
                            domain_vec_u[(dn + dtotal) as usize] = -A
                                * (2.0 * pif * fx / a).sin()
                                * (2.0 * pif * fy / b).cos()
                                * (2.0 * pif * fz / c).sin(); // y;
                            domain_vec_u[(dn + dtotal * 2) as usize] = A
                                * (2.0 * pif * fx / a).sin()
                                * (2.0 * pif * fy / b).sin()
                                * (2.0 * pif * fz / c).cos(); // z
                            domain_vec_rho[(dn) as usize] = 1.0
                                - (A * A) * 3.0 / 4.0 * (4.0 * pif * fx / a).cos()
                                + (4.0 * pif * fy / b).cos();
                        }
                    }
                }
            }

            // Write to domain buffers
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
            if self.config.ext_electric_force {
                // Only write to buffer if it is needed/initialized
                self.domains[d as usize]
                    .q
                    .as_mut()
                    .expect("q buffer used but not initialized")
                    .write(&domain_vec_q)
                    .enq()
                    .unwrap();
            }
            self.domains[d as usize].queue.finish().unwrap();
        }
        println!();
        println!("Finished setting up Taylor-Green vorticies");
    }
}
