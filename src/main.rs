extern crate ocl;
use ocl::ProQue;
use std::time::Duration;
use std::{sync::mpsc, thread};

//use image::{ImageBuffer, Rgb};
//use egui::ColorImage;
use eframe::*;

// Simulation
#[allow(unused)]
struct SimSize {
    /// Width of the simulation
    pub x: usize,
    /// Height of the simulation
    pub y: usize,
}

#[allow(unused)]
struct SimState {
    s: SimSize,
    //force: Array<Vec2, 2>,
    //force_prev: Array<Vec2, 2>,

    //dens: Array<f64, 2>,
    //dens_prev: Array<f64, 2>,
    visc: f64,
    diff: f64,
    dt: f64,
    dt_text: String,
    step: i32,

    //Control:
    paused: bool,
    save: bool,
}

#[allow(unused)]
impl SimState {
    pub fn new(s: SimSize, visc: f64, diff: f64, dt: f64) -> SimState {
        Self {
            //force: Array::new([s.x + 2, s.y + 2]),
            //force_prev: Array::new([s.x + 2, s.y + 2]),
            //dens: Array::new([s.x + 2, s.y + 2]),
            //dens_prev: Array::new([s.x + 2, s.y + 2]),
            visc,
            diff,
            dt,
            dt_text: "0.1".to_owned(),
            step: 0,
            s,
            paused: true,
            save: false,
        }
    }

    pub fn reset_sim(&mut self) {
        //self.force = Array::new([self.s.x + 2, self.s.y + 2]);
        //self.force_prev = Array::new([self.s.x + 2, self.s.y + 2]);
        //self.dens = Array::new([self.s.x + 2, self.s.y + 2]);
        //self.dens_prev = Array::new([self.s.x + 2, self.s.y + 2]);
        self.step = 0;
        self.paused = true;
    }

    pub fn update_sim(&mut self) {
        println!("Step {}", self.step);
        //TODO: Update function
        // Adds sources for the first 400 steps
        //if self.step < 400 {
        //    self.dens[[20, 20]] += 4.0;
        //    self.force_prev[[20, 20]].x += 5.0;
        //    self.dens[[20, 50]] += 4.0;
        //    self.force_prev[[20, 50]].y -= 5.0;
        //    self.dens[[50, 20]] += 4.0;
        //    self.force_prev[[50, 20]].y += 5.0;
        //    self.dens[[50, 50]] += 4.0;
        //    self.force_prev[[50, 50]].x -= 5.0;
        //}
        //vel_step(
        //    &self.s,
        //    &mut self.force,
        //    &mut self.force_prev,
        //    self.visc,
        //    self.dt,
        //);
        //dens_step(
        //    &self.s,
        //    &mut self.dens,
        //    &mut self.dens_prev,
        //    &mut self.force,
        //    self.diff,
        //    self.dt,
        //);
        self.step += 1;
    }
}

// UI+Control
#[allow(unused)]
struct SimControlTx {
    //SimControlTx is used to send information to the simulation thread
    paused: bool,
    save: bool,
    active: bool,
}

#[allow(unused)]

struct SimControl {
    //SimControl handles control/displaying of the Simulation over the UI
    paused: bool,
    save: bool,
    ctrl_tx: mpsc::Sender<SimControlTx>,
    sim_rx: mpsc::Receiver<SimState>,
}

impl SimControl {
    pub fn reset_sim(&mut self) {
        self.paused = true;
    }
}

impl App for SimControl {
    /// UI update loop, called repeatedly if sim is not paused
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // let spectrum_img = draw_spectrum(&self.s, &self.dens, self.step, "densRel", self.save);

        // Prepare images
        // let size = spectrum_img.dimensions();
        // let size = [size.0 as usize, size.1 as usize];

        //Update gui
        egui::TopBottomPanel::top("top_controls").show(ctx, |ui| {
            ui.heading("IonSolver Simulation");
            ui.separator();
            ui.horizontal(|ui| {
                if ui
                    .button(if self.paused { "Play" } else { "Pause" })
                    .clicked()
                {
                    self.paused = !self.paused;
                    self.ctrl_tx
                        .send(SimControlTx {
                            paused: self.paused,
                            save: self.save,
                            active: true,
                        })
                        .expect("GUI cannot communicate with sim thread");
                    let _z = self.paused;
                }
                if ui.button("Reset").clicked() {
                    self.reset_sim();
                }
            })
        });
        //TODO: Display the Simulation
        // egui::CentralPanel::default().show(ctx, |ui| {
        //     ui.image(
        //         ui.ctx()
        //             .load_texture(
        //                 "sim",
        //                 ColorImage::from_rgb(size, spectrum_img.into_raw().as_slice()),
        //                 Default::default(),
        //             )
        //             .id(),
        //         ui.available_size(),
        //     );
        // });

        ctx.request_repaint_after(Duration::from_millis(100))
    }
}

pub fn load_icon(icon_bytes: &Vec<u8>) -> Option<eframe::IconData> {
    if let Ok(image) = image::load_from_memory(icon_bytes) {
        let image = image.to_rgba8();
        let (width, height) = image.dimensions();
        Some(eframe::IconData {
            width,
            height,
            rgba: image.as_raw().to_vec(),
        })
    } else {
        None
    }
}

fn main() {
    // Simulation params
    let s = SimSize { x: 150, y: 100 };
    let visc = 0.000004;
    let diff = 0.00001;
    let dt = 0.05;

    // UI params
    let icon_bytes = include_bytes!("../icons/IonSolver.png");
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(320.0, 240.0)),
        icon_data: load_icon(&icon_bytes.to_vec()),
        ..Default::default()
    };

    let (sim_tx, sim_rx) = mpsc::channel();
    let (ctrl_tx, ctrl_rx) = mpsc::channel();
    let exit_ctrl_tx = ctrl_tx.clone();
    // Start Simulation loop
    let handle = thread::spawn(move || {
        //TODO: Setup mpsc messaging for data transmission
        let sim = SimState::new(s, visc, diff, dt);
        simloop(sim, sim_tx, ctrl_rx);
    });

    // setup simcontrol struc to pass params and channels to GUI
    let simcontrol = SimControl {
        paused: false,
        save: false,
        ctrl_tx: ctrl_tx,
        sim_rx: sim_rx,
    };

    // Start window loop
    eframe::run_native(
        "IonSolver",
        options,
        Box::new(|_| Box::<SimControl>::new(simcontrol)),
    )
    .expect("unable to open window");

    exit_ctrl_tx
        .send(SimControlTx {
            paused: true,
            save: false,
            active: false,
        })
        .expect("Exit Channel cannot reach Simulation thread");

    handle.join().unwrap();

    return;
}

#[allow(unused)]
fn simloop(sim: SimState, sim_tx: mpsc::Sender<SimState>, ctrl_rx: mpsc::Receiver<SimControlTx>) {
    //TODO: Simulation here
    //sim_tx.send(sim) sends data to the main window loop
    let mut state = SimControlTx {
        paused: false,
        save: false,
        active: true,
    };
    let mut i = 0;

    trivial().unwrap();

    let src = include_str!("kernels.cl");

    let pro_que = ProQue::builder().src(src).dims([1 << 6, 1 << 6]).build().unwrap();

    let buffer_a = pro_que.create_buffer::<f32>().unwrap();
    let buffer_b = pro_que.create_buffer::<f32>().unwrap();

    let setup_kernel = pro_que
        .kernel_builder("add")
        .arg(&buffer_a)
        .arg(1.0f32)
        .build().unwrap();

    let gauss_seidel_step_kernel = pro_que
        .kernel_builder("gauss_seidel_step")
        .arg(&buffer_a)
        .arg(&buffer_b)
        .arg(4.096f32)
        .arg(17.384f32)
        .build().unwrap();

    unsafe {
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
                    
                    println!("Step {}: The value at index [{}] is now '{}'!", i, 130, vec[130]);
                }

                //println!("Running {}", i);
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

fn trivial() -> ocl::Result<()> {
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
        for _i in 1..20 {
            gauss_seidel_step_kernel.enq()?;
        }
    }

    let mut vec = vec![0.0f32; buffer_a.len()];
    buffer_a.read(&mut vec).enq()?;

    println!("Gauss: The value at index [{}] is now '{}'!", 131, vec[131]);
    Ok(())
}
