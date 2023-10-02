//#![allow(non_snake_case)]
extern crate ocl;
use std::time::Duration;
use std::{sync::mpsc, thread};

mod info;
mod lbm;
mod opencl;
mod solver;
use solver::*;

//use image::{ImageBuffer, Rgb};
//use egui::ColorImage;
use eframe::*;

// Simulation
#[allow(unused)]
pub struct SimSize {
    /// Width of the simulation
    pub x: usize,
    /// Height of the simulation
    pub y: usize,
}

#[allow(unused)]
pub struct SimState {
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
pub struct SimControlTx {
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
    /// UI update loop
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // let spectrum_img = draw_spectrum(&self.s, &self.dens, self.step, "densRel", self.save);

        // Prepare images
        // let size = spectrum_img.dimensions();
        // let size = [size.0 as usize, size.1 as usize];

        //Update gui
        egui::TopBottomPanel::top("top_controls").show(ctx, |ui| {
            ui.heading("IonSolver Simulation");
            //ui.separator();
            ui.horizontal(|ui| {
                if ui
                    .add(
                        egui::Button::new(if self.paused { "Play" } else { "Pause" })
                            .rounding(0.0f32),
                    )
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
                if ui
                    .add(egui::Button::new("Reset").rounding(0.0f32))
                    .clicked()
                {
                    self.reset_sim();
                }
            })
        });
        //TODO: Display the Simulation
        egui::CentralPanel::default().show(ctx, |ui| {
            //ui.image(
            //    ui.ctx()
            //        .load_texture(
            //            "sim",
            //            ColorImage::from_rgb(size, spectrum_img.into_raw().as_slice()),
            //            Default::default(),
            //        )
            //        .id(),
            //    ui.available_size(),
            //);
        });

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

    println!("{}", info::LOGO_COLOR.to_string());

    // UI params
    let icon_bytes = include_bytes!("../icons/IonSolver.png");
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1000.0, 500.0)),
        icon_data: load_icon(&icon_bytes.to_vec()),
        follow_system_theme: false,
        default_theme: Theme::Light,
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
        paused: true,
        save: false,
        ctrl_tx,
        sim_rx,
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
