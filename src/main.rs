//#![allow(non_snake_case)]
extern crate ocl;
use std::time::Duration;
use std::{sync::mpsc, thread};

mod graphics;
mod info;
mod lbm;
mod opencl;
mod solver;
use egui::ColorImage;
use egui::TextBuffer;
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

    img: Vec<u8>,
}


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
            img: vec![],
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
    clear_images: bool,
    frame_spacing: u32,
    active: bool,
}

#[allow(unused)]

struct SimControl {
    //SimControl handles control/displaying of the Simulation over the UI
    paused: bool,
    save: bool,
    clear_images: bool,
    frame_spacing: u32,
    frame_spacing_str: String,
    ctrl_tx: mpsc::Sender<SimControlTx>,
    sim_rx: mpsc::Receiver<SimState>,
}

impl SimControl {
    pub fn reset_sim(&mut self) {
        self.paused = true;
    }

    
    pub fn send_control(&self){
        self.ctrl_tx.send(SimControlTx {
            paused: self.paused,
            save: self.save,
            clear_images: self.clear_images,
            frame_spacing: self.frame_spacing,
            active: true,
        })
        .expect("GUI cannot communicate with sim thread");
    }
}

impl App for SimControl {
    /// UI update loop
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let img: Vec<u8> = vec![1920*1080; 0];
        let recieve_result = self.sim_rx.try_recv();
        if let Ok(recieve) = recieve_result {
            img = recieve.img;
        }

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
                    self.send_control();
                    let _z = self.paused;
                }
                if ui
                    .add(egui::Button::new("Reset").rounding(0.0f32))
                    .clicked()
                {
                    self.reset_sim();
                }
                if ui
                    .add(egui::Button::new(if self.save {"Saving Enabled"} else {"Saving Disabled"})
                    .rounding(0.0f32))
                    .clicked()
                {
                    self.save = !self.save;
                    self.send_control();
                }
                if ui
                    .add(egui::Button::new(if self.clear_images {"Old Output Deleted"} else {"Old Output Kept"})
                    .rounding(0.0f32))
                    .clicked()
                {
                    self.clear_images = !self.clear_images;
                    self.send_control();
                }
                let mut text = self.frame_spacing_str.clone();
                if ui.add(egui::TextEdit::singleline(&mut text))
                    .changed()
                {
                    self.frame_spacing_str = text.clone(); // clone jik text is cbv here
                    let result = str::parse::<u32>(&text);
                    if let Ok(value) = result {
                        self.frame_spacing = value;
                        self.send_control();
                    }
                }
            })
        });
        //TODO: Display the Simulation
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.image(
                ui.ctx()
                    .load_texture(
                        "sim",
                        ColorImage::from_rgb([1920, 1080], spectrum_img.into_raw().as_slice()),
                        Default::default(),
                    )
                    .id(),
                ui.available_size(),
            );
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
        save: true,
        clear_images: true,
        frame_spacing: 10,
        frame_spacing_str: "1".to_string(),
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
            clear_images: true,
            frame_spacing: 10,
            active: false,
        })
        .expect("Exit Channel cannot reach Simulation thread");

    handle.join().unwrap();

    return;
}
