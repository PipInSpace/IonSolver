
use std::thread;
use std::time::Duration;

use image::{ImageBuffer, Rgb};
use egui::ColorImage;
use eframe::*;


// Simulation
struct SimSize {
    /// width of the simulation
    pub x: usize,
    /// height of the simulation
    pub y: usize,
}

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
struct SimControl {
    //SimControl handles control/displaying of the Simulation over the UI
    paused: bool,
    save: bool,
}

impl SimControl {

}

impl App for SimControl {
    /// UI update loop, called repeatedly if sim is not paused
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // let spectrum_img = draw_spectrum(&self.s, &self.dens, self.step, "densRel", self.save);

        if !self.paused {
            //TODO: self.update_sim();
        }

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
                }
                if ui.button("Reset").clicked() {
                    //TODO: self.reset_sim();
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


    // Start sim loop
    let handle = thread::spawn(move ||{
        //TODO: Setup mpsc messaging for data transmission
        let sim = SimState::new(s, visc, diff, dt);
        simloop(sim);
    });

    let simcontrol = SimControl { paused: false, save: false };

    // Start window loop
    eframe::run_native(
        "IonSolver",
        options,
        Box::new(|_| Box::<SimControl>::new(simcontrol)),
    )
    .expect("unable to open window");

    handle.join().unwrap();    

}

fn simloop(SimState { s, visc, diff, dt, dt_text, step, paused, save }: SimState) {
    //TODO: Simulation here

}