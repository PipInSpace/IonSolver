extern crate image;

use egui::ColorImage;
use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

use eframe::*;

mod solver;
use solver::*;

use std::time::Duration;

mod vector2;
use vector2::*;

use crate::debug::{print_maxval, print_sum};

mod debug;

// Visualisation functions and debug info

pub fn draw_spectrum(
    s: &SimSize,
    array: &Array<f64, 2>,
    step: i32,
    name: &'static str,
    save: bool,
) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // exports png image of a 2D float array with dynamic range in blue-green-red spectrum.
    // Highest value is red, 0 is blue.
    let f = 0.5;

    let mut img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = ((4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).clamp(0.0, 1.0)
            * 255.0) as u8;
        let g = ((-(4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).abs() + 2.0)
            .clamp(0.0, 1.0)
            * 255.0) as u8;
        let b = ((4.0 * -(array[[x as usize + 1, y as usize + 1]] * f) + 2.0).clamp(0.0, 1.0)
            * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    if save {
        img.save(format!(r"out/{name}{step}.png")).unwrap();
    }
    img
}

pub fn draw_spectrum_relative(
    s: &SimSize,
    array: &Array<f64, 2>,
    step: i32,
    name: &'static str,
    save: bool,
) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // exports png image of a 2D float array with dynamic range in blue-green-red spectrum.
    // Highest value is red, 0 is blue.
    let f = 1.0
        / array
            .iter()
            .map(|x| *x.1)
            .max_by(|a, b| f64::partial_cmp(a, b).unwrap())
            .expect("empty array should not be possible");

    let mut img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = ((4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).clamp(0.0, 1.0)
            * 255.0) as u8;
        let g = ((-(4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).abs() + 2.0)
            .clamp(0.0, 1.0)
            * 255.0) as u8;
        let b = ((4.0 * -(array[[x as usize + 1, y as usize + 1]] * f) + 2.0).clamp(0.0, 1.0)
            * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    if save {
        img.save(format!(r"out/{name}{step}.png")).unwrap();
    }
    img
}

pub fn draw_multichannel(
    s: &SimSize,
    r_channel: &Array<f64, 2>,
    g_channel: &Array<f64, 2>,
    b_channel: &Array<f64, 2>,
    step: i32,
    name: &'static str,
    save: bool,
) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // exports png image of 3 2D float arrays with range 0-1 in the three RGB channels
    let mut img = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = (r_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        let g = (g_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        let b = (b_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    if save {
        img.save(format!(r"out/{name}{step}.png")).unwrap();
    }
    img
}

struct SimState {
    s: SimSize,
    force: Array<Vec2, 2>,
    force_prev: Array<Vec2, 2>,

    dens: Array<f64, 2>,
    dens_prev: Array<f64, 2>,
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
            force: Array::new([s.x + 2, s.y + 2]),
            force_prev: Array::new([s.x + 2, s.y + 2]),
            dens: Array::new([s.x + 2, s.y + 2]),
            dens_prev: Array::new([s.x + 2, s.y + 2]),

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
        self.force = Array::new([self.s.x + 2, self.s.y + 2]);
        self.force_prev = Array::new([self.s.x + 2, self.s.y + 2]);
        self.dens = Array::new([self.s.x + 2, self.s.y + 2]);
        self.dens_prev = Array::new([self.s.x + 2, self.s.y + 2]);
        self.step = 0;
        self.paused = true;
    }
}

impl App for SimState {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        //draw_spectrum(n, &x, i, "dens", false);
        let spectrum_img = draw_spectrum(&self.s, &self.dens, self.step, "densRel", self.save);
        //draw_multichannel(n, &x, &x, &x, i, "densGrey", false);
        //draw_multichannel(n, &x, &u, &v, i, "combined", false);

        if !self.paused {
            println!("Step {}", self.step);
            if self.step < 400 {
                self.dens[[20, 20]] += 4.0;
                self.force_prev[[20, 20]].x += 5.0;
                self.dens[[20, 50]] += 4.0;
                self.force_prev[[20, 50]].y -= 5.0;
                self.dens[[50, 20]] += 4.0;
                self.force_prev[[50, 20]].y += 5.0;
                self.dens[[50, 50]] += 4.0;
                self.force_prev[[50, 50]].x -= 5.0;
            }
            vel_step(
                &self.s,
                &mut self.force,
                &mut self.force_prev,
                self.visc,
                self.dt,
            );
            dens_step(
                &self.s,
                &mut self.dens,
                &mut self.dens_prev,
                &mut self.force,
                self.diff,
                self.dt,
            );
            print_sum(&self.dens, "dens");
            self.step += 1;
        }

        // Prepare images
        let size = spectrum_img.dimensions();
        let size = [size.0 as usize, size.1 as usize];

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
                    self.reset_sim();
                }
                if ui.text_edit_singleline(&mut self.dt_text).changed() {
                    match self.dt_text.parse() {
                        Ok(dt) => self.dt = dt,
                        Err(_) => {
                            if self
                                .dt_text
                                .contains(|c| !((c >= '0' && c <= '9') || c == '.'))
                            {
                                self.dt_text = self.dt.to_string();
                            }
                        }
                    }
                }
                ui.checkbox(&mut self.save, "Save");
                ui.label(format!("Step {}", self.step));
            })
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.image(
                ui.ctx()
                    .load_texture(
                        "sim",
                        ColorImage::from_rgb(size, spectrum_img.into_raw().as_slice()),
                        Default::default(),
                    )
                    .id(),
                ui.available_size(),
            );
        });

        if !self.paused {
            ctx.request_repaint_after(Duration::from_millis(0))
        }
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

/// main method, starts ui/sim loop
fn main() {
    // s is the size of the simulation. Size: s.x * s.y
    // visc: viscosity, default 0.0
    // diff: diffusion rate, default 0.0
    // dt: delta-time, controls time step, default 0.1
    let _ = Vec2 { x: 0.0, y: 0.0 }.normalize();
    let s = SimSize { x: 228, y: 128 };
    let visc = 0.000004;
    let diff = 0.00001;
    let dt = 0.1;

    let sim = SimState::new(s, visc, diff, dt);

    //UI code

    let icon_bytes = include_bytes!("../icons/IonSolver.png");

    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(320.0, 240.0)),
        icon_data: load_icon(&icon_bytes.to_vec()),
        ..Default::default()
    };

    eframe::run_native(
        "IonSolver",
        options,
        Box::new(|_| Box::<SimState>::new(sim)),
    )
    .expect("unable to open window");
}
