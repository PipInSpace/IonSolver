extern crate image;

use egui::ColorImage;
use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

use eframe::*;

mod solver;
use solver::*;

use std::{cmp::Ordering, time::Duration};

mod vector2;
use vector2::*;

// Visualisation functions and debug info

pub fn draw_spectrum(
    s: &SimSize,
    array: &Array<f64, 2>,
    step: i32,
    name: &'static str,
    save: bool,
) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // exports png image of a 2D float array with range 0-1 in blue-green-red spectrum
    let mut img = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r =
            ((4.0 * array[[x as usize + 1, y as usize + 1]] - 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        let g = ((-(4.0 * array[[x as usize + 1, y as usize + 1]] - 2.0).abs() + 2.0)
            .clamp(0.0, 1.0)
            * 255.0) as u8;
        let b =
            ((4.0 * -array[[x as usize + 1, y as usize + 1]] + 2.0).clamp(0.0, 1.0) * 255.0) as u8;
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
    fix: &mut bool,
) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // exports png image of a 2D float array with dynamic range in blue-green-red spectrum.
    // Highest value is red, 0 is blue.
    let f = 1.0
        / array
            .iter()
            .map(|x| *x.1)
            .max_by(|a, b| match f64::partial_cmp(a, b) {
                Some(x) => x,
                None => {
                    *fix = true;
                    Ordering::Equal
                }
            })
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

pub fn print_maxval(x: &Array<f64, 2>, name: &'static str) {
    println!(
        "Max {name}: {:?}",
        x.iter()
            .map(|x| *x.1)
            .max_by(|a, b| f64::partial_cmp(a, b).expect("boom"))
    );
}

struct SimState {
    s: SimSize,
    force: Array<Vec2, 2>,
    force_prev: Array<Vec2, 2>,

    working_dens: Array<f64, 2>,
    working_dens_prev: Array<f64, 2>,
    dens: Array<f64, 2>,
    dens_prev: Array<f64, 2>,
    visc: f64,
    diff: f64,
    dt: f64,
    step: i32,
}

impl SimState {
    pub fn new(s: SimSize, visc: f64, diff: f64, dt: f64) -> SimState {
        Self {
            force: Array::new([s.x + 2, s.y + 2]),
            force_prev: Array::new([s.x + 2, s.y + 2]),
            working_dens: Array::new([s.x + 2, s.y + 2]),
            working_dens_prev: Array::new([s.x + 2, s.y + 2]),
            dens: Array::new([s.x + 2, s.y + 2]),
            dens_prev: Array::new([s.x + 2, s.y + 2]),

            visc,
            diff,
            dt,
            step: 0,
            s,
        }
    }
}

impl App for SimState {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        println!("Step {}", self.step);
        let mut fix_dens = false;
        //draw_spectrum(n, &x, i, "dens", false);
        let spectrum_relative_img = draw_spectrum_relative(
            &self.s,
            &mut self.dens,
            self.step,
            "densRel",
            false,
            &mut fix_dens,
        );
        //draw_multichannel(n, &x, &x, &x, i, "densGrey", false);
        //draw_multichannel(n, &x, &u, &v, i, "combined", false);
        // fix INFINITY values
        if fix_dens {
            println!("Fixing...");
            for (c, item) in self.dens_prev.iter_mut() {
                *item = self.working_dens_prev[c];
                *item /= 1.0E32;
            }
            for (c, item) in self.dens.iter_mut() {
                *item = self.working_dens[c];
                *item /= 1.0E32;
            }
            println!("All items have been divided by 1E10 and one step was skipped.");
        }
        self.working_dens = self.dens.clone();
        self.working_dens_prev = self.dens_prev.clone();

        if self.step < 200 {
            self.dens[[10, 20]] += 1.0;
            self.force_prev[[10, 20]] = Vec2 { x: 5.0, y: 0.0 };
            //v[[20, 50]] += 20.0;
            self.dens[[50, 50]] += 1.0;
            self.force_prev[[50, 50]] = Vec2 { x: 0.0, y: -5.0 };
        }
        print_maxval(&self.dens, "densBef");
        vel_step(
            &self.s,
            &mut self.force,
            &mut self.force_prev,
            self.visc,
            self.dt,
        );
        print_maxval(&self.dens, "densAft");
        dens_step(
            &self.s,
            &mut self.dens,
            &mut self.dens_prev,
            &mut self.force,
            self.diff,
            self.dt,
        );
        print_maxval(&self.dens, "densAft");

        let size = spectrum_relative_img.dimensions();
        let size = [size.0 as usize, size.1 as usize];

        //Update gui
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("IonSolver Simulation");
            ui.label(format!("Step {}", self.step));
            ui.image(
                ui.ctx()
                    .load_texture(
                        "sim",
                        ColorImage::from_rgb(size, spectrum_relative_img.into_raw().as_slice()),
                        Default::default(),
                    )
                    .id(),
                ui.available_size(),
            );
        });

        self.step += 1;

        ctx.request_repaint_after(Duration::from_millis(1000));
    }
}

/// main method
fn main() {
    // s is the size of the simulation. Size: s.x * s.y
    // visc: viscosity, default 0.0
    // diff: diffusion rate, default 0.0
    // dt: delta-time, controls time step, default 0.1
    let _ = Vec2 { x: 0.0, y: 0.0 }.normalize();
    let s = SimSize { x: 64, y: 64 };
    // Arrays need to be 1px wider on each side, therefor s is used + 2
    let visc = 0.0;
    let diff = 0.0;
    let dt = 0.1;

    let sim = SimState::new(s, visc, diff, dt);

    //UI code
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(320.0, 240.0)),
        ..Default::default()
    };
    eframe::run_native(
        "IonSolver",
        options,
        Box::new(|_cc| Box::<SimState>::new(sim)),
    )
    .expect("unable to open window");
}
