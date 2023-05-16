extern crate image;

use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

mod solver;
use solver::*;

// Visualisation functions and debug info

pub fn draw_spectrum(s: &SimSize, array: &Array<f64, 2>, step: i32, name: &'static str) {
    // exports png image of a 2D float array with range 0-1 in blue-green-red spectrum
    let mut img = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = ((4.0 * array[[x as usize + 1, y as usize + 1]] - 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        let g = ((-(4.0 * array[[x as usize + 1, y as usize + 1]] - 2.0).abs() + 2.0).clamp(0.0, 1.0) * 255.0) as u8; 
        let b = ((4.0 * -array[[x as usize + 1, y as usize + 1]] + 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    img.save(format!(r"out/{name}{step}.png")).unwrap();
}

pub fn draw_spectrum_relative(s: &SimSize, array: &Array<f64, 2>, step: i32, name: &'static str) {
    // exports png image of a 2D float array with dynamic range in blue-green-red spectrum.
    // Highest value is red, 0 is blue.
    let f = 1.0 / array
        .iter()
        .map(|x| *x.1)
        .max_by(|a, b| f64::partial_cmp(a, b).expect("boom"))
        .expect("empty array should not be possible");
    
    let mut img = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = ((4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        let g = ((-(4.0 * (array[[x as usize + 1, y as usize + 1]] * f) - 2.0).abs() + 2.0).clamp(0.0, 1.0) * 255.0) as u8; 
        let b = ((4.0 * -(array[[x as usize + 1, y as usize + 1]] * f) + 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    img.save(format!(r"out/{name}{step}.png")).unwrap();
}

pub fn draw_multichannel(
    s: &SimSize, 
    r_channel: &Array<f64, 2>, 
    g_channel: &Array<f64, 2>, 
    b_channel: &Array<f64, 2>, 
    step: i32, 
    name: &'static str) 
    {
    // exports png image of 3 2D float arrays with range 0-1 in the three RGB channels
    let mut img = ImageBuffer::new(s.x as u32, s.y as u32);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = (r_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        let g = (g_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        let b = (b_channel[[x as usize + 1, y as usize + 1]].clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    img.save(format!(r"out/{name}{step}.png")).unwrap();
}

pub fn print_maxval (x: &Array<f64, 2>, name: &'static str) {
    println!(
        "Max {name}: {:?}", x.iter().map(|x| *x.1).max_by(|a, b| f64::partial_cmp(a, b).expect("boom"))
    );
}

// main method

fn main() {
    // s is the size of the simulation. Size: s.x * s.y
    // visc: viscosity, default 0.0
    // diff: diffusion rate, default 0.0
    // dt: multiplier in add_source function, default 0.1
    let s = SimSize{x: 114, y: 64};
    // Arrays need to be 1px wider on each side, therefor s is used + 2
    let mut force_x: Array<f64, 2> = Array::new([s.x+2, s.y+2]);
    let mut force_y: Array<f64, 2> = Array::new([s.x+2, s.y+2]);
    let mut force_x_prev: Array<f64, 2> = Array::new([s.x+2, s.y+2]);
    let mut force_y_prev: Array<f64, 2> = Array::new([s.x+2, s.y+2]);

    let mut dens: Array<f64, 2> = Array::new_with([s.x+2, s.y+2], 0.0);
    let mut dens_prev: Array<f64, 2> = Array::new_with([s.x+2, s.y+2], 0.0);

    let visc = 0.0;
    let diff = 0.0;
    let dt = 0.1;

    for i in 0..400 {
        println!("Step {}", i);
        if i % 1 == 0 {
            //draw_spectrum(n, &x, i, "dens");
            draw_spectrum_relative(&s, &mut dens, i, "densRel")
            //draw_multichannel(n, &x, &x, &x, i, "densGrey");
            //draw_multichannel(n, &x, &u, &v, i, "combined");
        }

        if i < 200 {
            dens_prev[[10, 20]] += 1.0;
            force_x_prev[[10, 20]] = 5.0;
            //v[[20, 50]] += 20.0;
            dens_prev[[50, 50]] += 1.0;
            force_y_prev[[50, 50]] = -5.0;
        }

        vel_step(&s, &mut force_x, &mut force_y, &mut force_x_prev, &mut force_y_prev, visc, dt);
        dens_step(&s, &mut dens, &mut dens_prev, &mut force_x, &mut force_y, diff, dt);
    }
}
