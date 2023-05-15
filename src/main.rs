extern crate image;

use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

use std::{mem::swap};

struct Vec2 {
    //TODO: Implement
    x: f64,
    y: f64,
}

pub struct SimSize {
    // Struct that saves simulation sizes.
    // This type is used in all methods instead of the original "n", enabling arbitrary aspect ratios
    x: usize,
    y: usize,
}

pub fn add_source(x: &mut Array<f64, 2>, s: &mut Array<f64, 2>, dt: f64) {
    // Adds sources, used in the density and velocity steps
    // Can also be called independently of those
    for (c, item) in x.iter_mut() {
        *item += dt * s[c];
    }
}

pub fn diffuse(
    s: &SimSize,
    b: i32,
    x: &mut Array<f64, 2>,
    x0: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    let a = dt * diff * s.x as f64 * s.y as f64;
    lin_solve(&s, b, x, x0, a, 1.0 + 4.0 * a)
}

pub fn project(
    s: &SimSize,
    force_x: &mut Array<f64, 2>,
    force_y: &mut Array<f64, 2>,
    force_x_prev: &mut Array<f64, 2>,
    force_y_prev: &mut Array<f64, 2>,
) {
    // Variable names in this context:
    // u = u in global context
    // v = v in global context
    // p = u_prev
    // div = v_prev

    //HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
    let f = (s.x as f64 * s.y as f64).sqrt();

    for i in 1..=s.x {
        for j in 1..=s.y {
            force_y_prev[[i, j]] =
                -0.5 * (force_x[[i + 1, j]] - force_x[[i - 1, j]] + force_y[[i, j + 1]] - force_y[[i, j - 1]]) / f ;
                //HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
            force_x_prev[[i, j]] = 0.0;
        }
    }
    set_bnd(&s, 0, force_y_prev);
    set_bnd(&s, 0, force_x_prev);

    lin_solve(&s, 0, force_x_prev, force_y_prev, 1.0, 4.0);

    for xi in 1..=s.x {
        for yi in 1..=s.y {
            force_x[[xi, yi]] -= 0.5 * f * (force_x_prev[[xi + 1, yi]] - force_x_prev[[xi - 1, yi]]);
            force_y[[xi, yi]] -= 0.5 * f * (force_x_prev[[xi, yi + 1]] - force_x_prev[[xi, yi - 1]]);
        }
    }
    set_bnd(&s, 1, force_x);
    set_bnd(&s, 2, force_y);
}

fn lin_solve(s: &SimSize, b: i32, x: &mut Array<f64, 2>, x0: &mut Array<f64, 2>, a: f64, c: f64) {
    for _k in 0..20 {
        for xi in 1..=s.x {
            for yi in 1..=s.y {
                x[[xi, yi]] = (x0[[xi, yi]]
                    + a * (x[[xi - 1, yi]] + x[[xi + 1, yi]] + x[[xi, yi - 1]] + x[[xi, yi + 1]]))
                    / c;
            }
        }
        set_bnd(&s, b, x)
    }
}

pub fn advect(
    s: &SimSize,
    b: i32,
    d: &mut Array<f64, 2>,
    d0: &mut Array<f64, 2>,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    dt: f64,
) {
    //HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f

    let dt0x = dt * s.x as f64;
    let dt0y = dt * s.x as f64;
    for i in 1..=s.x {
        for j in 1..=s.y {
            let mut x = i as f64 - dt0x * u[[i, j]];
            let mut y = j as f64 - dt0y * v[[i, j]];
            if x < 0.5 {
                x = 0.5;
            }
            if x > s.x as f64 + 0.5 {
                x = s.x as f64 + 0.5;
            }
            let i0 = x as usize;
            let i1 = i0 + 1;
            if y < 0.5 {
                y = 0.5;
            };
            if y > s.y as f64 + 0.5 {
                y = s.y as f64 + 0.5
            };
            let j0 = y as usize;
            let j1 = j0 + 1;
            let s1 = x - i0 as f64;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f64;
            let t0 = 1.0 - t1;
            d[[i, j]] = s0 * (t0 * d0[[i0, j0]] + t1 * d0[[i0, j1]])
                + s1 * (t0 * d0[[i1, j0]] + t1 * d0[[i1, j1]]);
        }
    }
    set_bnd(&s, b, d);
}

pub fn set_bnd(s: &SimSize, b: i32, x: &mut Array<f64, 2>) {
    for i in 1..=s.x {
        x[[i, 0]] = if b == 2 { -x[[i, 1]] } else { x[[i, 1]] };
        x[[i, s.y + 1]] = if b == 2 { -x[[i, s.y]] } else { x[[i, s.y]] };
    }
    for i in 1..=s.y {
        x[[0, i]] = if b == 1 { -x[[1, i]] } else { x[[1, i]] };
        x[[s.x + 1, i]] = if b == 1 { -x[[s.x, i]] } else { x[[s.x, i]] };
    }

    x[[0, 0]] = 0.5 * (x[[1, 0]] + x[[0, 1]]);
    x[[0, s.y + 1]] = 0.5 * (x[[1, s.y + 1]] + x[[0, s.y]]);
    x[[s.x + 1, 0]] = 0.5 * (x[[s.x, 0]] + x[[s.x + 1, 1]]);
    x[[s.x + 1, s.y + 1]] = 0.5 * (x[[s.x, s.y + 1]] + x[[s.x + 1, s.y]]);
}

pub fn dens_step(
    s: &SimSize,
    dens: &mut Array<f64, 2>,
    dens_prev: &mut Array<f64, 2>,
    force_x: &mut Array<f64, 2>,
    force_y: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    add_source(dens, dens_prev, dt);
    swap(dens_prev, dens);
    diffuse(&s, 0, dens, dens_prev, diff, dt);
    swap(dens_prev, dens);
    advect(&s, 0, dens, dens_prev, force_x, force_y, dt);
}

pub fn vel_step(
    s: &SimSize,
    force_x: &mut Array<f64, 2>,
    force_y: &mut Array<f64, 2>,
    force_x_prev: &mut Array<f64, 2>,
    force_y_prev: &mut Array<f64, 2>,
    visc: f64,
    dt: f64,
) {
    add_source(force_x, force_x_prev, dt);
    add_source(force_y, force_y_prev, dt);
    swap(force_x_prev, force_x);
    diffuse(&s, 1, force_x, force_x_prev, visc, dt);
    swap(force_y_prev, force_y);
    diffuse(&s, 2, force_y, force_y_prev, visc, dt);
    project(&s, force_x, force_y, force_x_prev, force_y_prev);
    swap(force_x_prev, force_x);
    swap(force_y_prev, force_y);
    advect(&s, 1, force_x, force_x_prev, &mut force_x_prev.clone(), force_y_prev, dt);
    advect(&s, 2, force_y, force_y_prev, force_x_prev, &mut force_y_prev.clone(), dt);
    project(&s, force_x, force_y, force_x_prev, force_y_prev);
}

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
    // n is the size of the simulation. Size: n*n
    // visc: viscosity, default 0.0
    // diff: diffusion rate, default 0.0
    // dt: multiplier in add_source function, default 0.1
    let s = SimSize{x: 64, y: 128};
    // Arrays need to be 1px wider on each side, therefor n is used 2 bigger
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
