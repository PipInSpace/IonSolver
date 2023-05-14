extern crate image;

use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

use std::{mem::swap};

pub fn add_source(x: &mut Array<f64, 2>, s: &mut Array<f64, 2>, dt: f64) {
    for (c, item) in x.iter_mut() {
        *item += dt * s[c];
    }
}

pub fn diffuse(
    n: usize,
    b: i32,
    x: &mut Array<f64, 2>,
    x0: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    let a = dt * diff * n as f64 * n as f64;
    lin_solve(n, b, x, x0, a, 1.0 + 4.0 * a)
}

pub fn project(
    n: usize,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    p: &mut Array<f64, 2>,
    div: &mut Array<f64, 2>,
) {
    // variable names in this context:
    // u = u in global context
    // v = v in global context
    // p = u_prev
    // div = v_prev
    for i in 1..=n {
        for j in 1..=n {
            div[[i, j]] =
                -0.5 * (u[[i + 1, j]] - u[[i - 1, j]] + v[[i, j + 1]] - v[[i, j - 1]]) / n as f64;
            p[[i, j]] = 0.0;
        }
    }
    set_bnd(n, 0, div);
    set_bnd(n, 0, p);

    lin_solve(n, 0, p, div, 1.0, 4.0);

    for xi in 1..=n {
        for yi in 1..=n {
            u[[xi, yi]] -= 0.5 * n as f64 * (p[[xi + 1, yi]] - p[[xi - 1, yi]]);
            v[[xi, yi]] -= 0.5 * n as f64 * (p[[xi, yi + 1]] - p[[xi, yi - 1]]);
        }
    }
    set_bnd(n, 1, u);
    set_bnd(n, 2, v);
}

fn lin_solve(n: usize, b: i32, x: &mut Array<f64, 2>, x0: &mut Array<f64, 2>, a: f64, c: f64) {
    for _k in 0..20 {
        for xi in 1..=n {
            for yi in 1..=n {
                x[[xi, yi]] = (x0[[xi, yi]]
                    + a * (x[[xi - 1, yi]] + x[[xi + 1, yi]] + x[[xi, yi - 1]] + x[[xi, yi + 1]]))
                    / c;
            }
        }
        set_bnd(n, b, x)
    }
}

pub fn advect(
    n: usize,
    b: i32,
    d: &mut Array<f64, 2>,
    d0: &mut Array<f64, 2>,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    dt: f64,
) {
    let dt0 = dt * n as f64;
    for i in 1..=n {
        for j in 1..=n {
            let mut x = i as f64 - dt0 * u[[i, j]];
            let mut y = j as f64 - dt0 * v[[i, j]];
            if x < 0.5 {
                x = 0.5;
            }
            if x > n as f64 + 0.5 {
                x = n as f64 + 0.5;
            }
            let i0 = x as usize;
            let i1 = i0 + 1;
            if y < 0.5 {
                y = 0.5;
            };
            if y > n as f64 + 0.5 {
                y = n as f64 + 0.5
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
    set_bnd(n, b, d);
}

pub fn set_bnd(n: usize, b: i32, x: &mut Array<f64, 2>) {
    for i in 1..=n {
        x[[0, i]] = if b == 1 { -x[[1, i]] } else { x[[1, i]] };
        x[[n + 1, i]] = if b == 1 { -x[[n, i]] } else { x[[n, i]] };
        x[[i, 0]] = if b == 2 { -x[[i, 1]] } else { x[[i, 1]] };
        x[[i, n + 1]] = if b == 2 { -x[[i, n]] } else { x[[i, n]] };
    }
    x[[0, 0]] = 0.5 * (x[[1, 0]] + x[[0, 1]]);
    x[[0, n + 1]] = 0.5 * (x[[1, n + 1]] + x[[0, n]]);
    x[[n + 1, 0]] = 0.5 * (x[[n, 0]] + x[[n + 1, 1]]);
    x[[n + 1, n + 1]] = 0.5 * (x[[n, n + 1]] + x[[n + 1, n]]);
}

pub fn dens_step(
    n: usize,
    x: &mut Array<f64, 2>,
    x0: &mut Array<f64, 2>,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    add_source(x, x0, dt);
    swap(x0, x);
    diffuse(n, 0, x, x0, diff, dt);
    swap(x0, x);
    advect(n, 0, x, x0, u, v, dt);
}

pub fn vel_step(
    n: usize,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    u0: &mut Array<f64, 2>,
    v0: &mut Array<f64, 2>,
    visc: f64,
    dt: f64,
) {
    add_source(u, u0, dt);
    add_source(v, v0, dt);
    swap(u0, u);
    diffuse(n, 1, u, u0, visc, dt);
    swap(v0, v);
    diffuse(n, 2, v, v0, visc, dt);
    project(n, u, v, u0, v0);
    swap(u0, u);
    swap(v0, v);
    advect(n, 1, u, u0, &mut u0.clone(), v0, dt);
    advect(n, 2, v, v0, u0, &mut v0.clone(), dt);
    project(n, u, v, u0, v0);
}

pub fn draw_spectrum(n: usize, dens: &Array<f64, 2>, step: i32, name: &'static str) {
    // exports png image of a 2D float array with range 0-1 in blue-green-red spectrum
    let mut img = ImageBuffer::new(n as u32 + 2, n as u32 + 2);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = ((4.0 * dens[[x as usize, y as usize]] - 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        let g = ((-(4.0 * dens[[x as usize, y as usize]] - 2.0).abs() + 2.0).clamp(0.0, 1.0) * 255.0) as u8; 
        let b = ((4.0 * -dens[[x as usize, y as usize]] + 2.0).clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    img.save(format!(r"out/{name}{step}.png")).unwrap();
}

pub fn draw_multichannel(
    n: usize, 
    r_channel: &Array<f64, 2>, 
    g_channel: &Array<f64, 2>, 
    b_channel: &Array<f64, 2>, 
    step: i32, 
    name: &'static str) 
    {
    // exports png image of 3 2D float arrays with range 0-1 in the three RGB channels
    let mut img = ImageBuffer::new(n as u32 + 2, n as u32 + 2);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = (r_channel[[x as usize, y as usize]].clamp(0.0, 1.0) * 255.0) as u8;
        let g = (g_channel[[x as usize, y as usize]].clamp(0.0, 1.0) * 255.0) as u8;
        let b = (b_channel[[x as usize, y as usize]].clamp(0.0, 1.0) * 255.0) as u8;
        *pixel = Rgb([r, g, b]);
    }
    img.save(format!(r"out/{name}{step}.png")).unwrap();
}

fn main() {
    // N is the size of the simulation. Size: N*N
    // u stores x force, v stores y force
    // "Symbol"0 variables store previous values
    // visc: viscosity, default 0.0
    // diff: diffusion rate, default 0.0
    // dt: multiplier in add_source function, default 0.1
    let n: usize = 66;
    let mut u: Array<f64, 2> = Array::new([n, n]);
    let mut u0: Array<f64, 2> = Array::new([n, n]);
    let mut v: Array<f64, 2> = Array::new([n, n]);
    let mut v0: Array<f64, 2> = Array::new([n, n]);

    let mut x: Array<f64, 2> = Array::new_with([n, n], 0.0);
    let mut x0: Array<f64, 2> = Array::new_with([n, n], 0.0);

    let visc = 0.0;
    let diff = 0.0;
    let dt = 0.1;

    let n: usize = 64;

    for i in 0..200 {
        //println!(
        //    "MAX DENS: {:?}",
        //    x.iter()
        //        .map(|x| *x.1)
        //        .max_by(|a, b| f64::partial_cmp(a, b).expect("boom"))
        //);
        println!("Step {}", i);
        if i % 10 == 0 {
            draw_spectrum(n, &x, i, "dens");
            draw_multichannel(n, &x, &x, &x, i, "densGrey");
            draw_multichannel(n, &x, &u, &v, i, "combined");
        }

        x0[[20, 50]] += 30.0;
        u0[[20, 50]] += 50.0;
        //v[[20, 50]] += 20.0;

        vel_step(n, &mut u, &mut v, &mut u0, &mut v0, visc, dt);
        dens_step(n, &mut x, &mut x0, &mut u, &mut v, diff, dt);
    }
}
