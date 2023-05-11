extern crate image;

use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

use std::{io::stdin, mem::swap};

pub fn add_source(N: usize, x: &mut Array<f64, 2>, s: &mut Array<f64, 2>, dt: f64) {
    for (c, item) in x.iter_mut() {
        //*item += dt * s[c];
    }
}

pub fn diffuse(
    N: usize,
    b: i32,
    x: &mut Array<f64, 2>,
    x0: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    let a = dt * diff * N as f64 * N as f64;
    lin_solve(N, b, x, x0, a, 1.0 + 4.0 * a)
}

pub fn project(
    N: usize,
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

    let h = 1.0 / N as f64;
    for i in 1..=N {
        for j in 1..=N {
            div[[i, j]] =
                -0.5 * (u[[i + 1, j]] - u[[i - 1, j]] + v[[i, j + 1]] - v[[i, j - 1]]) / N as f64;
            p[[i, j]] = 0.0;
        }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    lin_solve(N, 0, p, div, 1.0, 4.0);

    for xi in 1..=N {
        for yi in 1..=N {
            u[[xi, yi]] -= 0.5 * N as f64 * (p[[xi + 1, yi]] - p[[xi - 1, yi]]);
            v[[xi, yi]] -= 0.5 * N as f64 * (p[[xi, yi + 1]] - p[[xi, yi - 1]]);
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
}

fn lin_solve(N: usize, b: i32, x: &mut Array<f64, 2>, x0: &mut Array<f64, 2>, a: f64, c: f64) {
    for _k in 0..20 {
        for xi in 1..=N {
            for yi in 1..=N {
                x[[xi, yi]] = (x0[[xi, yi]]
                    + a * (x[[xi - 1, yi]] + x[[xi + 1, yi]] + x[[xi, yi - 1]] + x[[xi, yi + 1]]))
                    / c;
            }
        }
        set_bnd(N, b, x)
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

pub fn set_bnd(N: usize, b: i32, x: &mut Array<f64, 2>) {
    for i in 1..=N {
        x[[0, i]] = if b == 1 { -x[[1, i]] } else { x[[1, i]] };
        x[[N + 1, i]] = if b == 1 { -x[[N, i]] } else { x[[N, i]] };
        x[[i, 0]] = if b == 2 { -x[[i, 1]] } else { x[[i, 1]] };
        x[[i, N + 1]] = if b == 2 { -x[[i, N]] } else { x[[i, N]] };
    }
    x[[0, 0]] = 0.5 * (x[[1, 0]] + x[[0, 1]]);
    x[[0, N + 1]] = 0.5 * (x[[1, N + 1]] + x[[0, N]]);
    x[[N + 1, 0]] = 0.5 * (x[[N, 0]] + x[[N + 1, 1]]);
    x[[N + 1, N + 1]] = 0.5 * (x[[N, N + 1]] + x[[N + 1, N]]);
}

pub fn dens_step(
    N: usize,
    x: &mut Array<f64, 2>,
    x0: &mut Array<f64, 2>,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    diff: f64,
    dt: f64,
) {
    add_source(N, x, x0, dt);
    swap(x0, x);
    diffuse(N, 0, x, x0, diff, dt);
    swap(x0, x);
    advect(N, 0, x, x0, u, v, dt);
}

pub fn vel_step(
    N: usize,
    u: &mut Array<f64, 2>,
    v: &mut Array<f64, 2>,
    u0: &mut Array<f64, 2>,
    v0: &mut Array<f64, 2>,
    visc: f64,
    dt: f64,
) {
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    swap(u0, u);
    diffuse(N, 1, u, u0, visc, dt);
    swap(v0, v);
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    swap(u0, u);
    swap(v0, v);
    advect(N, 1, u, u0, &mut u0.clone(), v0, dt);
    advect(N, 2, v, v0, u0, &mut v0.clone(), dt);
    project(N, u, v, u0, v0);
}

pub fn draw_dens(N: usize, dens: &Array<f64, 2>, step: i32, name: &'static str) {
    let mut img = ImageBuffer::new(N as u32 + 2, N as u32 + 2);
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let r = (dens[[x as usize, y as usize]] * 255.0) as u8;
        *pixel = Rgb([r, 255 - r, 255 - r]);
    }
    img.save(format!(r"out\ {name}{step}.png")).unwrap();
}

fn main() {
    //let mut img = ImageBuffer::new(640, 480);
    //
    //// Set every pixel to input color
    //for (_, _, pixel) in img.enumerate_pixels_mut() {
    //    *pixel = Rgb([num[0], num[1], num[2]]);
    //}

    //// Save the image as a PNG file
    //img.save("output.png").unwrap();
    let N: usize = 42;
    let mut u: Array<f64, 2> = Array::new([N, N]);
    let mut u0: Array<f64, 2> = Array::new([N, N]);
    let mut v: Array<f64, 2> = Array::new([N, N]);
    let mut v0: Array<f64, 2> = Array::new([N, N]);

    let mut x: Array<f64, 2> = Array::new_with([N, N], 0.5);
    let mut x0: Array<f64, 2> = Array::new([N, N]);

    let visc = 0.1;
    let diff = 0.05;
    let dt = 0.01;

    //x[[10, 10]] = 1.0;
    //x[[11, 10]] = 1.0;
    //x[[12, 10]] = 1.0;
    //x[[13, 10]] = 1.0;
    //x[[14, 10]] = 2.0;
    //v[[14, 10]] = 1.0;
    //u[[14, 10]] = 1.0;
    //x[[15, 10]] = 2.0;
    //v[[15, 10]] = 1.0;
    //u[[15, 10]] = 1.0;

    let N: usize = 40;
    for i in 0..10 {
        //println!(
        //    "MAX DENS: {:?}",
        //    x.iter()
        //        .map(|x| *x.1)
        //        .max_by(|a, b| f64::partial_cmp(a, b).expect("boom"))
        //);
        println!("Step {}", i);
        draw_dens(N, &x, i, "dens");

        x[[10, 20]] = 1.0;
        v[[10, 20]] = 0.0;
        u[[10, 20]] = 20.0;

        vel_step(N, &mut u, &mut v, &mut u0, &mut v0, visc, dt);
        dens_step(N, &mut x, &mut x0, &mut u, &mut v, diff, dt);
        //draw_dens(N, &u0, i, "velx");
        //draw_dens(N, &v0, i, "vely");
    }
}
