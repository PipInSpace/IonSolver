extern crate image;

use image::{ImageBuffer, Rgb};

use micro_ndarray::Array;

type Vec2 = (f64, f64);

use std::{io::stdin, mem::swap};

macro_rules! swap {
    ($o:expr, $a:tt, $b:tt) => {
        swap(&mut $o.$a, &mut $o.$b);
    };
}

struct Engine {
    n: usize,
    velo_x: Array<f64, 2>,
    velo_x_prev: Array<f64, 2>,
    velo_y: Array<f64, 2>,
    velo_y_prev: Array<f64, 2>,
    dens: Array<f64, 2>,
    dens_prev: Array<f64, 2>,
}

fn add_source(cur: &mut Array<f64, 2>, prev: &mut Array<f64, 2>, dt: f64) {
    for ([x, y], item) in cur.iter_mut() {
        *item += dt * prev[[x, y]];
    }
}
fn diffuse(n: usize, cur: &mut Array<f64, 2>, prev: &mut Array<f64, 2>, diff: f64, dt: f64) {
    let a = dt * diff * n as f64;
    for _k in 0..20 {
        let d = cur.clone();
        for (([x, y], item), (_, item_old)) in cur.iter_mut().zip(prev.iter()) {
            *item = *item_old
                + a * (d[[x - 1, y]] + d[[x + 1, y]] + d[[x, y - 1]] + d[[x, y + 1]])
                    / (1.0 + 4.0 * a);
        }
    }
}

impl Engine {
    /// dt = &delta;time
    /// b = ???
    #[cfg_attr(rustfmt, rustfmt_skip)]
    pub fn advect(&mut self, dt: f64) {
        let dt0 = dt * self.n as f64;
        for ([cx, cy], item) in self.dens.iter_mut() {
            let x = cx as f64 - dt0 * self.velo_x[[cx, cy]];
            let y = cy as f64 - dt0 * self.velo_y[[cx, cy]];
            let x = x.clamp(0.5, self.n as f64 + 0.5);
            let y = y.clamp(0.5, self.n as f64 + 0.5);

            // Interpolate: 
            //    xl  xh
            // yh X   X
            //      
            // yl X   X
            let factor_xh = x % 1.0;
            let factor_xl = 1.0 - factor_xh;
            let factor_yh = y % 1.0;
            let factor_yl = 1.0 - factor_yh;
            let x = x as usize;
            let y = y as usize;
            *item = 
            factor_xl * (
                factor_yl * self.dens_prev[[x, y]] + 
                factor_yh * self.dens_prev[[x, y + 1]]
            ) + 
            factor_xh * (
                factor_yl * self.dens_prev[[x + 1, y]] +
                factor_yh * self.dens_prev[[x + 1, y + 1]]
            );
        }
    }

    pub fn step_density(&mut self, diff: f64, dt: f64) {
        add_source(&mut self.dens, &mut self.dens_prev, dt);
        swap!(self, dens, dens_prev);
        diffuse(self.n, &mut self.dens, &mut self.dens_prev, diff, dt);
        swap!(self, dens, dens_prev);
        self.advect(dt);
    }

    pub fn step_velocity(&mut self, visc: f64, dt: f64) {
        add_source(&mut self.velo_x, &mut self.velo_x_prev, dt);
        add_source(&mut self.velo_y, &mut self.velo_y_prev, dt);
        swap!(self, velo_x, velo_x_prev);
        swap!(self, velo_y, velo_y_prev);
        diffuse(self.n, &mut self.velo_x, &mut self.velo_x_prev, visc, dt);
        diffuse(self.n, &mut self.velo_y, &mut self.velo_y_prev, visc, dt);
        // self.project();
    }
    

    pub fn step(&mut self) {}
}

pub fn project(N: usize, u: &mut Array<f32, 2>, v: &mut Array<f32, 2>, p: &mut Array<f32, 2>, div: &mut Array<f32, 2>) {
    // variable names in this context:
    // u = u in global context
    // v = v in global context
    // p = u_prev
    // div = v_prev

    let h = 1.0 / N as f32;
    for ([x, y], item) in div.iter_mut() {
        *item = -0.5
            * h
            * (u[[x + 1, y]] - u[[x - 1, y]] + v[[x, y + 1]]
                - v[[x, y - 1]]);
        p[[x, y]] = 0.0;
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    lin_solve(N, 0, p, div, 1.0, 4.0);

    for xi in 1..=N {
        for yi in 1..=N {
            u[[xi, yi]] -= 0.5 * N as f32 * (p[[xi + 1, yi]] - p[[xi - 1, yi]]);
            v[[xi, yi]] -= 0.5 * N as f32 * (p[[xi, yi + 1]] - p[[xi, yi - 1]]);
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);

}

fn lin_solve(N: usize, b: i32, x: &mut Array<f32, 2>, x0: &mut Array<f32, 2>, a: f32, c: f32) {
    for _k in 0..20 {
        for xi in 1..=N {
            for yi in 1..=N {
                x[[xi, yi]] = (x0[[xi, yi]] + a * (x[[xi-1, yi]] + x[[xi+1, yi]] + x[[xi, yi-1]] + x[[xi, yi+1]])) / c;
            }
        }
        set_bnd(N, b, x)
    }
}

pub fn set_bnd(N: usize, b: i32, x: &mut Array<f32, 2>) {
    for i in 1..=N {
        x[[0, i]] = if b == 1 {-x[[1,i]]} else {x[[1,i]]};
        x[[N+1, i]] = if b == 1 {-x[[N,i]]} else {x[[N,i]]};
        x[[i, 0]] = if b == 2 {-x[[i,1]]} else {x[[i,1]]};
        x[[i, N+1]] = if b == 2 {-x[[i,N]]} else {x[[i,N]]};
    }
    x[[0, 0]] = 0.5 * (x[[1, 0]] + x[[0, 1]]);
    x[[0, N+1]] = 0.5 * (x[[1, N+1]] + x[[0, N]]);
    x[[N+1, 0]] = 0.5 * (x[[N, 0]] + x[[N+1, 1]]);
    x[[N+1, N+1]] = 0.5 * (x[[N, N+1]] + x[[N+1, N]]);
}

fn main() {
    let mut img = ImageBuffer::new(640, 480);
    
    // Get space-seperated RGB input
    let num: Vec<u8> = stdin()
        .lines()
        .next()
        .unwrap()
        .unwrap()
        .split(" ")
        .map(str::parse)
        .map(Result::unwrap)
        .collect();

    // Set every pixel to input color
    for (_, _, pixel) in img.enumerate_pixels_mut() {
        *pixel = Rgb([num[0], num[1], num[2]]);
    }

    // Save the image as a PNG file
    img.save("output.png").unwrap();
}
