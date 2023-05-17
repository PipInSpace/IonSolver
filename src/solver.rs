use micro_ndarray::Array;

use std::mem::{self, swap};

use crate::debug::*;

/// Struct that saves simulation size as x and y.
/// This type is used in all methods instead of the original "n", enabling arbitrary aspect ratios
pub struct SimSize {
    /// width of the simulation
    pub x: usize,
    /// height of the simulation
    pub y: usize,
}

/// Adds sources, used in the density and velocity steps. aka does things
pub fn add_source(x: &mut Array<f64, 2>, s: &mut Array<f64, 2>, dt: f64) {
    // Adds sources, used in the density and velocity steps
    // Can also be called independently of those
    for (c, item) in x.iter_mut() {
        *item += dt * s[c];
    }
}

/// Primes the lin_solve() function for diffusing grid cells
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
    // HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
    // This seems to work for now, and produces the same results for the same resolutions, but further testing is needed
    let f = (s.x as f64 * s.y as f64).sqrt();

    for i in 1..=s.x {
        for j in 1..=s.y {
            force_y_prev[[i, j]] = -0.5
                * (force_x[[i + 1, j]] - force_x[[i - 1, j]] + force_y[[i, j + 1]]
                    - force_y[[i, j - 1]])
                / f;
            //HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
            force_x_prev[[i, j]] = 0.0;
        }
    }
    set_bnd(&s, 0, force_y_prev);
    set_bnd(&s, 0, force_x_prev);

    lin_solve(&s, 0, force_x_prev, force_y_prev, 1.0, 4.0);

    for xi in 1..=s.x {
        for yi in 1..=s.y {
            force_x[[xi, yi]] -=
                0.5 * f * (force_x_prev[[xi + 1, yi]] - force_x_prev[[xi - 1, yi]]);
            force_y[[xi, yi]] -=
                0.5 * f * (force_x_prev[[xi, yi + 1]] - force_x_prev[[xi, yi - 1]]);
        }
    }
    set_bnd(&s, 1, force_x);
    set_bnd(&s, 2, force_y);
}

/// Uses Gauss-Seidel relaxation to solve a system of linear equations.
fn lin_solve(s: &SimSize, b: i32, x: &mut Array<f64, 2>, x0: &mut Array<f64, 2>, a: f64, c: f64) {
    // Gauss-Seidel relaxation
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
    //HERE is a potential problem with the dynamic sizes: prev: dt * n, now: dt * s.x, dt * s.y

    let f = (s.x as f64 * s.y as f64).sqrt();
    let dt0 = dt * f;
    for i in 1..=s.x {
        for j in 1..=s.y {
            let mut x = i as f64 - dt0 * u[[i, j]];
            let mut y = j as f64 - dt0 * v[[i, j]];
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
    
    print_sum(dens, "dens1");
    //add_source(dens, dens_prev, dt);
    swap(dens_prev, dens);
    diffuse(&s, 0, dens, dens_prev, diff, dt);
    print_sum(dens, "dens4");
    swap(dens_prev, dens);
    print_sum(dens, "dens5");
    let mass = get_mass(dens);
    advect(&s, 0, dens, dens_prev, force_x, force_y, dt);
    print_sum(dens, "dens6");
    fix_mass(dens, mass);
    print_sum(dens, "dens7");
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
    advect(
        &s,
        1,
        force_x,
        &mut force_x_prev.clone(),
        force_x_prev,
        force_y_prev,
        dt,
    );
    advect(
        &s,
        2,
        force_y,
        &mut force_y_prev.clone(),
        force_x_prev,
        force_y_prev,
        dt,
    );
    project(&s, force_x, force_y, force_x_prev, force_y_prev);
}

// Mass conservation: The advection of the density field looses or adds mass uncontrollably.
// Correction is needed
#[inline]
pub fn get_mass(dens: &Array<f64, 2>) -> f64 {
    dens.iter().map(|x| *x.1).sum::<f64>()
}
pub fn fix_mass(dens: &mut Array<f64, 2>, prev_mass: f64) {
    let factor = prev_mass / get_mass(dens);
    for (_, item) in dens.iter_mut() {
        *item *= factor;
    }
}