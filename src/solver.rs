use micro_ndarray::Array;

use std::mem::swap;

use crate::{print_maxval, vector2::*};

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

/// Adds sources, used in the density and velocity steps. aka does things. Vec2 variant.
pub fn add_source_vec2(x: &mut Array<Vec2, 2>, s: &mut Array<Vec2, 2>, dt: f64) {
    // Adds sources, used in the density and velocity steps
    // Can also be called independently of those
    for (c, item) in x.iter_mut() {
        *item = item.add(&s[c].scale(dt));
        //x[c].x += dt * s[c].x;
        //x[c].y += dt * s[c].y;
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

/// Primes the lin_solve() function for diffusing grid cells
pub fn diffuse_vec2(
    s: &SimSize,
    x: &mut Array<Vec2, 2>,
    x0: &mut Array<Vec2, 2>,
    diff: f64,
    dt: f64,
) {
    let a = dt * diff * s.x as f64 * s.y as f64;
    lin_solve_vec2(&s, x, x0, a, 1.0 + 4.0 * a)
}

pub fn project(s: &SimSize, force: &mut Array<Vec2, 2>, force_prev: &mut Array<Vec2, 2>) {
    let f = (s.x as f64 * s.y as f64).sqrt();

    for i in 1..=s.x {
        for j in 1..=s.y {
            force_prev[[i, j]].y = -0.5
                * (force[[i + 1, j]].x - force[[i - 1, j]].x + force[[i, j + 1]].y
                    - force[[i, j - 1]].y)
                / f;
            force_prev[[i, j]].x = 0.0;
        }
    }
    // Combined
    //set_bnd(&s, 0, force_y_prev);
    //set_bnd(&s, 0, force_x_prev);
    //ONLY BOUNDRIES
    set_bnd_vec2(s, force_prev, true);

    //lin_solve(&s, 0, force_x_prev, force_y_prev, 1.0, 4.0);
    //Manual implementation of Gauss-Seidel relaxation. Code reuse would be too performance intensive.
    for _k in 0..20 {
        for xi in 1..=s.x {
            for yi in 1..=s.y {
                force_prev[[xi, yi]].x = (force_prev[[xi, yi]].y
                    + 1.0
                        * (force_prev[[xi - 1, yi]].x
                            + force_prev[[xi + 1, yi]].x
                            + force_prev[[xi, yi - 1]].x
                            + force_prev[[xi, yi + 1]].x))
                    / 4.0;
            }
        }
        //set_bnd(&s, 0, x)
        force_prev[[0, 0]].x = 0.5 * (force_prev[[1, 0]].x + force_prev[[0, 1]].x);
        force_prev[[0, s.y + 1]].x = 0.5 * (force_prev[[1, s.y + 1]].x + force_prev[[0, s.y]].x);
        force_prev[[s.x + 1, 0]].x = 0.5 * (force_prev[[s.x, 0]].x + force_prev[[s.x + 1, 1]].x);
        force_prev[[s.x + 1, s.y + 1]].x =
            0.5 * (force_prev[[s.x, s.y + 1]].x + force_prev[[s.x + 1, s.y]].x);
    }

    for xi in 1..=s.x {
        for yi in 1..=s.y {
            force[[xi, yi]].x -=
                0.5 * f * (force_prev[[xi + 1, yi]].x - force_prev[[xi - 1, yi]].x);
            force[[xi, yi]].y -=
                0.5 * f * (force_prev[[xi, yi + 1]].x - force_prev[[xi, yi - 1]].x);
        }
    }

    // Combined
    //set_bnd(&s, 1, force_x);
    //set_bnd(&s, 2, force_y);
    set_bnd_vec2(s, force, false);
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

fn lin_solve_vec2(s: &SimSize, x: &mut Array<Vec2, 2>, x0: &mut Array<Vec2, 2>, a: f64, c: f64) {
    // Gauss-Seidel relaxation
    for _k in 0..20 {
        for xi in 1..=s.x {
            for yi in 1..=s.y {
                x[[xi, yi]] = x0[[xi, yi]] // Center
                    .add(
                        &x[[xi - 1, yi]] // Left
                            .add(&x[[xi + 1, yi]]) // Right
                            .add(&x[[xi, yi - 1]]) // Up
                            .add(&x[[xi, yi + 1]]) // Down
                            .scale(a),
                    )
                    .scale(1.0 / c);
                //x[[xi, yi]].x = (x0[[xi, yi]].x + a * (x[[xi - 1, yi]].x + x[[xi + 1, yi]].x + x[[xi, yi - 1]].x + x[[xi, yi + 1]].x)) / c;
                //x[[xi, yi]].y = (x0[[xi, yi]].y + a * (x[[xi - 1, yi]].y + x[[xi + 1, yi]].y + x[[xi, yi - 1]].y + x[[xi, yi + 1]].y)) / c;
            }
        }
        set_bnd_vec2(&s, x, false)
    }
}

pub fn advect(
    s: &SimSize,
    b: i32,
    d: &mut Array<f64, 2>,
    d0: &mut Array<f64, 2>,
    force: &mut Array<Vec2, 2>,
    dt: f64,
) {
    //HERE is a potential problem with the dynamic sizes: prev: dt * n, now: dt * s.x, dt * s.y

    let dt0 = dt * (s.x as f64 * s.y as f64).sqrt();
    for i in 1..=s.x {
        for j in 1..=s.y {
            let mut x = i as f64 - dt0 * force[[i, j]].x;
            let mut y = j as f64 - dt0 * force[[i, j]].y;
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

/// This function advects velocitys along themselves (along force)
pub fn advect_vec2(
    s: &SimSize,
    new_force: &mut Array<Vec2, 2>, //Is overwritten, not read at all
    force: &mut Array<Vec2, 2>,     //Is accessed for all steps
    dt: f64,
) {
    //HERE is a potential problem with the dynamic sizes: prev: dt * n, now: dt * s.x, dt * s.y

    let dt0 = dt * (s.x as f64 * s.y as f64).sqrt();
    for i in 1..=s.x {
        for j in 1..=s.y {
            let mut x = i as f64 - dt0 * force[[i, j]].x;
            let mut y = j as f64 - dt0 * force[[i, j]].y;
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
            //new_force[[i, j]].x = s0 * (t0 * force[[i0, j0]].x + t1 * force[[i0, j1]].x)
            //    + s1 * (t0 * force[[i1, j0]].x + t1 * force[[i1, j1]].x);
            //new_force[[i, j]].y = s0 * (t0 * force[[i0, j0]].y + t1 * force[[i0, j1]].y)
            //    + s1 * (t0 * force[[i1, j0]].y + t1 * force[[i1, j1]].y);

            new_force[[i, j]] = force[[i0, j0]]
                .scale(t0)
                .add(&force[[i0, j1]].scale(t1))
                .scale(s0)
                .add(
                    &force[[i1, j0]]
                        .scale(t0)
                        .add(&force[[i1, j1]].scale(t1))
                        .scale(s1),
                )
        }
    }
    set_bnd_vec2(&s, new_force, false);
}

pub fn set_bnd(s: &SimSize, b: i32, x: &mut Array<f64, 2>) {
    // b: 1 means x is being processed, 2 means y
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

pub fn set_bnd_vec2(s: &SimSize, x: &mut Array<Vec2, 2>, no_negate: bool) {
    if !no_negate {
        // Top and bottom row
        for i in 1..=s.x {
            x[[i, 0]] = x[[i, 1]].flip_y();
            x[[i, s.y + 1]] = x[[i, s.y]].flip_y();
        }
        // Left and right column
        for i in 1..=s.y {
            x[[0, i]] = x[[1, i]].flip_x();
            x[[s.x + 1, i]] = x[[s.x, i]].flip_x();
        }
    } else {
        // Top and bottom row
        for i in 1..=s.x {
            x[[i, 0]] = x[[i, 1]].scale(1.0);
            x[[i, s.y + 1]] = x[[i, s.y]].scale(1.0);
        }
        // Left and right column
        for i in 1..=s.y {
            x[[0, i]] = x[[1, i]].scale(1.0);
            x[[s.x + 1, i]] = x[[s.x, i]].scale(1.0);
        }
    }

    x[[0, 0]] = x[[1, 0]].add(&x[[0, 1]]).scale(0.5);
    x[[0, s.y + 1]] = x[[1, s.y + 1]].add(&x[[0, s.y]]).scale(0.5);
    x[[s.x + 1, 0]] = x[[s.x, 0]].add(&x[[s.x + 1, 1]]).scale(0.5);
    x[[s.x + 1, s.y + 1]] = x[[s.x, s.y + 1]].add(&x[[s.x + 1, s.y]]).scale(0.5);
}

/// The density step of the simulation
pub fn dens_step(
    s: &SimSize,
    dens: &mut Array<f64, 2>,
    dens_prev: &mut Array<f64, 2>,
    force: &mut Array<Vec2, 2>,
    diff: f64,
    dt: f64,
) {
    //add_source(dens, dens_prev, dt);
    print_maxval(&dens, "dens1");
    swap(dens_prev, dens);
    diffuse(&s, 0, dens, dens_prev, diff, dt);
    print_maxval(&dens, "dens1");
    swap(dens_prev, dens);
    advect(&s, 0, dens, dens_prev, force, dt);
    print_maxval(&dens, "dens1");
}

/// The velocity step of the simulation
pub fn vel_step(
    s: &SimSize,
    force: &mut Array<Vec2, 2>,
    force_prev: &mut Array<Vec2, 2>,
    visc: f64,
    dt: f64,
) {
    add_source_vec2(force, force_prev, dt);
    swap(force_prev, force);
    diffuse_vec2(&s, force, force_prev, visc, dt);

    project(&s, force, force_prev);

    // Combined
    swap(force_prev, force);

    // Combined.
    //advect(
    //    &s,
    //    1,
    //    force_x, //This thing get outputted, is never accessed, contains garbage at first
    //    force_x_prev, //This temporarily stores the current force, is accessed
    //    &mut force_x_prev.clone(), // These store the current force and are
    //    force_y_prev,              // used for advection
    //    dt,
    //);
    //advect(
    //    &s,
    //    2,
    //    force_y,
    //    force_y_prev,
    //    force_x_prev,
    //    &mut force_y_prev.clone(),
    //    dt,
    //);

    advect_vec2(&s, force, force_prev, dt);
    project(&s, force, force_prev);
}
