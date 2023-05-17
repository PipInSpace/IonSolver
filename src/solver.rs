use micro_ndarray::Array;

use std::mem::{self, swap};

use crate::debug::*;
use crate::{print_maxval, vector2::*};

/// Struct that saves simulation size as x and y.
/// This type is used in all methods instead of the original "n", enabling arbitrary aspect ratios
pub struct SimSize {
    /// width of the simulation
    pub x: usize,
    /// height of the simulation
    pub y: usize,
}

/// Adds sources, used in the density step. aka does things
pub fn add_source(x: &mut Array<f64, 2>, s: &mut Array<f64, 2>, dt: f64) {
    for (c, item) in x.iter_mut() {
        *item += dt * s[c];
    }
}

/// Adds sources, used in the velocity step. aka does things. Vec2 variant.
pub fn add_source_vec2(x: &mut Array<Vec2, 2>, s: &mut Array<Vec2, 2>, dt: f64) {
    for (c, item) in x.iter_mut() {
        *item = item.add(&s[c].scale(dt));
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
    set_bnd_vec2(s, force_prev, true);

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
        for i in 1..=s.x {
            force_prev[[i, 0]].x = force_prev[[i, 1]].x;
            force_prev[[i, s.y + 1]].x = force_prev[[i, s.y]].x;
        }
        for i in 1..=s.y {
            force_prev[[0, i]].x = force_prev[[1, i]].x;
            force_prev[[s.x + 1, i]].x = force_prev[[s.x, i]].x;
        }

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
    new_force: &mut Array<Vec2, 2>,
    force: &mut Array<Vec2, 2>,
    dt: f64,
) {
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
    // WARNING: 0 is also possible, special case
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
    swap(dens_prev, dens);
    diffuse(&s, 0, dens, dens_prev, diff, dt);
    swap(dens_prev, dens);
    let mass = get_mass(dens);
    advect(&s, 0, dens, dens_prev, force, dt);
    fix_mass(dens, mass);
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
    advect_vec2(&s, force, force_prev, dt);
    project(&s, force, force_prev);
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
