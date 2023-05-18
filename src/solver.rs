use micro_ndarray::{
    vec_split::{
        accessors::{Accessor, AccessorMut, IterateAccessorMut},
        SizedVectorArray,
    },
    Array,
};

use std::mem::swap;

use crate::vector2::Vec2;

type Idx = [usize; 2];

/// Struct that saves simulation size as x and y.
/// This type is used in all methods instead of the original "n", enabling arbitrary aspect ratios
pub struct SimSize {
    /// width of the simulation
    pub x: usize,
    /// height of the simulation
    pub y: usize,
}

/// Adds sources, used in the density and velocity steps. aka does things
pub fn add_source<'a>(
    size: &SimSize,
    x: &mut dyn AccessorMut<f64, Idx>,
    s: &dyn Accessor<f64, Idx>,
    dt: f64,
) {
    for i in 0..size.x + 2 {
        for j in 0..size.y + 2 {
            x[[i, j]] += dt * s[[i, j]];
        }
    }
}

/// Primes the lin_solve() function for diffusing grid cells
pub fn diffuse(
    s: &SimSize,
    b: i32,
    x: &mut dyn AccessorMut<f64, Idx>,
    x0: &dyn Accessor<f64, Idx>,
    diff: f64,
    dt: f64,
) {
    let a = dt * diff * s.x as f64 * s.y as f64;
    lin_solve(&s, b, x, x0, a, 1.0 + 4.0 * a)
}

pub fn project(s: &SimSize, force: &mut Array<Vec2, 2>, force_prev: &mut Array<Vec2, 2>) {
    // HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
    // This seems to work for now, and produces the same results for the same resolutions, but further testing is needed
    let f = (s.x as f64 * s.y as f64).sqrt();

    for i in 1..=s.x {
        for j in 1..=s.y {
            force_prev[[i, j]].y = -0.5
                * (force[[i + 1, j]].x - force[[i - 1, j]].x + force[[i, j + 1]].y
                    - force[[i, j - 1]].y)
                / f;
            //HERE is a potential problem with the dynamic sizes: prev: / n, now: / (s.x * s.y).sqrt() represented by f
            force_prev[[i, j]].x = 0.0;
        }
    }
    {
        let [mut force_x_prev, mut force_y_prev] = force_prev.vec_split_fast_mut();

        set_bnd(&s, 0, &mut force_y_prev);
        set_bnd(&s, 0, &mut force_x_prev);

        lin_solve(&s, 0, &mut force_x_prev, &force_y_prev, 1.0, 4.0);
    }

    for xi in 1..=s.x {
        for yi in 1..=s.y {
            force[[xi, yi]].x -=
                0.5 * f * (force_prev[[xi + 1, yi]].x - force_prev[[xi - 1, yi]].x);
            force[[xi, yi]].y -=
                0.5 * f * (force_prev[[xi, yi + 1]].x - force_prev[[xi, yi - 1]].x);
        }
    }
    {
        let [mut force_x, mut force_y] = force.vec_split_fast_mut();
        set_bnd(&s, 1, &mut force_x);
        set_bnd(&s, 2, &mut force_y);
    }
}

/// Uses Gauss-Seidel relaxation to solve a system of linear equations.
fn lin_solve(
    s: &SimSize,
    b: i32,
    x: &mut dyn AccessorMut<f64, Idx>,
    x0: &dyn Accessor<f64, Idx>,
    a: f64,
    c: f64,
) {
    // Gauss-Seidel relaxation
    for _k in 0..10 {
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
    d: &mut dyn AccessorMut<f64, Idx>,
    d0: &dyn Accessor<f64, Idx>,
    force: &Array<Vec2, 2>,
    dt: f64,
) {
    //HERE is a potential problem with the dynamic sizes: prev: dt * n, now: dt * s.x, dt * s.y

    let f = (s.x as f64 * s.y as f64).sqrt();
    let dt0 = dt * f;
    for i in 1..=s.x {
        for j in 1..=s.y {
            let f = force[[i, j]];
            let mut x = i as f64 - dt0 * f.x;
            let mut y = j as f64 - dt0 * f.y;
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

pub fn set_bnd(s: &SimSize, b: i32, x: &mut dyn AccessorMut<f64, Idx>) {
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

pub fn vel_step(
    s: &SimSize,
    force: &mut Array<Vec2, 2>,
    force_prev: &mut Array<Vec2, 2>,
    visc: f64,
    dt: f64,
) {
    {
        let [mut force_x, mut force_y] = force.vec_split_fast_mut();
        let [mut force_x_prev, mut force_y_prev] = force_prev.vec_split_fast_mut();
        add_source(s, &mut force_x, &force_x_prev, dt);
        add_source(s, &mut force_y, &force_y_prev, dt);
        swap(&mut force_x_prev, &mut force_x);
        swap(&mut force_y_prev, &mut force_y);
        diffuse(&s, 1, &mut force_x, &mut force_x_prev, visc, dt);
        diffuse(&s, 2, &mut force_y, &mut force_y_prev, visc, dt);
    }
    // swapping the split vars doesnt swap the actual data, so we do that here.
    swap(force, force_prev);
    project(&s, force, force_prev);
    swap(force_prev, force);
    {
        let [mut force_x, mut force_y] = force.vec_split_fast_mut();
        let [force_x_prev, force_y_prev] = force_prev.vec_split_fast();
        advect(&s, 1, &mut force_x, &force_x_prev, force_prev, dt);
        advect(&s, 2, &mut force_y, &force_y_prev, force_prev, dt);
    }
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
