use std::ops::Mul;

use crate::{info, lbm::Lbm};
use rayon::prelude::*;

/// calculates electric field vector at a cell with index n
/// from a vector of charges.
fn calculate_e(
    n: u64,
    charges: &[([f32; 3], f32)],
    lengths: (u32, u32, u32),
    def_ke: f32,
) -> [f32; 3] {
    let convert = move |n| -> [f32; 3] {
        [
            (n % lengths.1 as u64) as f32,
            (n / lengths.1 as u64 % lengths.2 as u64) as f32,
            (n / lengths.1 as u64 / lengths.2 as u64) as f32,
        ]
    };

    let coord = convert(n);
    let mut e_at_cell = [0.0; 3];
    for &(coord_charge, charge) in charges.iter() {
        let coord_diff = [
            coord[0] - coord_charge[0],
            coord[1] - coord_charge[1],
            coord[2] - coord_charge[2],
        ];
        let normalized = fast_normalize(coord_diff);
        let length_sq = len_sq(coord_diff);
        let length_sq_inv = 1.0 / length_sq;
        if length_sq != 0.0 {
            e_at_cell = [
                e_at_cell[0] + (charge * length_sq_inv * normalized[0]),
                e_at_cell[1] + (charge * length_sq_inv * normalized[1]),
                e_at_cell[2] + (charge * length_sq_inv * normalized[2]),
            ]
        }
    }
    e_at_cell[0] *= def_ke;
    e_at_cell[1] *= def_ke;
    e_at_cell[2] *= def_ke;
    e_at_cell
}

/// Converts charge u64 index positions to float3 positions
fn charge_float_pos(charges: Vec<(u64, f32)>, lengths: (u32, u32, u32)) -> Vec<([f32; 3], f32)> {
    let mut charges_vector_pos: Vec<([f32; 3], f32)> = Vec::with_capacity(charges.len());
    let convert = move |n| -> [f32; 3] {
        [
            (n % lengths.1 as u64) as f32,
            (n / lengths.1 as u64 % lengths.2 as u64) as f32,
            (n / lengths.1 as u64 / lengths.2 as u64) as f32,
        ]
    };

    // Precompute position vectors
    for &(i, charge) in charges.iter() {
        let coord_charge = convert(i);
        charges_vector_pos.push((coord_charge, charge))
    }

    charges_vector_pos
}

#[allow(unused_mut)]
// variables are mutated with deborrow
/// precomputes the electric field from a Vector of charges
pub fn precompute_E(lbm: &Lbm, charges: Vec<(u64, f32)>) {
    // TODO: Make multi-domain compatible
    println!("Precomputing E for {} charges", charges.len());

    // Set variables
    let n = lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64;
    let mut e_field: Vec<f32> = vec![0.0; (n * 3) as usize];
    let lengths: (u32, u32, u32) = (lbm.config.n_x, lbm.config.n_y, lbm.config.n_z);
    let def_ke = lbm.config.units.si_to_ke(8.987552E9);

    fn deborrow<'b, T>(r: &T) -> &'b mut T {
        // Neded to access e_field in parallel.
        // This is safe, because no indecies are accessed multiple times
        unsafe {
            #[allow(mutable_transmutes)]
            std::mem::transmute(r)
        }
    }

    let charges_float_pos: Vec<([f32; 3], f32)> = charge_float_pos(charges, lengths);

    let mut count: u32 = 0;
    (0..n).into_par_iter().for_each(|i| {
        let e_at = calculate_e(i, &charges_float_pos, lengths, def_ke);
        deborrow(&e_field)[i as usize] = e_at[0];
        deborrow(&e_field)[(i + n) as usize] = e_at[1];
        deborrow(&e_field)[(i + (n * 2)) as usize] = e_at[2];
        *deborrow(&count) += 1;
        print!("\r{}", info::progressbar(count as f32 / n as f32));
    });
    println!();

    lbm.domains[0]
        .e
        .as_ref()
        .expect("E buffer used but not initialized")
        .write(&e_field)
        .enq()
        .unwrap();
}

#[inline]
fn len_sq(v: [f32; 3]) -> f32 {
    v[0].sq() + v[1].sq() + v[2].sq()
}

// Fast vector normalization
fn fast_normalize(v: [f32; 3]) -> [f32; 3] {
    let len = fast_inv_sqrt(len_sq(v));
    v.map(|x| x * len)
}

// Fast inverse square root algorithm
fn fast_inv_sqrt(x: f32) -> f32 {
    let i = x.to_bits();
    let i = 0x5f3759df - (i >> 1);
    let y = f32::from_bits(i);

    y * (1.5 - 0.5 * x * y * y)
}

trait Sq: Copy + Sized + Mul<Self, Output = Self> {
    #[inline]
    fn sq(self) -> Self {
        self * self
    }
}

impl<T> Sq for T where T: Copy + Sized + Mul<Self, Output = Self> {}
