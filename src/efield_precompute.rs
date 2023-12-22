use std::ops::Mul;

use crate::*;
use rayon::prelude::*;

/// calculates electric field vector at a cell with index n
/// from a vector of charges.
fn calculate_e(
    n: u64,
    charges: &[([u32; 3], f32)],
    lengths: (u32, u32, u32),
    def_ke: f32,
) -> [f32; 3] {

    let coord_cell = [
        (n % lengths.1 as u64) as u32,
        (n / lengths.1 as u64 % lengths.2 as u64) as u32,
        (n / lengths.1 as u64 / lengths.2 as u64) as u32,
    ];
    let mut e_at_cell = [0.0; 3];

    for &(coord_charge, charge) in charges.iter() {
        let coord_diff = [
            coord_cell[0] - coord_charge[0],
            coord_cell[1] - coord_charge[1],
            coord_cell[2] - coord_charge[2],
        ];
        let length_sq = len_sq_u32(coord_diff);
        if length_sq != 0.0 {
            let charge_length_sq_inv = charge * (1.0 / length_sq);
            let normalized = fast_normalize_u32(coord_diff);
            e_at_cell = [
                e_at_cell[0] + (charge_length_sq_inv * normalized[0]),
                e_at_cell[1] + (charge_length_sq_inv * normalized[1]),
                e_at_cell[2] + (charge_length_sq_inv * normalized[2]),
            ]
        }
    }
    e_at_cell[0] *= def_ke;
    e_at_cell[1] *= def_ke;
    e_at_cell[2] *= def_ke;
    e_at_cell
}

/// Converts charge u64 index positions to u32 coords
fn charge_u32_pos(charges: Vec<(u64, f32)>, lengths: (u32, u32, u32)) -> Vec<([u32; 3], f32)> {
    let mut charges_vector_pos: Vec<([u32; 3], f32)> = Vec::with_capacity(charges.len());

    // Precompute position vectors
    for &(i, charge) in charges.iter() {
        let coord_charge = [
            (i % lengths.1 as u64) as u32,
            (i / lengths.1 as u64 % lengths.2 as u64) as u32,
            (i / lengths.1 as u64 / lengths.2 as u64) as u32,
        ];
        charges_vector_pos.push((coord_charge, charge))
    }

    charges_vector_pos
}

#[allow(unused_mut)]
// variables are mutated with deborrow
/// precomputes the electric field from a Vector of charges
pub fn precompute_E(lbm: &Lbm, charges: Vec<(u64, f32)>) {
    // TODO: Make multi-domain compatible

    // Set variables
    let n = lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64;
    let mut e_field: Vec<f32> = vec![0.0; (n * 3) as usize];
    let lengths: (u32, u32, u32) = (lbm.config.n_x, lbm.config.n_y, lbm.config.n_z);
    let def_ke = lbm.config.units.si_to_ke();
    println!(
        "Precomputing electric field for {} charges and {} cells. (This may take a while)",
        charges.len(),
        n
    );

    fn deborrow<'b, T>(r: &T) -> &'b mut T {
        // Needed to access e_field in parallel.
        // This is safe, because no indecies are accessed multiple times
        unsafe {
            #[allow(mutable_transmutes)]
            std::mem::transmute(r)
        }
    }

    let charges_float_pos: Vec<([u32; 3], f32)> = charge_u32_pos(charges, lengths);

    //let mut count: u32 = 0;
    (0..n).into_par_iter().for_each(|i| {
        let e_at = calculate_e(i, &charges_float_pos, lengths, def_ke);
        deborrow(&e_field)[i as usize] = e_at[0];
        deborrow(&e_field)[(i + n) as usize] = e_at[1];
        deborrow(&e_field)[(i + (n * 2)) as usize] = e_at[2];
        //*deborrow(&count) += 1;
        //print!("\r{}", info::progressbar(count as f32 / n as f32)); (This is really slow)
    });
    //println!();

    lbm.domains[0]
        .e
        .as_ref()
        .expect("E buffer used but not initialized")
        .write(&e_field)
        .enq()
        .unwrap();
}

#[inline]
#[allow(unused)]
fn len_sq(v: [f32; 3]) -> f32 {
    v[0].sq() + v[1].sq() + v[2].sq()
}

#[inline]
fn len_sq_u32(v: [u32; 3]) -> f32 {
    (v[0] as f32).sq() + (v[1] as f32).sq() + (v[2] as f32).sq()
}

// Fast vector normalization
#[inline]
#[allow(unused)]
fn fast_normalize(v: [f32; 3]) -> [f32; 3] {
    let len = fast_inv_sqrt(len_sq(v));
    v.map(|x| x * len)
}

// Fast vector normalization
fn fast_normalize_u32(v: [u32; 3]) -> [f32; 3] {
    let len = fast_inv_sqrt(len_sq_u32(v));
    v.map(|x| x as f32 * len)
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
