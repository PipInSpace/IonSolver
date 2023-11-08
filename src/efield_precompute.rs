use std::ops::Mul;

use crate::{lbm::Lbm, info};
use rayon::prelude::*;

fn calculate_e(
    n: u64,
    q: &Vec<(u64, f32)>,
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
    for &(i, charge) in q.iter() {
        let coord_charge = convert(i);
        let coord_diff = [
            coord[0] - coord_charge[0],
            coord[1] - coord_charge[1],
            coord[2] - coord_charge[2],
        ];
        let normalized = normalize(coord_diff);
        let length_sq = len_sq(coord_diff);
        e_at_cell = [
            e_at_cell[0] + charge / length_sq * normalized[0],
            e_at_cell[1] + charge / length_sq * normalized[1],
            e_at_cell[2] + charge / length_sq * normalized[2],
        ]
    }

    e_at_cell.map(|x| x * def_ke)
}

#[allow(unused_mut)]
// variables are mutated with deborrow
pub fn precompute_E(lbm: &Lbm, charges: Vec<(u64, f32)>) {
    println!("Precomputing E for {} charges", charges.len());

    let n = lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64;
    let mut e_field: Vec<f32> = vec![0.0; (n*3) as usize];
    let lenghts: (u32, u32, u32) = (lbm.config.n_x, lbm.config.n_y, lbm.config.n_z);
    let def_ke = lbm.config.units.si_to_ke(8.987552E9);

    fn deborrow<'a, 'b, T>(r: &'a T) -> &'b mut T {
        // This is safe, because no indecies are accessed multiple times
        unsafe { #[allow(mutable_transmutes)] std::mem::transmute(r) }
    }

    let mut count: u32 = 0;

    (0..n).into_par_iter().for_each(|i|{
        let e_at = calculate_e(i, &charges, lenghts, def_ke);
        deborrow(&e_field)[i       as usize] = e_at[0];
        deborrow(&e_field)[(i+n)   as usize] = e_at[1];
        deborrow(&e_field)[(i+n+2) as usize] = e_at[2];
        *deborrow(&count) += 1;
        print!("\r{}", info::progressbar(count as f32 / n as f32));
    });
    println!();

    lbm.domains[0].e.as_ref().expect("E buffer used but not initialized").write(&e_field).enq().unwrap();
}

#[inline]
fn len_sq(v: [f32; 3]) -> f32 {
    v[0].sq() + v[1].sq() + v[2].sq()
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = len_sq(v).sqrt();
    v.map(|x| x / len)
}

trait Sq: Copy + Sized + Mul<Self, Output = Self> {
    fn sq(self) -> Self {
        self * self
    }
}

impl<T> Sq for T where T: Copy + Sized + Mul<Self, Output = Self> {}
