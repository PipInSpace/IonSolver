use crate::*;
use rayon::prelude::*;

/*
lbm.domains[0]
        .e
        .as_ref()
        .expect("E buffer used but not initialized")
        .write(&e_field)
        .enq()
        .unwrap();
*/

// source: https://www.db-thueringen.de/servlets/MCRFileNodeServlet/dbt_derivate_00017776/IWK_2006_3_3_11.pdf

fn calculate_b(
    n: u64,
    magnets: &[([u32; 3], [f32; 3])],
    lengths: (u32, u32, u32),
    def_mu0: f32,
) -> [f32; 3] {

    // Compute current cell coordinates
    let coord_cell = [
        (n % lengths.1 as u64) as u32,
        (n / lengths.1 as u64 % lengths.2 as u64) as u32,
        (n / lengths.1 as u64 / lengths.2 as u64) as u32,
    ];

    let psi_at_cell = [0.0f32; 3];

    // loop over all magnet cells
    for &(coord_magnet, magnetization) in magnets.iter() {
        // Compute difference vector from cell to current magnet
        let coord_diff = [
            (coord_cell[0] as i32) - (coord_magnet[0] as i32),
            (coord_cell[1] as i32) - (coord_magnet[1] as i32),
            (coord_cell[2] as i32) - (coord_magnet[2] as i32),
        ];
        if !coord_diff.eq(&[0 as i32; 3]) {
            let pre_psi = dotp_f32_i32(magnetization, coord_diff) / cb(length(coord_diff));
            psi_at_cell = [
                psi_at_cell[0] + coord_diff[0] as f32 * pre_psi,
                psi_at_cell[1] + coord_diff[1] as f32 * pre_psi,
                psi_at_cell[2] + coord_diff[2] as f32 * pre_psi,
            ];
        }
    }
    
    // TODO: finish b field

    return [0.0; 3];
}

#[inline]
fn dotp_f32_i32(a: [f32; 3], b: [i32; 3]) {
    a[0] * (b[0] as f32) + a[1] * (b[1] as f32) + a[2] * (b[2] as f32)
}

#[inline]
/// Returns the lenght of a u32 vector with all components squared
fn len_sq_i32(v: [i32; 3]) -> f32 {
    sq(v[0] as f32) + sq(v[1] as f32) + sq(v[2] as f32)
}

#[inline]
fn length(v: [i32; 3]) -> f32 {
    f32::sqrt(len_sq_i32(v))
}

#[inline]
/// square
fn sq(x: f32) -> f32 {
    x * x
}
