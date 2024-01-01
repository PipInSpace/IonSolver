
use crate::*;
use rayon::prelude::*;


#[allow(unused_mut)] // Variables are mutated with deborrow
/// Precomputes the electric field from a Vector of charges
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

    let charges_u32_3_pos: Vec<([u32; 3], f32)> = u32_3_pos(charges, lengths);

    (0..n).into_par_iter().for_each(|i| {
        let e_at = calculate_e_at(i, &charges_u32_3_pos, lengths, def_ke);
        deborrow(&e_field)[i as usize] = e_at[0];
        deborrow(&e_field)[(i + n) as usize] = e_at[1];
        deborrow(&e_field)[(i + (n * 2)) as usize] = e_at[2];
    });

    // Write to device
    lbm.domains[0]
        .e
        .as_ref()
        .expect("E buffer used but not initialized")
        .write(&e_field)
        .enq()
        .unwrap();
}

/// Calculates electric field vector at a cell with index n
/// from a vector of charges.
fn calculate_e_at(
    n: u64,
    charges: &[([u32; 3], f32)],
    lengths: (u32, u32, u32),
    def_ke: f32,
) -> [f32; 3] {

    // Compute current cell coordinates
    let coord_cell = coord(n, lengths);
    // Initialize field vector
    let mut e_at_cell = [0.0; 3];

    // Loop over all charges in the simulation
    for &(coord_charge, charge) in charges.iter() {
        // Compute difference vector from cell to current charge
        let coord_diff = [
            (coord_cell[0] as i32) - (coord_charge[0] as i32),
            (coord_cell[1] as i32) - (coord_charge[1] as i32),
            (coord_cell[2] as i32) - (coord_charge[2] as i32),
        ];

        // if not the cell we are checking
        if !coord_diff.eq(&[0 as i32; 3]) {
            let pre_e = charge / cb(length(coord_diff));
            e_at_cell = [
                e_at_cell[0] + coord_diff[0] as f32 * pre_e,
                e_at_cell[1] + coord_diff[1] as f32 * pre_e,
                e_at_cell[2] + coord_diff[2] as f32 * pre_e,
            ];
        }
    }
    // Scale field vector with coulombs constant
    e_at_cell[0] *= def_ke;
    e_at_cell[1] *= def_ke;
    e_at_cell[2] *= def_ke;
    e_at_cell
}

#[allow(unused_mut)] // Variables are mutated with deborrow
/// Precomputes the magnetic field from a Vector of magnetic scalar potentials
pub fn precompute_B(lbm: &Lbm, magnets: Vec<(u64, [f32; 3])>) {
    // TODO: Make multi-domain compatible

    let vec_psi = precompute::calculate_psi_field(&lbm, magnets);

    // Set variables
    let n = lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64;
    let mut b_field: Vec<f32> = vec![0.0; (n * 3) as usize];
    let lengths: (u32, u32, u32) = (lbm.config.n_x, lbm.config.n_y, lbm.config.n_z);
    let def_mu0 = lbm.config.units.si_to_mu0();
    println!(
        "Precomputing magnetic field for {} cells. (This may take a while)",
        n
    );

    fn deborrow<'b, T>(r: &T) -> &'b mut T {
        // Needed to access b_field in parallel.
        // This is safe, because no indecies are accessed multiple times
        unsafe {
            #[allow(mutable_transmutes)]
            std::mem::transmute(r)
        }
    }

    (0..n).into_par_iter().for_each(|i| {
        let b_at = calculate_b_at(i, &vec_psi, lengths, def_mu0);
        deborrow(&b_field)[i as usize] = b_at[0];
        deborrow(&b_field)[(i + n) as usize] = b_at[1];
        deborrow(&b_field)[(i + (n * 2)) as usize] = b_at[2];
    });

    // Write to device
    lbm.domains[0]
        .b
        .as_ref()
        .expect("B buffer used but not initialized")
        .write(&b_field)
        .enq()
        .unwrap();
}

#[allow(unused)]
fn calculate_b_at(
    n: u64,
    psi: &[f32],
    lengths: (u32, u32, u32),
    def_mu0: f32,
) -> [f32; 3] {
    let mut b = [0.0f32; 3];

    let coord = coord(n, lengths);
    let pre_b = nabla(psi, (lengths.0 + 2, lengths.1 + 2, lengths.2 + 2), coord[0] + 1, coord[1] + 1, coord[2] + 1);

    b[0] = -def_mu0 * pre_b[0];
    b[1] = -def_mu0 * pre_b[1];
    b[2] = -def_mu0 * pre_b[2];

    b
}

#[allow(unused_mut)]
/// Calculate a psi field with one cell padding on all sides (Needed for nabla operator in B field)
pub fn calculate_psi_field(
    lbm: &Lbm,
    magnets: Vec<(u64, [f32; 3])>,
) -> Vec<f32> {
    // Set variables
    let lengths: (u32, u32, u32) = (lbm.config.n_x, lbm.config.n_y, lbm.config.n_z);
    let n = ((lengths.0 + 2) as u64) * ((lengths.1 + 2) as u64) * ((lengths.2 + 2) as u64);
    // 1 padding on each side
    let mut psi_field = vec![0.0f32; n as usize];

    println!(
        "Precomputing magnetic scalar potential for {} magnets and {} cells. (This may take a while)",
        magnets.len(),
        n,
    );

    let magnets_u32_3_pos: Vec<([u32; 3], [f32; 3])> = u32_3_pos(magnets, lengths);

    // get psi for all including padding
    psi_field.par_iter_mut().enumerate().for_each(|(i, item)| {
        *item = calculate_psi_at(i as u64, &magnets_u32_3_pos, lengths);
    });

    psi_field
}

/// Calculate psi field at a cell
fn calculate_psi_at(
    n: u64,
    magnets: &[([u32; 3], [f32; 3])],
    lengths: (u32, u32, u32),
) -> f32 {
    // Compute current cell coordinates with padding
    let coord_cell = coord(n, (lengths.0 + 2, lengths.1 + 2, lengths.2 + 2));

    let mut psi_at_cell = 0.0f32;

    // loop over all magnet cells
    for &(coord_magnet, magnetization) in magnets.iter() {
        // Compute difference vector from cell to current magnet
        let coord_diff = [
            (coord_cell[0] as i32) - (coord_magnet[0] as i32 + 1),
            (coord_cell[1] as i32) - (coord_magnet[1] as i32 + 1),
            (coord_cell[2] as i32) - (coord_magnet[2] as i32 + 1),
        ];
        if !coord_diff.eq(&[0 as i32; 3]) {
            let pre_psi = dotp_f32_i32(magnetization, coord_diff) / cb(length(coord_diff));
            psi_at_cell += pre_psi;
        }
    }

    psi_at_cell / (4.0 * PI)
}

/// Converts u64 index positions to u32 coords in a vector
fn u32_3_pos<T: Copy>(gen_vector: Vec<(u64, T)>, lengths: (u32, u32, u32)) -> Vec<([u32; 3], T)> {
    let mut gen_vector_u32_3: Vec<([u32; 3], T)> = Vec::with_capacity(gen_vector.len());

    // Precompute position vectors
    for &(i, gent) in gen_vector.iter() {
        let coord = coord(i, lengths);
        gen_vector_u32_3.push((coord, gent))
    }

    gen_vector_u32_3
}

#[inline]
/// Converts u64 index positions to u32 coords
fn coord(n: u64, lengths: (u32, u32, u32)) -> [u32; 3] {
    [
        (n % lengths.0 as u64) as u32,
        (n / lengths.0 as u64 % lengths.1 as u64) as u32,
        (n / lengths.0 as u64 / lengths.1 as u64) as u32,
    ]
}

#[inline]
fn index_of(x: u32, y: u32, z:u32, lengths: (u32, u32, u32)) -> u64 {
    (x as u64) + (y as u64 + z as u64 * lengths.1 as u64) * (lengths.0 as u64)
}

/// Nabla operator. Returns gradient vector of a f32 field
fn nabla(
    f: &[f32],
    lengths: (u32, u32, u32),
    x: u32,
    y: u32,
    z: u32,
) -> [f32; 3] {
    // assume padding goes to -1, shifting done outside
    let minus_x = index_of(x - 1, y, z, lengths);
    let plus_x = index_of(x + 1, y, z, lengths);
    let minus_y = index_of(x, y - 1, z, lengths);
    let plus_y = index_of(x, y + 1, z, lengths);
    let minus_z = index_of(x, y, z - 1, lengths);
    let plus_z = index_of(x, y, z + 1, lengths);

    let mut grad = [0.0f32; 3];
    grad[0] = (f[plus_x as usize] - f[minus_x as usize]) / 2.0;
    grad[1] = (f[plus_y as usize] - f[minus_y as usize]) / 2.0;
    grad[2] = (f[plus_z as usize] - f[minus_z as usize]) / 2.0;

    grad
}

#[inline]
/// Dot product between f32 and i32 vectors
fn dotp_f32_i32(a: [f32; 3], b: [i32; 3]) -> f32 {
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

#[inline]
fn cb(x: f32) -> f32 {
    x * x * x
}