use crate::*;

/// `setup()` is called at simulation start. Edit this function to change simulation parameters.
///
/// Usage:
///
/// Start by requesting a new `LbmConfig` with `LbmConfig::new()`.
/// You can then set individual arguments by setting fields of the `LbmConfig`.
/// If you are doing something that requires real-world units/scales, you can define them using `your_lbm_config.units.set()` and pass in your desired units.
///
/// After setting your config struct, request a new `Lbm` struct with `Lbm::new(your_lbm_config)`.
/// Domain setup is handled automatically, you might now directly set specific cells in your domains, for an example look at `setup_taylor_green()`.
///
/// Run `your_lbm.initialize()` and return it with the config.
pub fn setup() -> Lbm {
    //let now = Instant::now();
    let mut lbm_config = LbmConfig::new();
    lbm_config.units.print();
    lbm_config.units.set(128.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.2250, 1.0);
    lbm_config.units.print();
    lbm_config.n_x = 128;
    lbm_config.n_y = 128;
    lbm_config.n_z = 256;
    lbm_config.d_x = 1;
    lbm_config.nu = lbm_config.units.si_to_nu(1.48E-5);
    println!("    nu in LU is: {}", lbm_config.units.si_to_nu(1.48E-3));
    lbm_config.velocity_set = VelocitySet::D3Q19;
    // Extensions
    lbm_config.ext_volume_force = true;
    lbm_config.ext_electro_hydro = true;
    // Graphics
    lbm_config.graphics_config.graphics_active = true;
    lbm_config.graphics_config.background_color = 0x1c1b22;
    lbm_config.graphics_config.streamline_every = 4;
    lbm_config.graphics_config.vec_vis_mode = graphics::VecVisMode::U;
    lbm_config.graphics_config.streamline_mode = true;
    lbm_config.graphics_config.flags_surface_mode = true;
    lbm_config.graphics_config.flags_mode = true;

    //let mut lbm = setup_from_file("b-field_spin.ion", lbm_config);

    let mut lbm = Lbm::new(lbm_config);

    let mut flags: Vec<u8> = vec![0; 128 * 128 * 256];
    let len = flags.len() - 1;
    // Solid
    for i in 0..(128 * 128 * 8) {
        flags[i] = 0x01;
        flags[len - i] = 0x01;
    }
    lbm.domains[0].flags.write(&flags).enq().unwrap();

    //// electric field
    //let mut vec_q: Vec<(u64, f32)> = vec![];
    //vec_q.push((1056824, 0.000000001));
    //vec_q.push((1056840, -0.000000001));
    //precompute::precompute_E(&lbm, vec_q);

    // magnetic field
    let mut vec_m: Vec<(u64, [f32; 3])> = vec![];
    //vec_m.push((1056824, [1.0, 0.0, 0.0]));
    //vec_m.push((1056832, [1.0, 0.0, 0.0]));
    //vec_m.push((1056833, [1.0, 0.0, 0.0]));
    //vec_m.push((1056840, [1.0, 0.0, 0.0]));
    for i in 0..243 {
        vec_m.push((i * 68, [0.0, 0.0, 2000000.0]));
        vec_m.push((i * 68 + 2097152 * 2, [0.0, 0.0, 2000000.0]));
    }
    precompute::precompute_B(&lbm, vec_m);

    lbm.setup_velocity_field((0.1, 0.01, 0.0), 0.5);

    lbm
}

#[allow(unused)]
/// Set up Lbm from file. Requires LbmConfig for additional customization at compilation (Graphics etc.).
/// Use in setup();
///
/// Format documented under https://github.com/PipInSpace/ionsolver-files/blob/main/src/FILE.txt
pub fn setup_from_file(path: &str, lbm_config: LbmConfig) -> Lbm {
    println!("Parsing from \"{}\"", path);
    let buffer: Vec<u8> = std::fs::read(path).expect("Location should exist");
    let vals = parse::decode(&buffer).unwrap();

    // Extension is disabled when not needed
    let electro_hydro = vals.charges.len() > 0 || vals.magnets.len() > 0;

    let mut lbm_config = LbmConfig {
        n_x: vals.n_x,
        n_y: vals.n_y,
        n_z: vals.n_z,
        units: units::Units {
            m: vals.units_m,
            kg: vals.units_kg,
            s: vals.units_s,
            c: vals.units_c,
        },
        ext_volume_force: electro_hydro,
        ext_electro_hydro: electro_hydro,
        ..lbm_config
    };
    lbm_config.units.set(128.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.225, 1.0);

    let lbm = Lbm::new(lbm_config);
    if vals.charges.len() > 0 {
        precompute::precompute_E(&lbm, vals.charges);
    }
    if vals.magnets.len() > 0 {
        precompute::precompute_B(&lbm, vals.magnets);
    }

    lbm.domains[0].flags.write(&vals.flags).enq().unwrap();

    lbm
}

impl Lbm {
    #[allow(unused)]
    /// 3D Taylor-Green vorticies setup.
    pub fn setup_taylor_green(&mut self) {
        println!("Setting up Taylor-Green vorticies");
        let nx = self.config.n_x;
        let ny = self.config.n_y;
        let nz = self.config.n_z;

        let pif = std::f32::consts::PI;
        #[allow(non_snake_case)]
        let A = 0.25f32;
        let periodicity = 2u32;
        let a = nx as f32 / periodicity as f32;
        let b = ny as f32 / periodicity as f32;
        let c = nz as f32 / periodicity as f32;

        let domain_numbers: u32 = self.config.d_x * self.config.d_y * self.config.d_z;
        let dx = self.config.d_x;
        let dy = self.config.d_y;
        let dz = self.config.d_z;
        let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64 * 2; // Domain size on each axis
        let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64 * 2; // Needs to account for halo offsets
        let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64 * 2;
        let dtotal = dsx * dsy * dsz;

        print!("");
        for d in 0..domain_numbers {
            let x = (d % (dx * dy)) % dx; // Current Domain coordinates
            let y = (d % (dx * dy)) / dx;
            let z = d / (dx * dy);

            let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
            let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
            for zi in 0..dsz {
                // iterates over every cell in the domain, filling it with  Taylor-Green-vortex
                print!(
                    "\r{}",
                    info::progressbar(
                        ((zi as f32 + 1.0) / dsz as f32) * (1.0 / domain_numbers as f32)
                            + ((d as f32) / domain_numbers as f32)
                    )
                );
                for yi in 0..dsy {
                    for xi in 0..dsx {
                        if !(((xi == 0 || xi == dsx - 1) && dx > 1)
                            || ((yi == 0 || yi == dsy - 1) && dy > 1)
                            || ((zi == 0 || zi == dsz - 1) && dz > 1))
                        {
                            //do not set at halo offsets
                            let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                            let gx =
                                xi - (dx > 1u32) as u64 + x as u64 * (dsx - (dx > 1u32) as u64 * 2); // Global coordinates
                            let gy =
                                yi - (dy > 1u32) as u64 + y as u64 * (dsy - (dy > 1u32) as u64 * 2);
                            let gz =
                                zi - (dz > 1u32) as u64 + z as u64 * (dsz - (dz > 1u32) as u64 * 2);

                            let fx = gx as f32 + 0.5 - 0.5 * nx as f32;
                            let fy = gy as f32 + 0.5 - 0.5 * ny as f32;
                            let fz = gz as f32 + 0.5 - 0.5 * nz as f32;

                            domain_vec_u[(dn) as usize] = A
                                * (2.0 * pif * fx / a).cos()
                                * (2.0 * pif * fy / b).sin()
                                * (2.0 * pif * fz / c).sin(); // x
                            domain_vec_u[(dn + dtotal) as usize] = -A
                                * (2.0 * pif * fx / a).sin()
                                * (2.0 * pif * fy / b).cos()
                                * (2.0 * pif * fz / c).sin(); // y;
                            domain_vec_u[(dn + dtotal * 2) as usize] = A
                                * (2.0 * pif * fx / a).sin()
                                * (2.0 * pif * fy / b).sin()
                                * (2.0 * pif * fz / c).cos(); // z
                            domain_vec_rho[(dn) as usize] = 1.0
                                - (A * A) * 3.0 / 4.0 * (4.0 * pif * fx / a).cos()
                                + (4.0 * pif * fy / b).cos();
                        }
                    }
                }
            }

            // Write to domain buffers
            self.domains[d as usize]
                .u
                .write(&domain_vec_u)
                .enq()
                .unwrap();
            self.domains[d as usize]
                .rho
                .write(&domain_vec_rho)
                .enq()
                .unwrap();
            self.domains[d as usize].queue.finish().unwrap();
        }
        println!();
        println!("Finished setting up Taylor-Green vorticies");
    }

    #[allow(unused)]
    /// Sets all fluid cells to the specified velocity
    fn setup_velocity_field(&mut self, velocity: (f32, f32, f32), density: f32) {
        println!("Setting up velocity field");
        let nx = self.config.n_x;
        let ny = self.config.n_y;
        let nz = self.config.n_z;
        let domain_numbers: u32 = self.config.d_x * self.config.d_y * self.config.d_z;
        let dx = self.config.d_x;
        let dy = self.config.d_y;
        let dz = self.config.d_z;
        let dsx = nx as u64 / dx as u64 + (dx > 1u32) as u64 * 2; // Domain size on each axis
        let dsy = ny as u64 / dy as u64 + (dy > 1u32) as u64 * 2; // Needs to account for halo offsets
        let dsz = nz as u64 / dz as u64 + (dz > 1u32) as u64 * 2;
        let dtotal = dsx * dsy * dsz;

        print!("");
        for d in 0..domain_numbers {
            let mut domain_vec_u: Vec<f32> = vec![0.0; (dsx * dsy * dsz * 3) as usize];
            let mut domain_vec_rho: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize];
            for zi in 0..dsz {
                // iterates over every cell in the domain, filling it with the velocity field
                print!(
                    "\r{}",
                    info::progressbar(
                        ((zi as f32 + 1.0) / dsz as f32) * (1.0 / domain_numbers as f32)
                            + ((d as f32) / domain_numbers as f32)
                    )
                );
                for yi in 0..dsy {
                    for xi in 0..dsx {
                        if !(((xi == 0 || xi == dsx - 1) && dx > 1)
                            || ((yi == 0 || yi == dsy - 1) && dy > 1)
                            || ((zi == 0 || zi == dsz - 1) && dz > 1))
                        {
                            //do not set at halo offsets
                            let dn = (zi * dsx * dsy) + (yi * dsx) + xi; // Domain 1D index
                            domain_vec_u[(dn) as usize] = velocity.0; // x
                            domain_vec_u[(dn + dtotal) as usize] = velocity.1; // y;
                            domain_vec_u[(dn + dtotal * 2) as usize] = velocity.2; // z
                            domain_vec_rho[(dn) as usize] = density;
                        }
                    }
                }
            }

            // Write to domain buffers
            self.domains[d as usize]
                .u
                .write(&domain_vec_u)
                .enq()
                .unwrap();
            self.domains[d as usize]
                .rho
                .write(&domain_vec_rho)
                .enq()
                .unwrap();
            self.domains[d as usize].queue.finish().unwrap();
        }
        println!();
        println!("Finished setting up velocity field");
    }
}
