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
    lbm_config.units.set(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    lbm_config.n_x = 128;
    lbm_config.n_y = 128;
    lbm_config.n_z = 128;
    lbm_config.d_x = 1;
    lbm_config.nu = 0.1;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    lbm_config.graphics_config.graphics = true;
    lbm_config.graphics_config.streamline_every = 16;
    lbm_config.ext_volume_force = true;
    lbm_config.ext_electric_force = true;
    let mut lbm = Lbm::new(lbm_config);

    // Setup test charges
    let mut vec_q: Vec<(u64, f32)> = vec![];
    //vec_q.push((1056824, 0.00000001));
    //vec_q.push((1056840, -0.00000001));
    
    vec_q.push((1048576+8192, 0.00000001));
    vec_q.push((1048576+8192+127, -0.00000001));
    // Pseudorandom generation 
    //let mut s: u64 = now.elapsed().as_nanos() as u64;
    //for _i in 0..20 {
    //    s = s.wrapping_add(0xA0761D6478BD642F);
    //    let t = u128::from(s) * u128::from(s ^ 0xE7037ED1A0B428DB);
    //    let r = (t as u64) ^ (t >> 64) as u64;
    //    let f = r as f64 / u64::MAX as f64;
    //    vec_q.push(((2097152.0 * f) as u64, -0.0000000001));
    //}
    efield_precompute::precompute_E(&lbm, vec_q);

    //lbm.setup_taylor_green();
    lbm.domains[0]
        .graphics
        .as_mut()
        .expect("grapics not enabled")
        .streamline_mode = true;
    lbm.domains[0].graphics.as_mut().expect("grapics not enabled").vector_e_mode = true;
    lbm.domains[0].graphics.as_mut().expect("grapics not enabled").field_mode = false;
    lbm.domains[0].graphics.as_mut().expect("grapics not enabled").q_mode = false;

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
}
