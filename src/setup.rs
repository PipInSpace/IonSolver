use crate::{lbm::Lbm, *};

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
pub fn setup() -> (Lbm, LbmConfig) {
    let mut lbm_config = LbmConfig::new();
    lbm_config.units.set(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    lbm_config.n_x = 128;
    lbm_config.n_y = 128;
    lbm_config.n_z = 128;
    lbm_config.d_x = 1;
    lbm_config.nu = 0.1;
    lbm_config.velocity_set = VelocitySet::D3Q19;
    lbm_config.ext_volume_force = true;
    lbm_config.fx = 0.0001;
    let mut lbm = Lbm::new(lbm_config);
    lbm.setup_taylor_green();
    lbm.domains[0].graphics.streamline_mode = true;
    lbm.domains[0].graphics.q_mode = true;
    lbm.initialize();
    (lbm, lbm_config)
}

impl Lbm {
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
            #[allow(unused_mut)]
            let mut domain_vec_q: Vec<f32> = vec![0.0; (dsx * dsy * dsz) as usize]; // temporary q (all 0)
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
            if self.config.ext_electric_force {
                // Only write to buffer if it is needed/initialized
                self.domains[d as usize]
                    .q
                    .as_mut()
                    .expect("q buffer used but not initialized")
                    .write(&domain_vec_q)
                    .enq()
                    .unwrap();
            }
            self.domains[d as usize].queue.finish().unwrap();
        }
        println!();
        println!("Finished setting up Taylor-Green vorticies");
    }
}
