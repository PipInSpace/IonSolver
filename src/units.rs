#![allow(unused)]

/// Unit conversion factors (Lattice unit * factor = si unit eg. lengthLU * m = length in meters)
#[derive(Clone, Copy, Default)]
#[allow(non_snake_case)]
pub struct Units {
    /// meter
    m: f32,
    /// kilogram
    kg: f32,
    /// second
    s: f32,
    /// coulomb
    c: f32,
}

impl Units {
    pub fn new() -> Units {
        Units {
            m: 1.0,
            kg: 1.0,
            s: 1.0,
            c: 1.0,
        }
    }

    pub fn set(
        &mut self,
        lbm_length: f32,
        lbm_velocity: f32,
        lbm_rho: f32,
        c: f32,
        si_lenght: f32,
        si_velocity: f32,
        si_rho: f32,
        si_c: f32,
    ) {
        self.m = si_lenght / lbm_length;
        self.kg = si_rho / lbm_rho * cb(self.m);
        self.s = lbm_velocity / si_velocity * self.m;
        self.c = si_c / c;
    }

    // to si units from lattice units (need to be called after .set();)
    pub fn len_to_si(&self, l: f32) -> f32 {
        l * self.m
    }

    pub fn mass_to_si(&self, m: f32) -> f32 {
        m * self.kg
    }

    pub fn dens_to_si(&self, rho: f32) -> f32 {
        rho * (self.kg / cb(self.m))
    }

    pub fn time_to_si(&self, t: f32) -> f32 {
        t * self.s
    }

    pub fn speed_to_si(&self, v: f32) -> f32 {
        v * (self.m / self.s)
    }

    pub fn charge_to_si(&self, q: f32) -> f32 {
        q * self.c
    }

    pub fn force_to_si(&self, f: f32) -> f32 {
        f * (self.kg * self.m / (self.s * self.s))
    }

    // to lattice units from si units (need to be called after .set();)
    pub fn si_to_len(&self, l: f32) -> f32 {
        l / self.m
    }

    pub fn si_to_mass(&self, m: f32) -> f32 {
        m / self.kg
    }

    pub fn si_to_dens(&self, rho: f32) -> f32 {
        rho / (self.kg / cb(self.m))
    }

    pub fn si_to_time(&self, t: f32) -> f32 {
        t / self.s
    }

    pub fn si_to_speed(&self, v: f32) -> f32 {
        v / (self.m / self.s)
    }

    pub fn si_to_charge(&self, q: f32) -> f32 {
        q / self.c
    }

    pub fn si_to_force(&self, f: f32) -> f32 {
        f / (self.kg * self.m / (self.s * self.s))
    }

    pub fn si_to_nu(&self, nu: f32) -> f32 {
        nu * self.s / (self.m * self.m)
    }

    pub fn si_to_ke(&self) -> f32 {
        // 1 / (4 * pi * epsilon_0) = k_e
        // epsilon_0 has the unit Farad/meter and needs to be converted to lattice units
        // 1 / (4 * 3.14159 * (8.8541878128 * 10^-12)) = 8.987552E9 (k_e)
        8.987552E9 / (self.kg * cb(self.m) / (self.c * self.c * self.s * self.s))
    }
}

#[inline]
fn cb(x: f32) -> f32 {
    x * x * x
}
