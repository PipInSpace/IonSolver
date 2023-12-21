#![allow(unused)]

/// SI Base units (Lattice unit * unit = si unit eg. lengthLU * m = length in meters)
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
        x: f32,
        u: f32,
        rho: f32,
        c: f32,
        si_x: f32,
        si_u: f32,
        si_rho: f32,
        si_c: f32,
    ) {
        self.m = si_x / x;
        self.kg = si_rho / rho * self.m * self.m * self.m;
        self.s = u / si_u * self.m;
        self.c = si_c / c;
    }

    // to si from LU
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

    // to LU from si

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

    pub fn si_to_ke(&self) -> f32 {
        8.987552E9 / (self.kg * cb(self.m) / (self.c * self.c * self.s * self.s))
    }
}

#[inline]
fn cb(x: f32) -> f32 {
    x * x * x
}
