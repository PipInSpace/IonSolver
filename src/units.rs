#![allow(unused)]

/// Unit conversion factors
/// 
/// Lattice unit * factor = SI unit eg. lengthLU * m = length in meters
/// 
/// SI unit / factor = Lattice unit eg. 1 meter / m = x lenghtLU
#[derive(Clone, Copy, Default)]
#[allow(non_snake_case)]
pub struct Units {
    /// meter
    pub m: f32,
    /// kilogram
    pub kg: f32,
    /// second
    pub s: f32,
}

impl Units {
    pub fn new() -> Units {
        Units {
            m: 1.0,
            kg: 1.0,
            s: 1.0,
        }
    }

    pub fn set(
        &mut self,
        lbm_length: f32,
        lbm_velocity: f32,
        lbm_rho: f32,
        si_lenght: f32,
        si_velocity: f32,
        si_rho: f32,
    ) {
        self.m = si_lenght / lbm_length;
        self.kg = si_rho / lbm_rho * cb(self.m);
        self.s = lbm_velocity / si_velocity * self.m;
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

    pub fn force_to_si(&self, f: f32) -> f32 {
        f * (self.kg * self.m / (self.s * self.s))
    }

    pub fn charge_to_si(&self, cl: f32) -> f32 {
        // Unit: A/s, A fac is 1
        cl / self.s
    }

    pub fn mag_flux_to_si(&self, b: f32) -> f32 {
        // b unit is Tesla (T) = V * s / m^2 = ((kg * m^2 / (s^3 * A)) * s) / m^2 = kg / s^2
        b * self.kg / sq(self.s)
    }

    pub fn e_field_to_si(&self, e: f32) -> f32 {
        // E unit: V/m = (kg * m^2 / (s^3 * A)) / m = (m * kg) / s^3
        e * ((self.m * self.kg) / cb(self.s))
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

    pub fn si_to_force(&self, f: f32) -> f32 {
        f / (self.kg * self.m / (self.s * self.s))
    }

    pub fn si_to_nu(&self, nu: f32) -> f32 {
        nu * self.s / (self.m * self.m)
    }

    pub fn si_to_epsilon_0(&self) -> f32 {
        // 8.8541878128E-12 F/m
        // Unit: F/m  ==  A * s / V * m  ==  A * s / (kg*m^2/s^3 * A) * m  ==  s^4 * A^2 / kg * m^3
        // Ampere can be ignored. Conversion factor is 1.
        8.8541878128E-12 / to4(self.s) * (self.kg / cb(self.m))
    }

    pub fn si_to_ke(&self) -> f32 {
        // 1 / (4 * pi * epsilon_0) = k_e (Coulombs constant)
        // epsilon_0 has the unit Farad/meter and needs to be converted to lattice units
        // 1 / (4 * 3.14159 * (8.8541878128 * 10^-12)) = 8.987552E9 (k_e)
        //8.987552E9 / (self.kg * cb(self.m) / (self.c * self.c * self.s * self.s)) -- OLD
        1.0 / (4.0 * std::f32::consts::PI * self.si_to_epsilon_0())
    }

    pub fn si_to_mu_0(&self) -> f32 {
        // 1 / (epsilon_0 * c²) = mu_0 (Magnetic Field Constant)
        //1.25663706212E-6 / (self.kg * self.m / (self.c * self.c)) -- OLD
        1.0 / (self.si_to_epsilon_0() * sq(2.99792458E8) / sq(self.m) * sq(self.s))
    }

    pub fn si_to_charge_per_dens(&self, cpd: f32) -> f32 {
        // unit: (A/s)/(kg/m^3) = m^3/(kg*s)
        cpd * (cb(self.m) * (self.kg / self.s))
    }

    /// From lbm.n_x and velocity u
    pub fn nu_from_Re(&self, Re: f32, x: f32, u: f32) -> f32 {
        x * u / Re
    }

    pub fn print(&self) {
        println!("Units:\n    1 meter = {} lenght LU\n    1 kg = {} dens LU\n    1 second = {} time steps", 1.0/self.m, 1.0/self.kg, 1.0/self.s);
        println!("    epsilon_0 in LU is: {}", self.si_to_epsilon_0());
        println!("    mu_0 in LU is: {}", self.si_to_mu_0());
        println!("    ke in LU is: {}", self.si_to_ke());
    }
}

#[inline]
fn sq(x: f32) -> f32 {
    x * x
}

#[inline]
fn cb(x: f32) -> f32 {
    x * x * x
}

#[inline]
fn to4(x: f32) -> f32 {
    x * x * x * x
}
