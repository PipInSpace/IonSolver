//! # units
//! 
//! Holds the `Units` struct and associated unit conversion functions.

#![allow(unused)]

/// Unit conversion factors
///
/// Lattice unit * factor = SI unit eg. lengthLU * m = length in meters
///
/// SI unit / factor = Lattice unit eg. 1 meter / m = x lenghtLU
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
#[allow(non_snake_case)]
pub struct Units {
    /// meter
    pub m: f32,
    /// kilogram
    pub kg: f32,
    /// second
    pub s: f32,
    /// ampere
    pub a: f32,
}

impl Units {
    pub fn new() -> Units {
        Units {
            m: 1.0,
            kg: 1.0,
            s: 1.0,
            a: 1.0,
        }
    }

    pub fn set(
        &mut self,
        lbm_length: f32,
        lbm_velocity: f32,
        lbm_rho: f32,
        lbm_charge: f32,
        si_length: f32,
        si_velocity: f32,
        si_rho: f32,
        si_charge: f32,
    ) {
        self.m = si_length / lbm_length;
        self.kg = si_rho / lbm_rho * cb(self.m);
        self.s = self.m / (si_velocity / lbm_velocity);
        self.a = si_charge / lbm_charge / self.s;
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
        f * (self.kg * (self.m / (self.s * self.s)))
    }

    pub fn charge_to_si(&self, q: f32) -> f32 {
        // Unit: As
        q * (self.a * self.s)
    }

    pub fn mag_flux_to_si(&self, b: f32) -> f32 {
        // b unit is Tesla (T) = V * s / m^2 = ((kg * m^2 / (s^3 * A)) * s) / m^2 = kg / (A * s²)
        b * (self.kg / (self.a * sq(self.s)))
    }

    pub fn e_field_to_si(&self, e: f32) -> f32 {
        // E unit: V/m = (kg * m^2 / (s^3 * A)) / m = (kg * m) / (A * s³)
        e * ((self.kg * self.m) / (self.a * cb(self.s)))
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
        f / (self.kg * (self.m / sq(self.s)))
    }

    pub fn si_to_nu(&self, nu: f32) -> f32 {
        nu / (sq(self.m) / self.s)
    }

    pub fn si_to_epsilon_0(&self) -> f32 {
        // 8.8541878128E-12 F/m
        // Unit: F/m  ==  A * s / V * m  ==  A * s / (kg*m^2/s^3 * A) * m  ==  s^4 * A^2 / kg * m^3
        8.8541878128E-12 / ((sq(self.a) * to4(self.s)) / (self.kg * cb(self.m)))
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
        1.256637062E-6 / ((self.kg * self.m) / (sq(self.a) * sq(self.s)))
    }

    pub fn si_to_charge(&self, q: f32) -> f32 {
        // unit: (A*s)
        q / (self.a * self.s)
    }

    pub fn si_to_k_charge_expansion(&self) -> f32 {
        // TODO: q advection uses an expansion coefficient (org. thermal)
        // This value is set to 1.0 in test simulations, the resulting def_w_T is 0.4
        1.0 //todo!()
    }

    pub fn si_to_kkge(&self) -> f32 { // Mass per charge constant for electrons
        (9.109_383_713_9E-31_f64/-1.602_176_634E-19_f64) as f32 / (self.kg / (self.a * self.s))
    }

    /// From lbm.n_x and velocity u
    pub fn nu_from_Re(&self, Re: f32, x: f32, u: f32) -> f32 {
        x * u / Re
    }

    pub fn print(&self) {
        println!("Units:\n    1 meter = {} length LU\n    1 kg = {} mass LU\n    1 second = {} time steps\n    1 coulomb = {} charge LU", self.si_to_len(1.0), self.si_to_mass(1.0), self.si_to_time(1.0), self.si_to_charge(1.0));
        println!("    epsilon_0 in LU is: {}", self.si_to_epsilon_0());
        println!("    mu_0 in LU is: {}", self.si_to_mu_0());
        println!("    ke in LU is: {}", self.si_to_ke());
    }
}

// Enum holding different gas types and their ionization energies
enum IonE {
    H,
    He,
    Ne,
    Ar,
    Kr,
    Xe,
}

impl IonE {
    /// Returns the minimal energy value in eV for successfull first ionization
    /// https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
    fn val(&self) -> f32 {
        match self {
            IonE::H => 13.598434599702,
            IonE::He =>  24.587389011,
            IonE::Ne =>  21.564541,
            IonE::Ar =>  15.7596119,
            IonE::Kr => 13.9996055,
            IonE::Xe =>  12.1298437,
        }
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
