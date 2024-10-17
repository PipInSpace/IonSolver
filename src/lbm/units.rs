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
    /// kelvin
    pub k: f32,
    // propellant gas
    pub prop: Propellant,
}

impl Units {
    pub fn new() -> Units {
        Units {
            m: 1.0,
            kg: 1.0,
            s: 1.0,
            a: 1.0,
            k: 1.0,
            prop: Propellant::default(),
        }
    }

    pub fn set(
        &mut self,
        lbm_length: f32,
        lbm_velocity: f32,
        lbm_rho: f32,
        lbm_charge: f32,
        lbm_temp: f32,
        si_length: f32,
        si_velocity: f32,
        si_rho: f32,
        si_charge: f32,
        si_temp: f32,
    ) {
        self.m = si_length / lbm_length;
        self.kg = si_rho / lbm_rho * cb(self.m);
        self.s = self.m / (si_velocity / lbm_velocity);
        self.a = si_charge / lbm_charge / self.s;
        self.k = si_temp / lbm_temp;
    }

    // from lattice units to si units (need to be called after .set();)
    pub fn len_lu_si(&self, l: f32) -> f32 {
        l * self.m
    }

    pub fn mass_lu_si(&self, m: f32) -> f32 {
        m * self.kg
    }

    pub fn dens_lu_si(&self, rho: f32) -> f32 {
        rho * (self.kg / cb(self.m))
    }

    pub fn time_lu_si(&self, t: f32) -> f32 {
        t * self.s
    }

    pub fn speed_lu_si(&self, v: f32) -> f32 {
        v * (self.m / self.s)
    }

    pub fn force_lu_si(&self, f: f32) -> f32 {
        f * (self.kg * (self.m / (self.s * self.s)))
    }

    pub fn charge_lu_si(&self, q: f32) -> f32 {
        // Unit: As
        q * (self.a * self.s)
    }

    pub fn mag_flux_lu_si(&self, b: f32) -> f32 {
        // b unit is Tesla (T) = V * s / m^2 = ((kg * m^2 / (s^3 * A)) * s) / m^2 = kg / (A * s²)
        b * (self.kg / (self.a * sq(self.s)))
    }

    pub fn e_field_lu_si(&self, e: f32) -> f32 {
        // E unit: V/m = (kg * m^2 / (s^3 * A)) / m = (kg * m) / (A * s³)
        e * ((self.kg * self.m) / (self.a * cb(self.s)))
    }

    // from si units to lattice units (need to be called after .set();)
    pub fn len_si_lu(&self, l: f32) -> f32 {
        l / self.m
    }

    pub fn mass_si_lu(&self, m: f32) -> f32 {
        m / self.kg
    }

    pub fn dens_si_lu(&self, rho: f32) -> f32 {
        rho / (self.kg / cb(self.m))
    }

    pub fn time_si_lu(&self, t: f32) -> f32 {
        t / self.s
    }

    pub fn speed_si_lu(&self, v: f32) -> f32 {
        v / (self.m / self.s)
    }

    pub fn force_si_lu(&self, f: f32) -> f32 {
        f / (self.kg * (self.m / sq(self.s)))
    }

    pub fn nu_si_lu(&self, nu: f32) -> f32 {
        nu / (sq(self.m) / self.s)
    }

    pub fn charge_si_lu(&self, q: f32) -> f32 {
        // unit: (A*s)
        q / (self.a * self.s)
    }

    pub fn epsilon_0_lu(&self) -> f32 {
        // 8.8541878128E-12 F/m
        // Unit: F/m  ==  A * s / V * m  ==  A * s / (kg*m^2/s^3 * A) * m  ==  s^4 * A^2 / kg * m^3
        8.8541878128E-12 / ((sq(self.a) * to4(self.s)) / (self.kg * cb(self.m)))
    }

    pub fn ke_lu(&self) -> f32 {
        // 1 / (4 * pi * epsilon_0) = k_e (Coulombs constant)
        // epsilon_0 has the unit Farad/meter and needs to be converted to lattice units
        // 1 / (4 * 3.14159 * (8.8541878128 * 10^-12)) = 8.987552E9 (k_e)
        //8.987552E9 / (self.kg * cb(self.m) / (self.c * self.c * self.s * self.s)) -- OLD
        1.0 / (4.0 * std::f32::consts::PI * self.epsilon_0_lu())
    }

    pub fn mu_0_lu(&self) -> f32 {
        // 1 / (epsilon_0 * c²) = mu_0 (Magnetic Field Constant)
        //1.25663706212E-6 / (self.kg * self.m / (self.c * self.c)) -- OLD
        1.256637062E-6 / ((self.kg * self.m) / (sq(self.a) * sq(self.s)))
    }

    pub fn k_charge_expansion_lu(&self) -> f32 {
        // TODO: q advection uses an expansion coefficient (org. thermal)
        // This value is set to 1.0 in test simulations, the resulting def_w_T is 0.4
        1.0 //todo!()
    }

    pub fn kkge_lu(&self) -> f32 { // Mass per charge constant for electrons
        (9.109_383_713_9E-31_f64/-1.602_176_634E-19_f64) as f32 / (self.kg / (self.a * self.s))
    }

    pub fn kimg_lu(&self) -> f32 { // Inverse of mass of a propellant gas atom, scaled by 10^20
        ((1.0/(self.prop.atom_mass() * 1e20)) / self.kg as f64) as f32
    }

    pub fn kveV_lu(&self) -> f32 { // 9.10938356e-31kg / (2*1.6021766208e-19), velocity to eV for electrons
        (9.109_383_713_9E-31_f64/(2.0_f64 * 1.602_176_634E-19_f64) / self.kg as f64) as f32
    }

    pub fn kkBme_lu(&self) -> f32 { // -1.5 * 1.38064852e-23J/K / 9.10938356e-31kg -- m^2 * s^-2 * K^-1
        (-22734499.72063751808909449412_f64 / (sqd(self.m as f64) / (sqd(self.s as f64) * self.k as f64)) ) as f32
    }

    /// From lbm.n_x and velocity u
    pub fn nu_from_Re(&self, Re: f32, x: f32, u: f32) -> f32 {
        x * u / Re
    }

    pub fn print(&self) {
        println!("Units:\n    1 meter = {} length LU\n    1 kg = {} mass LU\n    1 second = {} time steps\n    1 coulomb = {} charge LU", self.len_si_lu(1.0), self.mass_si_lu(1.0), self.time_si_lu(1.0), self.charge_si_lu(1.0));
        println!("    epsilon_0 in LU is: {}", self.epsilon_0_lu());
        println!("    mu_0 in LU is: {}", self.mu_0_lu());
        println!("    ke in LU is: {}", self.ke_lu());
    }
}

/// Enum holding different gas types and their ionization energies
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub enum Propellant {
    #[default]
    H,
    He,
    Ne,
    Ar,
    Kr,
    Xe,
}

impl Propellant {
    /// Returns the minimal energy value in eV for successfull first ionization
    /// https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
    fn ion_energy(&self) -> f32 {
        match self {
            Propellant::H => 13.598434599702,
            Propellant::He =>  24.587389011,
            Propellant::Ne =>  21.564541,
            Propellant::Ar =>  15.7596119,
            Propellant::Kr => 13.9996055,
            Propellant::Xe =>  12.1298437,
        }
    }

    fn atom_mass(&self) -> f64 {
        match self {
            Propellant::H =>  1.6735575e-27_f64,
            Propellant::He => 6.6464731e-27_f64,
            Propellant::Ne => 3.3509177e-26_f64,
            Propellant::Ar => 6.6335209e-26_f64,
            Propellant::Kr => 1.3914984e-25_f64,
            Propellant::Xe => 2.1801714e-25_f64,
        }
    }
}

#[inline]
fn sq(x: f32) -> f32 {
    x * x
}

#[inline]
fn sqd(x: f64) -> f64 {
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
