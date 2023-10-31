/// SI Base units
pub struct Units {
    /// meter
    pub m : f32,
    /// kilogram
    pub kg : f32,
    /// second
    pub s : f32,
    /// coulomb
    pub c : f32,
    /// Newton
    pub f : f32,
    /// coulomb's constant
    pub ke : f32,
}

impl Units {
    pub fn new() -> Units {
        Units {
            m : 1.0,
            kg : 1.0,
            s : 1.0,
            c : 1.0,
            f : 1.0,
            ke : 1.0,
        }
    }

    pub fn set(&mut self, x : f32, u : f32, rho : f32, c : f32, si_x : f32, si_u : f32, si_rho : f32, si_c : f32) {
        self.m = si_x / x;
        self.kg = si_rho / rho * self.m * self.m * self.m;
        self.s = u / si_u * self.m;
        self.c = si_c / c;

        self.f = self.kg * self.m / (self.s * self.s);
        self.ke = self.f * self.m * self.m / (self.c * self.c);
    }
}


