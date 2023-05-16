pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    /// Returns a new Vector with the coordinates specified.
    pub fn new(x: f64, y: f64) -> Vec2 {
        Self { x: x, y: y }
    }

    /// Returns a Vector on which the function is called + the Vector specified.
    pub fn add(&self, a: &Vec2) -> Vec2 {
        Vec2::new(self.x + a.x, self.y + a.y)
    }

    /// Returns a scaled Vector on which the function is called.
    pub fn scale(&self, scalar: f64) -> Vec2 {
        Vec2::new(self.x * scalar, self.y * scalar)
    }

    /// Returns a negated Vector on which the function is called.
    pub fn negate(&self) -> Vec2 {
        self.scale(-1.0)
    }

    /// Returns the magnitude of the Vector on which the function is called.
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Returns the distance from the Vector on which the function is called to the Vector specified.
    pub fn distance(&self, to: &Vec2) -> f64 {
        to.add(&self.negate()).magnitude()
    }

    /// Returns a new normalized Vector of the one on which the function is called.
    pub fn normalize(&self) -> Vec2 {
        self.scale(1.0 / self.magnitude())
    }

    /// Returns a new normalized Vector pointing from the one on which the function is called to the Vector specified.
    pub fn direction(&self, to: &Vec2) -> Vec2 {
        to.add(&self.negate()).normalize()
    }

    /// Returns a String with the format "(x,y)".
    pub fn to_string(&self) -> String {
        format!("({},{})", self.x, self.y)
    }
}
