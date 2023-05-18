use micro_ndarray::vec_split::{RawVector, Vector};

#[derive(Default, Clone, Copy, Debug)]
#[repr(C)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vector<f64, 2> for Vec2 {
    fn get<'a>(&'a self, i: usize) -> Option<&'a f64> {
        match i {
            0 => Some(&self.x),
            1 => Some(&self.y),
            _ => None,
        }
    }

    fn get_mut<'a>(&'a mut self, i: usize) -> Option<&'a mut f64> {
        match i {
            0 => Some(&mut self.x),
            1 => Some(&mut self.y),
            _ => None,
        }
    }
}

// I don't know if this is OK, but it extremely likely is. Remove if weirdness happens.
unsafe impl RawVector<f64, 2> for Vec2 {}

#[cfg(test)]
mod test {
    use micro_ndarray::vec_split::{accessors::IterateAccessorMut, SizedVectorArray};

    use super::Vec2;

    // Test passes, I'll just assume RawVector is ok to implement.
    #[test]
    fn raw_vector() {
        let mut array = [Vec2::default(); 5];
        {
            let [mut x_array, mut y_array] = array.vec_split_fast_mut();
            for item in x_array.iter_mut().chain(y_array.iter_mut()) {
                *item = f64::MAX;
            }
        }
        for item in array {
            println!("{item:?}");
            assert_eq!(item.x, f64::MAX);
            assert_eq!(item.y, f64::MAX);
        }
    }
}

#[allow(dead_code)]
impl Vec2 {
    /// Returns a new Vector with the coordinates specified.
    pub fn new(x: f64, y: f64) -> Vec2 {
        Self { x, y }
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

    /// Returns a new Vector with inverted x component
    pub fn flip_x(&self) -> Vec2 {
        Vec2::new(-self.x, self.y)
    }

    /// Returns a new Vector with inverted y component
    pub fn flip_y(&self) -> Vec2 {
        Vec2::new(self.x, -self.y)
    }

    /// Returns a String with the format "(x,y)".
    pub fn to_string(&self) -> String {
        format!("({},{})", self.x, self.y)
    }
}
