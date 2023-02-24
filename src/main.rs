use micro_ndarray::Array;


type Vec2 = (f64, f64);

struct Engine {
    particles: Array<f64, 2>,
    velocity: Array<Vec2, 2>,
    gravity: f64,
}

impl Engine {
    pub fn step(&mut self) {
        
    }
}

fn main() {
    Engine {
        particles: Array::new([200, 200]),
        velocity: Array::new([200, 200]),
        gravity: 1.0,
    }.step();
}
