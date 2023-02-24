use micro_ndarray::Array;

type Vec2 = (f64, f64);

struct Engine {
    u: Array<f64, 2>,
    velo: Array<f64, 2>,
    u_prev: Array<f64, 2>,
    velo_prev: Array<f64, 2>,
    dens: Array<f64, 2>,
    dens_prev: Array<f64, 2>,
}

impl Engine {
    pub fn add_source(&mut self, s: Array<f64, 2>, dt: f64) {
        for ([x, y], item) in self.dens.iter_mut() {
            *item += dt * s[[x, y]];
        }
    }

    pub fn diffuse(&mut self, diff: f64, dt: f64) {
        let a = dt * diff * (self.dens.size().iter().map(|x| *x as f64).product::<f64>());
        let d = self.dens.clone();
        for k in 0..20 {
            for (([x, y], item), (_, item_old)) in self.dens.iter_mut().zip(self.dens_prev.iter()) {
                *item = *item_old
                    + a * (d[[x - 1, y]] + d[[x + 1, y]] + d[[x, y - 1]] + d[[x, y + 1]])
                        / (1.0 + 4.0 * a);
            }
        }
    }

    pub fn step_density(&mut self) {}

    pub fn step(&mut self) {}
}

fn main() {}
