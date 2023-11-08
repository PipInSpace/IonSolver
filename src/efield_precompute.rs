use std::ops::Mul;

pub fn calculate_e(
    n: usize,
    q: Vec<(usize, [f32; 3])>,
    lengths: (usize, usize, usize),
    def_ke: f32,
) -> [f32; 3] {
    let convert = move |n| -> [f32; 3] {
        [
            (n % lengths.1) as f32,
            (n / lengths.1 % lengths.2) as f32,
            (n / lengths.1 / lengths.2) as f32,
        ]
    };

    let coord = convert(n);
    let mut e_at_cell = [0.0; 3];
    for &(i, charge) in q.iter() {
        let coord_charge = convert(i);
        let coord_diff = [
            coord[0] - coord_charge[0],
            coord[1] - coord_charge[1],
            coord[2] - coord_charge[2],
        ];
        let normalized = normalize(coord_diff);
        let length_sq = len_sq(coord_diff);
        e_at_cell = [
            e_at_cell[0] + charge[0] / length_sq * normalized[0],
            e_at_cell[1] + charge[1] / length_sq * normalized[1],
            e_at_cell[2] + charge[2] / length_sq * normalized[2],
        ]
    }

    e_at_cell.map(|x| x * def_ke)
}

#[inline]
fn len_sq(v: [f32; 3]) -> f32 {
    v[0].sq() + v[1].sq() + v[2].sq()
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = len_sq(v).sqrt();
    v.map(|x| x / len)
}

trait Sq: Copy + Sized + Mul<Self, Output = Self> {
    fn sq(self) -> Self {
        self * self
    }
}

impl<T> Sq for T where T: Copy + Sized + Mul<Self, Output = Self> {}
