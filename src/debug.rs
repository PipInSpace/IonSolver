use std::cmp::Ordering;

use micro_ndarray::Array;

pub fn print_maxval(x: &Array<f64, 2>, name: &'static str) {
    println!(
        "Max {name}: {:?}",
        x.iter()
            .map(|x| *x.1)
            .max_by(|a, b| f64::partial_cmp(a, b).unwrap_or(Ordering::Equal))
    );
}

/// Sums all elements of the provided f64 array
pub fn print_sum(x: &Array<f64, 2>, name: &'static str) {
    println!("Sum {name}: {:?}", x.iter().map(|x| *x.1).sum::<f64>());
}
