//! Testing utilities

use std::iter;

use itertools::chain;

use crate::{Poly, Poly64, Scalar};

fn binary_coeffs_inner(max_len: usize) -> Box<dyn Iterator<Item = Vec<f64>>> {
    if max_len == 0 {
        return Box::new(iter::once(vec![]));
    }
    Box::new(binary_coeffs_inner(max_len - 1).flat_map(|v| {
        let mut v1 = v.clone();
        let mut v2 = v;
        v1.push(0.0);
        v2.push(1.0);
        [v1, v2].into_iter()
    }))
}

pub fn binary_coeffs(min_degree: i32, max_degree: usize) -> impl Iterator<Item = Poly<f64>> {
    binary_coeffs_inner(max_degree + 1)
        .map(|v| Poly64::from_real_vec(v))
        .filter(move |p| p.degree() >= min_degree)
}
