//! Testing utilities, do not depend on any of these in production!

use std::iter;

use fastrand::Rng;
use itertools::Itertools;
use num::{complex::Complex64, Complex};

use crate::{util::float::f64_make_safe, Poly, Poly64, RealScalar};

use super::float::f64_make_nonzero;

struct RandStreamF64 {
    state: Rng,
}

impl RandStreamF64 {
    fn new(seed: u64) -> Self {
        Self {
            state: Rng::with_seed(seed),
        }
    }
}

impl Iterator for RandStreamF64 {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        // NOTE: I think fastrand::f64 does not emit subnormals, but just in case
        Some(f64_make_nonzero(self.state.f64()))
    }
}

pub struct RandStreamR64 {
    real_stream: RandStreamF64,
    min: f64,
    max: f64,
}

impl RandStreamR64 {
    #[must_use]
    pub fn new(seed: u64, min: f64, max: f64) -> Self {
        assert!(min <= max, "minimum should be smaller or equal to maximum");
        let real_stream = RandStreamF64::new(seed);
        Self {
            real_stream,
            min,
            max,
        }
    }
}

impl Iterator for RandStreamR64 {
    type Item = Complex64;

    fn next(&mut self) -> Option<Self::Item> {
        let re = (self.real_stream.next()?).mul_add(self.max - self.min, self.min);
        let im = 0.0;
        Some(Complex64 {
            re: f64_make_safe(re),
            im,
        })
    }
}

pub struct RandStreamC64Cartesian {
    real_stream: RandStreamF64,
    min_re: f64,
    max_re: f64,
    min_im: f64,
    max_im: f64,
}

impl RandStreamC64Cartesian {
    #[must_use]
    pub fn new(seed: u64, min_re: f64, max_re: f64, min_im: f64, max_im: f64) -> Self {
        assert!(
            min_re <= max_re && min_im <= max_im,
            "minimum should be smaller or equal to maximum"
        );
        let real_stream = RandStreamF64::new(seed);
        Self {
            real_stream,
            min_re,
            max_re,
            min_im,
            max_im,
        }
    }
}

impl Iterator for RandStreamC64Cartesian {
    type Item = Complex64;

    fn next(&mut self) -> Option<Self::Item> {
        let re = (self.real_stream.next()?).mul_add(self.max_re - self.min_re, self.min_re);
        let im = (self.real_stream.next()?).mul_add(self.max_im - self.min_im, self.min_im);
        Some(Complex::new(f64_make_safe(re), f64_make_safe(im)))
    }
}

pub struct RandStreamC64Polar {
    real_stream: RandStreamF64,
    min_radius: f64,
    max_radius: f64,
    min_angle: f64,
    max_angle: f64,
}

impl RandStreamC64Polar {
    #[must_use]
    pub fn new(
        seed: u64,
        min_radius: f64,
        max_radius: f64,
        min_angle: f64,
        max_angle: f64,
    ) -> Self {
        assert!(
            0.0 <= min_angle && max_angle <= 1.0,
            "angles should be specified in the range [0,1]"
        );
        assert!(
            min_angle <= max_angle,
            "min_angle should be smaller or equal to max_angle"
        );
        assert!(0.0 <= min_radius, "radius should be non-negative");
        assert!(
            min_radius <= max_radius,
            "min_radius should be smaller or equal to max_radius"
        );
        Self {
            real_stream: RandStreamF64::new(seed),
            min_radius,
            max_radius,
            min_angle,
            max_angle,
        }
    }
}

impl Iterator for RandStreamC64Polar {
    type Item = Complex64;

    fn next(&mut self) -> Option<Self::Item> {
        let r =
            (self.real_stream.next()?).mul_add(self.max_radius - self.min_radius, self.min_radius);
        let a = (self.real_stream.next()?).mul_add(self.max_angle - self.min_angle, self.min_angle);
        debug_assert!(r >= 0.0);
        debug_assert!((0.0..=1.0).contains(&a));
        let c = Complex::from_polar(r, a * std::f64::consts::TAU);
        Some(Complex::new(f64_make_safe(c.re), f64_make_safe(c.im)))
    }
}

pub struct RandStreamConjugate64<I: Iterator<Item = Complex<f64>>> {
    upstream: I,
}

impl<I: Iterator<Item = Complex<f64>>> RandStreamConjugate64<I> {
    pub const fn new(upstream: I) -> Self {
        Self { upstream }
    }
}

impl<I: Iterator<Item = Complex<f64>>> Iterator for RandStreamConjugate64<I> {
    type Item = (Complex64, Complex64);

    fn next(&mut self) -> Option<Self::Item> {
        let c = self.upstream.next()?;
        Some((c, c.conj()))
    }
}

pub struct PolyStream<T: RealScalar> {
    max_degree: usize,
    root_stream: Box<dyn Iterator<Item = Complex<T>>>,
}

impl<T: RealScalar> PolyStream<T> {
    pub fn new(max_degree: usize, root_stream: impl Iterator<Item = Complex<T>> + 'static) -> Self {
        Self {
            max_degree,
            root_stream: Box::new(root_stream),
        }
    }
}

impl<T: RealScalar> Iterator for PolyStream<T> {
    type Item = (Vec<Complex<T>>, Poly<T>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut roots = vec![];
        for _ in 0..self.max_degree {
            roots.push(self.root_stream.next()?);
        }
        let poly = Poly::from_roots(&roots);
        Some((roots, poly))
    }
}

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

pub fn binary_coeffs(min_degree: usize, max_degree: usize) -> impl Iterator<Item = Poly<f64>> {
    binary_coeffs_inner(max_degree + 1)
        .map(Poly64::from_real_vec)
        .filter(move |p| p.degree_raw() >= min_degree)
}

/// Generate one test case where the roots are known and can be compared
pub fn test_case_roots(
    roots_stream: impl Iterator<Item = Complex64>,
    mut scale_stream: impl Iterator<Item = Complex64>,
    degree: usize,
) -> (Poly64, Vec<Complex64>) {
    let roots = roots_stream.take(degree).collect_vec();
    let poly = Poly64::from_roots(&roots)
        .scaled(&scale_stream.next().expect("rng stream should be infinite"));
    (poly, roots)
}

/// Generate one test case where the roots are known and can be compared, this
/// makes conjugate roots.
pub fn test_case_conj_roots(
    roots_stream: impl Iterator<Item = Complex64>,
    mut scale_stream: impl Iterator<Item = Complex64>,
    degree: usize,
) -> (Poly64, Vec<Complex64>) {
    let roots_stream = RandStreamConjugate64::new(roots_stream);
    let roots = roots_stream
        .take((degree + 1) / 2)
        .flat_map(|(r1, r2)| [r1, r2])
        .collect_vec();
    let poly = Poly64::from_roots(&roots)
        .scaled(&scale_stream.next().expect("rng stream should be infinite"));
    (poly, roots)
}

pub fn test_case_multiple_roots(
    roots_stream: impl Iterator<Item = Complex64>,
    mut scale_stream: impl Iterator<Item = Complex64>,
    degree: usize,
    multiplicity: usize,
) -> (Poly64, Vec<Complex64>) {
    let mut roots = roots_stream.take(degree - multiplicity).collect_vec();
    let first_root = roots[0];
    for _ in 0..multiplicity {
        roots.push(first_root);
    }
    let poly = Poly64::from_roots(&roots)
        .scaled(&scale_stream.next().expect("rng stream should be infinite"));
    (poly, roots)
}

/// Check that all roots have been found
#[must_use]
pub fn check_roots(roots1: Vec<Complex64>, mut roots2: Vec<Complex64>, tol: f64) -> bool {
    if roots1.len() != roots2.len() {
        return false;
    }

    for r1 in roots1 {
        let mut best_idx = 0;
        let mut best_d = f64::MAX;
        for (i, r2) in roots2.iter().enumerate() {
            let d = (r1 - r2).norm();
            if d < best_d {
                best_idx = i;
                best_d = d;
            }
        }
        if best_d > tol {
            return false;
        }
        roots2.remove(best_idx);
    }
    true
}
