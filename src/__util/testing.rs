//! Testing utilities

use std::{fs, iter, path::Path};

use fastrand::Rng;
use itertools::{chain, Itertools};
use num::{complex::Complex64, Complex};
use plotters::style::SizeDesc;

use crate::{Poly, Poly64, Scalar, __util::float::f64_make_safe};

use super::float::{f64_make_nonzero, F32_BIG_NUM, F64_BIG_NUM};

struct RandStreamBool {
    state: Rng,
}

impl RandStreamBool {
    fn new(seed: u64) -> Self {
        Self {
            state: Rng::with_seed(seed),
        }
    }
}

impl Iterator for RandStreamBool {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.state.bool())
    }
}

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

struct RandStreamC64Cartesian {
    real_stream: RandStreamF64,
    min_re: f64,
    max_re: f64,
    min_im: f64,
    max_im: f64,
}

impl RandStreamC64Cartesian {
    fn new(seed: u64, min_re: f64, max_re: f64, min_im: f64, max_im: f64) -> Self {
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
        let re = self.real_stream.next()? * (self.max_re - self.min_re) + self.min_re;
        let im = self.real_stream.next()? * (self.max_im - self.min_im) + self.min_im;
        Some(Complex::new(f64_make_safe(re), f64_make_safe(im)))
    }
}

struct RandStreamC64Polar {
    real_stream: RandStreamF64,
    min_radius: f64,
    max_radius: f64,
    min_angle: f64,
    max_angle: f64,
}

impl RandStreamC64Polar {
    fn new(seed: u64, min_radius: f64, max_radius: f64, min_angle: f64, max_angle: f64) -> Self {
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
        let r = self.real_stream.next()? * (self.max_radius - self.min_radius) + self.min_radius;
        let a = self.real_stream.next()? * (self.max_angle - self.min_angle) + self.min_angle;
        debug_assert!(r >= 0.0);
        debug_assert!(0.0 <= a && a <= 1.0);
        let c = Complex::from_polar(r, a * std::f64::consts::TAU);
        Some(Complex::new(f64_make_safe(c.re), f64_make_safe(c.im)))
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

pub fn binary_coeffs(min_degree: i32, max_degree: usize) -> impl Iterator<Item = Poly<f64>> {
    binary_coeffs_inner(max_degree + 1)
        .map(|v| Poly64::from_real_vec(v))
        .filter(move |p| p.degree() >= min_degree)
}

struct PolyStream<T: Scalar> {
    max_degree: usize,
    root_stream: Box<dyn Iterator<Item = Complex<T>>>,
}

impl<T: Scalar> PolyStream<T> {
    fn new(max_degree: usize, root_stream: impl Iterator<Item = Complex<T>> + 'static) -> Self {
        Self {
            max_degree,
            root_stream: Box::new(root_stream),
        }
    }
}

impl<T: Scalar + PartialOrd> Iterator for PolyStream<T> {
    type Item = (Vec<Complex<T>>, Poly<T>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut roots = vec![];
        for _ in 0..self.max_degree {
            roots.push(self.root_stream.next()?)
        }
        let poly = Poly::from_roots(&roots);
        Some((roots, poly))
    }
}

/// Make polynomials from random complex roots, sampled uniformly in the span
/// `[-span, span]`, without subnormals or zeros.
///
/// The constant is also random
pub fn random_uniform_non_zero_roots(
    seed: u64,
    max_degree: usize,
    min_re: f64,
    max_re: f64,
    min_im: f64,
    max_im: f64,
) -> impl Iterator<Item = (Vec<Complex64>, Poly<f64>)> {
    let root_stream = RandStreamC64Cartesian::new(seed, min_re, max_re, min_im, max_im);
    PolyStream::new(max_degree, root_stream)
}

/// Run this manually to visually inspect the random generators
#[test]
#[ignore]
fn plot_random_roots() -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let file_id: String = chain!(
        "temp/".chars(),
        iter::from_fn(|| Some(fastrand::alphanumeric())).take(22),
        ".png".chars()
    )
    .collect();

    // if fails it already exists and that's fine
    let _ = fs::create_dir("./temp/");

    let root = BitMapBackend::new(&file_id, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root).build_cartesian_2d(-2.0..2.0, -2.0..2.0)?;
    chart.draw_series(
        RandStreamC64Polar::new(1, 0.5, 1.5, 0.125, 0.5)
            .take(1000)
            .map(|z| Circle::new((z.re, z.im), 1, BLACK.filled())),
    )?;
    root.present()?;
    Ok(())
}
