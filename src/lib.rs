// TODO(version: v1.0.0): license/author header project-wide, see MIT guidelines
#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]

extern crate nalgebra as na;
use std::ops::Index;

use na::Normed;
pub use num;

use num::{Complex, One, Zero};

/// A more convenient way to write `Complex::new(...)`.
///
/// # Examples
///
/// ```
/// use rust_poly::complex;
/// use num::Complex;
///
/// let c1: Complex<f32> = complex!();
/// let c2 = Complex::new(0.0, 0.0);
/// let c3 = complex!(1.0f32, 2.0);
/// let c4 = Complex::new(1.0, 2.0);
///
/// assert_eq!(c1, c2);
/// assert_eq!(c3, c4);
/// assert_eq!(complex!(4.20), complex!(4.20, 0.0));
/// ```
#[macro_export]
macro_rules! complex {
    () => {{
        <$crate::num::Complex<_> as $crate::num::Zero>::zero()
    }};
    ($re:expr) => {{
        $crate::num::Complex::new(
            $re,
            <$crate::num::Complex<_> as $crate::num::Zero>::zero().im,
        )
    }};
    ($re:expr, $im: expr) => {{
        $crate::num::Complex::new($re, $im)
    }};
}

/// A more convenient way of writing `Poly::new(&[Complex::new(...)...])`
///
/// It takes ownership of its arguments.
///
/// It can take a list of `Scalar` or `Complex<Scalar>`. If left empty, it is
/// equivalent to `Poly::zero()`.
///
/// # Examples
///
/// Basic syntax
/// ```
/// use rust_poly::{poly, Poly};
/// use num::Zero;
/// use num::Complex;
///
/// let p1: Poly<f32> = poly![];
/// let p2 = poly![1.0f32, 2.0, 3.0];
/// let p3 = poly![Complex::from(1.0), Complex::from(2.0), Complex::from(3.0)];
///
/// assert_eq!(p1, Poly::zero());
/// assert_eq!(p2, p3);
/// ```
///
/// Similarly to `vec!`, you can initialize a large polynomial where all coefficients
/// are equal like so:
/// ```
/// # use rust_poly::{poly, Poly};
/// use num::Complex;
///
/// let p1 = poly![2.0; 16];
/// let p2 = poly![Complex::from(2.0); 16];
///
/// assert_eq!(p1, p2);
/// ```
///
/// You can also express complex numbers as a tuple of two scalars, mixing and matching
/// this syntax with the other syntax rules:
/// ```
/// use rust_poly::{poly, Poly};
/// use num::Complex;
///
/// let p1 = poly![(1.0, 2.0), (1.0, 2.0)];
/// let p2 = poly![(1.0, 2.0); 2];
/// let p3 = poly![Complex::new(1.0, 2.0); 2];
/// let p4 = poly![Complex::new(1.0, 2.0), Complex::new(1.0, 2.0)];
///
/// assert_eq!(p1, p2);
/// assert_eq!(p1, p3);
/// assert_eq!(p1, p4);
/// ```
#[macro_export]
macro_rules! poly {
    () => {{
        $crate::Poly::zero()
    }};
    (($re:expr, $im:expr); $n:expr) => {{
        $crate::Poly::from_complex_vec(vec![$crate::complex!($re, $im); $n])
    }};
    ($elem:expr; $n:expr) => {{
        $crate::Poly::from(vec![$elem; $n])
    }};
    ($(($re:expr, $im:expr)),+ $(,)?) => {{
        $crate::Poly::from(vec![$($crate::complex!($re, $im)),*])
    }};
    ($($elems:expr),+ $(,)?) => {{
        $crate::Poly::from(vec![$($elems),*])
    }};
}

mod scalar;
pub use scalar::{FloatScalar, Scalar, ScalarOps};

mod complex_util;
use complex_util::{c_neg, complex_sort_mut};
mod bessel;
mod casting_util;
mod linalg_util;
use casting_util::usize_to_u32;
mod poly;
// TODO: wildcard bad
pub use poly::*;
mod util;

pub type Poly32 = Poly<f32>;
pub type Poly64 = Poly<f64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn macro_complex() {
        assert_eq!(complex!(), Complex::<f64>::zero());
        assert_eq!(complex!(1.0, 2.0), Complex::<f64>::new(1.0, 2.0));
    }

    #[test]
    fn macro_poly() {
        assert_eq!(poly!(), Poly::<f64>::zero());
        assert_eq!(poly!(1.0), Poly::<f64>::one());
        assert_eq!(
            poly!(1.0, 2.0, 3.0),
            Poly::<f64>::new(&[
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
            ])
        );
        assert_eq!(poly!((1.0, 0.0)), Poly::<f64>::one());
        assert_eq!(
            poly!((1.0, 1.0), (2.0, 2.0), (3.0, 3.0)),
            Poly::<f64>::new(&[
                Complex::new(1.0, 1.0),
                Complex::new(2.0, 2.0),
                Complex::new(3.0, 3.0)
            ])
        );
        assert_eq!(
            poly!(2.0; 3),
            Poly::<f64>::new(&[
                Complex::new(2.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(2.0, 0.0)
            ])
        );
        assert_eq!(
            poly!((1.0, -1.0); 3),
            Poly::<f64>::new(&[
                Complex::new(1.0, -1.0),
                Complex::new(1.0, -1.0),
                Complex::new(1.0, -1.0)
            ])
        );
    }

    #[test]
    fn poly_new() {
        // trivial, but here for completeness
        Poly::new(&[Complex::new(2.0, -2.0)]);
    }

    #[test]
    fn poly_from_complex_slice() {
        let p = Poly::from_complex_slice(&[Complex::new(1.0, 2.0), Complex::new(3.0, 4.0)]);
        let e = poly!((1.0, 2.0), (3.0, 4.0));
        assert_eq!(p, e);
    }

    // TODO: test the rest of the "boring" functions

    #[test]
    fn poly_line() {
        let p = Poly::<f64>::line(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(2.0, 0.0));
        let e = poly!(1.0, 2.0);
        assert_eq!(p, e);
    }

    #[test]
    fn poly_term() {
        let p = Poly64::term(complex!(2.0), 2);
        let e = poly!(0.0, 0.0, 2.0);
        assert_eq!(p, e);
    }

    #[test]
    fn poly_bessel() {
        assert_eq!(Poly64::bessel(0).unwrap(), poly![1.0]);
        assert_eq!(Poly64::bessel(1).unwrap(), poly![1.0, 1.0]);
        assert_eq!(Poly64::bessel(2).unwrap(), poly![1.0, 3.0, 3.0]);
        assert_eq!(Poly64::bessel(3).unwrap(), poly![1.0, 6.0, 15.0, 15.0]);
    }

    #[test]
    fn poly_reverse_bessel() {
        assert_eq!(Poly64::reverse_bessel(0).unwrap(), poly![1.0]);
        assert_eq!(Poly64::reverse_bessel(1).unwrap(), poly![1.0, 1.0]);
        assert_eq!(Poly64::reverse_bessel(2).unwrap(), poly![3.0, 3.0, 1.0]);
        assert_eq!(
            Poly64::reverse_bessel(3).unwrap(),
            poly![15.0, 15.0, 6.0, 1.0]
        );
    }
}
