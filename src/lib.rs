#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]



extern crate nalgebra as na;
pub use num_complex;

use num_complex::Complex;
use num_traits::{One, Zero};

macro_rules! complex {
    ($re:expr, $im:expr) => {
        $crate::num_complex::Complex::new($re, $im)
    };
}

mod scalar;
pub use scalar::Scalar;

// mod roots;
// pub use roots::Roots;

mod complex_util;
use complex_util::{c_neg, complex_sort_mut};
mod impl_num;
mod num_util;

mod linalg_util;


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<T: Scalar>(na::DVector<Complex<T>>);

impl<T: Scalar> Poly<T> {
    pub fn new(coeffs: &[Complex<T>]) -> Self {
        Self(na::DVector::from_row_slice(coeffs))
    }

    #[must_use] pub fn from_roots(roots: na::DVector<Complex<T>>) -> Self {
        if roots.is_empty() {
            return Self::one();
        }

        let mut roots: na::DVector<Complex<T>> = roots;
        complex_sort_mut(&mut roots);

        roots
            .as_slice()
            .iter()
            .map(|e| Self::line(c_neg(e.clone()), Complex::<T>::one()))
            .fold(Self::one(), |acc, x| acc * x)
    }

    pub fn line(offset: Complex<T>, scale: Complex<T>) -> Self {
        if scale.is_zero() {
            return Self::new(&[offset]);
        }
        Self::new(&[offset, scale])
    }

    fn len_raw(&self) -> usize {
        self.0.len()
    }

    #[must_use] pub fn len(&self) -> usize {
        self.normalize().len_raw()
    }

    fn is_normalized(&self) -> bool {
        let n = self.len_raw();
        !self.0.index(n - 1).is_zero()
    }

    fn normalize(&self) -> Self {
        if self.len_raw() == 0 {
            return self.clone();
        }
        // while self.0.iter().last().unwrap().is_zero() {
        //     self.0.remove_row(self.len_raw() - 1);
        // }
        let mut end = self.len_raw();
        loop {
            if !self.0.iter().last().unwrap().is_zero() {
                break;
            }
            end -= 1;
        }
        Self(na::DVector::from_column_slice(&self.0.as_slice()[0..end]))
    }

    #[must_use] pub fn pow(&self, pow: u32) -> Self {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        if pow == 0 {
            return Self::one();
        }

        if pow == 1 {
            return self.clone();
        }

        // TODO: divide and conquer with powers of 2
        let mut res = self.clone();
        for _ in 2..=pow {
            res = res * self;
        }
        res
    }

    /// ```
    /// use rust_poly::Poly;
    /// use rust_poly::num_complex::Complex;
    ///
    /// let p = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0), Complex::new(0.0, -1.5)]);
    /// ```
    #[must_use] pub fn companion(&self) -> na::DMatrix<Complex<T>> {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        // pre-condition: poly has degree 1 or more
        assert!(
            self.len_raw() >= 2,
            "Poly must have maximum degree of at least 1"
        );

        if self.len_raw() == 2 {
            return na::DMatrix::from_row_slice(
                1,
                1,
                &[c_neg(self.0[0].clone()) / self.0[1].clone()],
            );
        }

        let n = self.len_raw() - 1;
        let mut mat: na::DMatrix<Complex<T>> = na::DMatrix::<Complex<T>>::zeros(n, n);

        // fill sub-diagonal with 1
        mat.view_mut((1, 0), (n - 1, n - 1))
            .fill_diagonal(Complex::<T>::one());

        // fill the rightmost column with the coefficients of the associated
        // monic polynomial
        let monic = self
            .0
            .view((0, 0), (n, 1))
            .map(|x| c_neg(x) / self.0[n].clone());
        for i in 0..n {
            mat.column_mut(n - 1)[i] = monic[i].clone();
        }
        mat
    }

    /// ```
    /// use rust_poly::Poly;
    /// use rust_poly::num_complex::Complex;
    ///
    /// let p = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)]);
    /// dbg!(p.roots());
    /// assert!(false);
    /// ```
    #[must_use] pub fn roots(&self) -> Option<na::DVector<Complex<T>>> {
        // invariant: polynomial is normalized
        debug_assert!(self.is_normalized());

        if self.len_raw() < 2 {
            return Some(na::dvector![]);
        }

        if self.len_raw() == 2 {
            return Some(na::dvector![c_neg(self.0[0].clone()) / self.0[1].clone()]);
        }

        // rotated companion matrix reduces error
        let mut comp = self.companion();
        let n = comp.shape().0;
        for i in 0..n / 2 {
            comp.swap_rows(i, n - i - 1);
            comp.swap_columns(i, n - i - 1);
        }

        let mut r: na::DVector<Complex<T>> = comp.eigenvalues()?;
        complex_sort_mut(&mut r);
        Some(r)
    }

    /// Compose two polynomials, returning a new polynomial.
    ///
    /// Substitute the given polynomial `x` into `self` and expand the
    /// result into a new polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_poly::Poly;
    /// use num_complex::Complex;
    /// use num_traits::identities::One;
    ///
    /// let f = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)]);
    /// let g = Poly::one();
    ///
    /// assert_eq!(f.compose(g), f);
    #[must_use] pub fn compose(&self, x: Self) -> Self {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(x.is_normalized());

        // TODO begin: are these checks actually making things faster?
        if self.is_zero() || x.is_zero() {
            return Self::zero();
        }

        if self.is_one() {
            return x;
        }

        if x.is_one() {
            return self.clone();
        }
        // end

        // TODO: prove that composing two normalized polynomials always results
        //       in a normalized polynomial or else disprove and call .normalize()
        (0..self.len_raw())
            .map(|i| Self::new(&[self.0[i].clone()]) * x.pow(i as u32))
            .sum()
    }
}
