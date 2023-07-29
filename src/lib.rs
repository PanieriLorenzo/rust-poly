#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]

use std::cmp::Ordering;

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
use complex_util::c_neg;
mod impl_num;
mod num_util;
use num_util::neg;
mod linalg_util;
use linalg_util::reverse_mut;

#[derive(Clone, Debug)]
pub struct Poly<T: Scalar>(na::DVector<Complex<T>>);

impl<T: Scalar> Poly<T> {
    pub fn new(coeffs: &[Complex<T>]) -> Self {
        Self(na::DVector::from_row_slice(coeffs))
    }

    fn len_raw(&self) -> usize {
        self.0.len()
    }

    pub fn len(&self) -> usize {
        self.normalize().len_raw()
    }

    fn is_normalized(&self) -> bool {
        let n = self.len_raw();
        !self.0.index(n - 1).is_zero()
    }

    fn normalize(&self) -> Self {
        let mut ret = self.clone();
        ret.normalize_mut();
        ret
    }

    fn normalize_mut(&mut self) {
        todo!()
    }

    /// ```
    /// use rust_poly::Poly;
    /// use rust_poly::num_complex::Complex;
    ///
    /// let p = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0), Complex::new(0.0, -1.5)]);
    /// ```
    pub fn companion(&self) -> na::DMatrix<Complex<T>> {
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
        dbg!(monic.clone());
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
    pub fn roots(&self) -> Option<na::DVector<Complex<T>>> {
        // invariant: poly is normalized
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
        dbg!(comp.clone());
        for i in 0..n / 2 {
            comp.swap_rows(i, n - i - 1);
            comp.swap_columns(i, n - i - 1);
        }
        dbg!(comp.clone());

        let mut r: na::DVector<Complex<T>> = comp.eigenvalues()?;
        r.as_mut_slice()
            .sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap_or(Ordering::Equal));
        Some(r)
    }
}
