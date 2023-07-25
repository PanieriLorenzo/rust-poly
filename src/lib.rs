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
use complex_util::c_neg;
mod impl_num;
mod num_util;
use num_util::neg;

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
    /// dbg!(p.companion());
    /// assert!(false);
    /// ```
    pub fn companion(&self) -> na::DMatrix<Complex<T>> {
        // pre-condition: poly is normalized
        assert!(self.is_normalized());

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
        let mut mat = na::DMatrix::<Complex<T>>::zeros(n, n);
        mat.view_mut((1, 0), (n - 1, n - 1))
            .fill_diagonal(Complex::<T>::one());
        let monic = self.0.view((0, 0), (n, 1)).map(|x| x / self.0[0].clone());
        for i in 0..(n - 1) {
            mat.column_mut(n - 1)[i] = mat.column(n - 1)[i].clone() - monic[i].clone();
        }
        mat
    }
}
