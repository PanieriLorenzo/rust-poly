use crate::{util::complex::c_neg, Scalar, ScalarOps};

use super::Poly;
use num::{traits::float::FloatCore, Complex, One, Zero};

impl<T> Poly<T> {
    /// The length of the polynomial without checking pre-conditions
    pub(crate) fn len_raw(&self) -> usize {
        self.0.len()
    }

    /// The degree of the polynomial without checking pre-conditions
    pub(crate) fn degree_raw(&self) -> i32 {
        self.len_raw() as i32 - 1
    }
}

impl<T: Scalar> Poly<T> {
    pub(crate) fn is_normalized(&self) -> bool {
        let n = self.len_raw();
        if n == 0 {
            return true;
        }
        !self.0.index(n - 1).is_zero()
    }

    pub(crate) fn normalize(self) -> Self {
        if self.is_normalized() {
            return self;
        }
        let mut end = self.len_raw();
        loop {
            if end == 0 {
                return Self::zero();
            }
            if !self.0.as_slice()[end - 1].is_zero() {
                break;
            }
            end -= 1;
        }
        let ret = Self(na::DVector::from_column_slice(&self.0.as_slice()[0..end]));

        // post-condition: polynomial is now normalized
        debug_assert!(ret.is_normalized());
        ret
    }

    pub(crate) fn companion(&self) -> na::DMatrix<Complex<T>> {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        // pre-condition: poly has degree 1 or more
        assert!(
            self.len_raw() >= 2,
            "polynomials of degree 0 or less do not have a companion matrix"
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
}

impl<T: ScalarOps> Poly<T> {
    pub(crate) fn checked_div_impl(self, rhs: &Self) -> Option<Self> {
        Some(self.div_rem(rhs)?.0)
    }

    pub(crate) fn checked_rem_impl(self, rhs: &Self) -> Option<Self> {
        Some(self.div_rem(rhs)?.1)
    }
}

impl<T: ScalarOps + FloatCore> Poly<T> {
    // Check that the polynomial does not contain `NaN` or infinite values.
    pub(crate) fn is_well_formed(&self) -> bool {
        self.0.iter().all(|x| !x.is_nan() && x.is_finite())
    }
}
