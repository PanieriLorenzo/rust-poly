use na::Complex;
use num::{One, Zero};

use crate::{Poly, Scalar, __util::complex::c_neg, scalar::SafeConstants};

use super::indexing::Get;

impl<T: Scalar> Poly<T> {
    /// Applies a closure to each coefficient in-place
    pub(crate) fn apply(&mut self, f: impl FnMut(&mut Complex<T>)) {
        self.0.apply(f);
    }

    /// The length of the polynomial without checking pre-conditions
    pub(crate) fn len_raw(&self) -> usize {
        self.0.len()
    }

    /// Scale a polynomial in-place
    pub(crate) fn scale(&mut self, factor: Complex<T>) {
        self.apply(|z| *z = *z * factor);
    }

    /// Moving version of `scale`
    pub(crate) fn scaled(mut self, factor: Complex<T>) -> Self {
        self.scale(factor);
        self
    }

    /// The degree of the polynomial without checking pre-conditions
    #[inline]
    pub(crate) fn degree_raw(&self) -> usize {
        self.len_raw() - 1
    }

    pub(crate) fn is_normalized(&self) -> bool {
        let n = self.len_raw();
        if n == 0 {
            // zero-polynomials not allowed
            // TODO: ZeroPoly type, requires Poly trait
            return false;
        }
        // a constant is always normalized, as it may be just a constant zero
        if n == 1 {
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

    // TODO: move this somewhere else, it's not a "base" algorithm
    pub(crate) fn companion(&self) -> na::DMatrix<Complex<T>> {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        // pre-condition: poly has degree 1 or more
        assert!(
            self.degree_raw() >= 1,
            "polynomials of degree 0 or less do not have a companion matrix"
        );

        if self.len_raw() == 2 {
            return na::DMatrix::from_row_slice(1, 1, &[c_neg(self.0[0]) / self.0[1]]);
        }

        let n = self.len_raw() - 1;
        let mut mat: na::DMatrix<Complex<T>> = na::DMatrix::zeros(n, n);

        // fill sub-diagonal with 1
        mat.view_mut((1, 0), (n - 1, n - 1))
            .fill_diagonal(Complex::<T>::one());

        // fill the rightmost column with the coefficients of the associated
        // monic polynomial
        let mut monic = self.clone();
        monic.make_monic();
        for i in 0..n {
            mat.column_mut(n - 1)[i] = c_neg(monic[i]);
        }
        mat
    }

    /// The last coefficient
    pub(crate) fn last(&self) -> Complex<T> {
        self.0[self.len_raw() - 1]
    }

    /// Make the polynomial monic in-place.
    ///
    /// Monic polynomials are scaled such that the last coefficient is 1, and
    /// the roots are preserved
    pub(crate) fn make_monic(&mut self) {
        debug_assert!(self.is_normalized());
        let last_coeff = self.last();
        if last_coeff.is_one() {
            // already monic
            return;
        }
        self.apply(|x| *x = *x / last_coeff);
    }

    /// Make sure trailing almost-zero coefficients are removed
    pub(crate) fn trim(&mut self) {
        let last_coeff = self.last();
        if last_coeff.is_small() {
            let mut res = self.get(0..self.len_raw() - 1).expect("infallible");
            res.trim();
            *self = res;
        }
        *self = self.clone().normalize();
    }

    // TODO: deflate_downward would be a better name
    /// Factor out one root of the polynomial, by scaling coefficients from
    /// the one with the highest degree down, then discarding the smallest
    /// coefficient.
    pub(crate) fn deflate_forward(mut self, r: Complex<T>) -> Self {
        // TODO: it is possible to use FFT for forward deflation
        let mut z0 = Complex::<T>::zero();
        // FIXME: I think this does exactly one wasted iteration at the end
        for j in 0..self.len_raw() {
            z0 = z0 * r + *self.coeff_descending(j);
            *self.coeff_descending_mut(j) = z0;
        }
        self.shift_down(1).normalize()
    }

    // TODO: deflate_upward would be a better name
    // TODO: single letter names
    /// Factor out one root of the polynomial, by scaling coefficients from
    /// the one with the lowest degree upwards, then discarding the largest
    /// coefficient.
    #[allow(clippy::many_single_char_names)]
    pub(crate) fn deflate_backward(mut self, r: Complex<T>) -> Self {
        let n = self.degree_raw();
        if n == 0 {
            return self;
        }
        let mut z0 = Complex::zero();
        if r != z0 {
            let mut i = n - 1;
            let mut t = *self.coeff_descending(n.try_into().expect("overflow"));
            let mut s;
            loop {
                s = t;
                t = *self.coeff_descending(i.try_into().expect("overflow"));
                z0 = (z0 - s) / r;
                *self.coeff_descending_mut(i.try_into().expect("overflow")) = z0;
                i -= 1;
                if i == 0 {
                    break;
                }
            }
        }
        self.shift_down(1).normalize()
    }

    /// Synthetic division that reduces numeric error by fusing the results
    /// of forward deflation and backward deflation
    pub(crate) fn deflate_composite(mut self, r: Complex<T>) -> Self {
        // TODO: should take a mutable reference instead
        let n = self.degree_raw();
        let fwd = self.clone().deflate_forward(r);
        let bwd = self.clone().deflate_backward(r);
        // TODO: in order to drop the Bounded trait bound, this should be
        //       done without explicit reference to max value
        let mut ra = T::max_value();
        let mut ua;
        let mut k = 0;
        for i in 0..n {
            ua = fwd.coeff_descending(i.try_into().expect("overflow")).norm()
                + bwd.coeff_descending(i.try_into().expect("overflow")).norm();
            if !ua.is_zero() {
                ua = (fwd.coeff_descending(i.try_into().expect("overflow"))
                    - bwd.coeff_descending(i.try_into().expect("overflow")))
                .norm()
                    / ua;
                if ua < ra {
                    ra = ua;
                    k = i;
                }
            }
        }
        let mut i = k.saturating_sub(1);
        loop {
            *self.coeff_descending_mut(i.try_into().expect("overflow")) =
                *fwd.coeff_descending(i.try_into().expect("overflow"));
            if i == 0 {
                break;
            }
            i -= 1;
        }
        *self.coeff_descending_mut(k.try_into().expect("overflow")) = (fwd
            .coeff_descending(k.try_into().expect("overflow"))
            + bwd.coeff_descending(k.try_into().expect("overflow")))
        .scale(T::from_u8(2).expect("should fit").recip());
        for i in (k + 1)..n {
            *self.coeff_descending_mut(i.try_into().expect("overflow")) =
                *bwd.coeff_descending(i.try_into().expect("overflow"));
        }

        self.shift_down(1).normalize()
    }
}

// impl<T: Scalar + PartialOrd> Poly<T> {
//     /// Check if the polynomial is actually a monomial
//     pub(crate) fn is_monomial(&self, tol: T) -> bool {
//         debug_assert!(self.is_normalized());
//         let degree = self.degree_raw();
//         if degree < 0 {
//             return false;
//         }
//         let degree: usize = degree.try_into().expect("overflow");
//         self.iter().take(degree).all(|z| z.norm_sqr() < tol)
//     }
// }

#[cfg(test)]
mod test {
    use na::{DVector, Matrix3};
    use num::{complex::Complex64, Zero};

    use crate::Poly;

    /// This was a bug
    #[test]
    fn normalize0() {
        let p = Poly(DVector::from_column_slice(&[Complex64::zero()]));
        assert_eq!(p.normalize().0.as_slice(), &[Complex64::zero()]);
    }

    /// This was a bug
    #[test]
    fn is_normalized0() {
        let p = Poly(DVector::from_column_slice(&[Complex64::zero()]));
        assert!(p.is_normalized());
    }

    #[test]
    fn companion() {
        let p = poly![1.0, 2.0, 3.0, 4.0];
        let c = p.companion();
        let c_expected =
            Matrix3::new(0., 0., -0.25, 1., 0., -0.5, 0., 1., -0.75).cast::<Complex64>();
        assert_eq!(c, c_expected);
    }

    #[test]
    fn companion_tiny() {
        let p = poly![1.0, 2.0];
        let c = p.companion();
        assert_eq!(c[0], complex!(-0.5));
    }

    #[test]
    fn monic() {
        let mut p = poly![1.0, 3.0, 2.0];
        p.make_monic();
        assert_eq!(p, poly![0.5, 3.0 / 2.0, 1.0]);
    }
}
