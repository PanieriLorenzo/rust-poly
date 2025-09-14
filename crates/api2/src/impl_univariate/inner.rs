//! Private API that are composed to make public APIs

use std::mem;

use crate::{
    aliases::C,
    base::{BasePoly, UnivariateMarker},
    errors::{CAST_OVERFLOW, DEGREE_TOO_LARGE},
    scalar_traits::{BasicScalar, NonIntegerScalar, RealScalar},
    storage_traits::{BasicStorage, MutStorage, OwnedStorage},
};

use num_complex::Complex;
use num_traits::{FromPrimitive, One, Zero};

impl<S: BasicStorage> BasePoly<S, UnivariateMarker> {
    pub fn into_owned_poly_inner(&self) -> BasePoly<Vec<S::T>, UnivariateMarker> {
        let owned = self.storage.borrow().to_vec();
        // SAFETY: preservig original length should always be safe
        unsafe { BasePoly::from_parts(owned, self.vars) }
    }

    pub fn len_inner(&self) -> usize {
        self.storage.borrow().len()
    }

    /// Degree that saturates at 0, giving incorrect degree if polynomial is
    /// of degree -1.
    pub fn degree_raw_inner(&self) -> usize {
        self.len_inner().saturating_sub(1)
    }

    pub fn degree_inner(&self) -> i32 {
        if self.len_inner() == 0 {
            -1
        } else {
            self.degree_raw_inner().try_into().expect(DEGREE_TOO_LARGE)
        }
    }

    pub fn iter_inner(&self) -> impl Iterator<Item = &S::T> {
        self.storage.borrow().iter()
    }

    pub fn is_normalized_inner(&self) -> bool {
        !self.storage.borrow()[self.degree_raw_inner()].is_zero()
    }

    pub fn eval_complex_inner(
        &self,
        x: <S::T as NonIntegerScalar>::ToComplexScalar,
    ) -> <S::T as NonIntegerScalar>::ToComplexScalar
    where
        S::T: NonIntegerScalar,
    {
        // NOTE: partial duplicate of `eval_inner`, please keep in sync
        //
        // use Horner's method: https://en.wikipedia.org/wiki/Horner%27s_method
        // in theory, Estrin's method is more parallelizable, but benchmarking
        // against fast_polynomial crate shows no significant difference, this
        // is close to optimal in terms of performance. You may get some small
        // speedups by dividing large polynomials into 4 or 8 evaluations and
        // computing them in parallel using SIMD, Rayon or GPU.
        let mut eval = self.storage.borrow()[self.degree_raw_inner()].to_complex();
        let n = self.len_inner();
        for i in 1..n {
            // SAFETY: index n - i - 1 is always in bounds
            // PROOF:
            // 1. the safe bounds for the index are [0, n - 1]
            // 2. i is always in the range [1, n - 1]
            // 3. the range of the index n - i - 1 is given by:
            //    n - i - 1 => n - [1, n - 1] - 1
            //             <=> [n - 1 - 1, n - 1 - (n - 1)]
            //             <=> [n - 2, 0]
            //             <=> [0, n - 2]   reversing bounds is equivalent
            // 4. the range of the index [0, n - 2] is a subset of [0, n - 1],
            //    therefore the index is in bounds. QED
            let c = (unsafe { self.storage.borrow().get_unchecked(n - i - 1) }).to_complex();
            eval = eval * x.clone() + c;
        }
        eval
    }

    pub fn eval_inner(&self, x: S::T) -> S::T {
        // NOTE: partial duplicate of `eval_complex_inner`, please keep in sync
        //
        // use Horner's method: https://en.wikipedia.org/wiki/Horner%27s_method
        // in theory, Estrin's method is more parallelizable, but benchmarking
        // against fast_polynomial crate shows no significant difference, this
        // is close to optimal in terms of performance. You may get some small
        // speedups by dividing large polynomials into 4 or 8 evaluations and
        // computing them in parallel using SIMD, Rayon or GPU.
        let mut eval = self.storage.borrow()[self.degree_raw_inner()].clone();
        let n = self.len_inner();
        for i in 1..n {
            // SAFETY: index n - i - 1 is always in bounds
            // PROOF:
            // 1. the safe bounds for the index are [0, n - 1]
            // 2. i is always in the range [1, n - 1]
            // 3. the range of the index n - i - 1 is given by:
            //    n - i - 1 => n - [1, n - 1] - 1
            //             <=> [n - 1 - 1, n - 1 - (n - 1)]
            //             <=> [n - 2, 0]
            //             <=> [0, n - 2]   reversing bounds is equivalent
            // 4. the range of the index [0, n - 2] is a subset of [0, n - 1],
            //    therefore the index is in bounds. QED
            let c = (unsafe { self.storage.borrow().get_unchecked(n - i - 1) }).clone();
            eval = eval * x.clone() + c;
        }
        eval
    }

    pub fn eval_multiple_complex_inner(
        &self,
        points: &[<S::T as NonIntegerScalar>::ToComplexScalar],
        out: &mut [<S::T as NonIntegerScalar>::ToComplexScalar],
    ) where
        S::T: NonIntegerScalar,
    {
        // TODO: conditionally parallelize this loop
        for (y, x) in out.iter_mut().zip(points) {
            *y = self.eval_complex_inner(x.clone());
        }
    }
}

impl<S: MutStorage> BasePoly<S, UnivariateMarker> {
    pub fn iter_mut_inner(&mut self) -> impl Iterator<Item = &mut S::T> {
        self.storage.borrow_mut().iter_mut()
    }

    pub fn scale_inner(&mut self, factor: &S::T) {
        for x in self.iter_mut_inner() {
            *x = x.clone() * factor.clone();
        }
    }

    pub fn make_monic_inner(&mut self) {
        let last_coeff = self.storage.borrow()[self.degree_raw_inner()].clone();
        if last_coeff.is_one() {
            // already monic
            return;
        }
        for x in self.iter_mut_inner() {
            *x = x.clone() / last_coeff.clone();
        }
    }
}

impl<S: OwnedStorage> BasePoly<S, UnivariateMarker> {
    pub fn normalize_inner(&mut self) {
        // trim trailing zero coefficients
        let mut end = self.len_inner();
        loop {
            if end == 0 {
                return;
            }
            if !self.storage.borrow()[end - 1].is_zero() {
                break;
            }
            end -= 1;
        }
        self.storage.truncate(end);
    }

    /// Remove `n` leading coefficients, reducing the degree by `n`.
    pub fn shift_down_inner(&mut self, n: usize) {
        self.storage
            .truncate_start(self.len_inner().saturating_sub(n));
    }

    /// Remove trailing coefficients that are almost zero.
    pub fn trim_inner(&mut self) {
        let mut end = self.len_inner();
        loop {
            if end == 0 {
                return;
            }
            if !self.storage.borrow()[end - 1].taxicab_norm().is_tiny() {
                break;
            }
            end -= 1;
        }
        self.storage.truncate(end);
    }

    pub fn diff_inner(&mut self) {
        let len = self.len_inner();
        self.iter_mut_inner()
            .zip(0..len)
            .skip(1)
            .for_each(|(c, n)| *c = c.clone() * S::T::from_usize(n).expect(CAST_OVERFLOW));
        self.shift_down_inner(1);
    }

    /// Find all roots where `r = 0`, factoring them out of the polynomial.
    /// This is much faster than using root finders, and polynomials with zero
    /// roots are a quite common case, so they are worth checking for.
    pub fn zero_roots_inner(&mut self, roots_out: &mut Vec<C<S>>)
    where
        S::T: NonIntegerScalar,
    {
        for _ in 0..self.degree_raw_inner() {
            if self.eval_inner(S::T::zero()).taxicab_norm()
                <= <S::T as BasicScalar>::RealPartScalar::TINY
            {
                roots_out.push(C::<S>::zero());
                // deflating zero roots can be accomplished simply by shifting
                self.shift_down_inner(1);
            } else {
                break;
            }
        }
    }

    /// Check if the polynomial is of degree one, if it is, finding roots for
    /// it is trivial and saves us running an expensive root finder.
    pub fn linear_roots_inner(&mut self, roots_out: &mut Vec<C<S>>)
    where
        S::T: NonIntegerScalar,
    {
        if self.degree_raw_inner() != 1 {
            return;
        }
        let a = &self.storage.borrow()[1];
        let b = &self.storage.borrow()[0];
        roots_out.push((C::<S>::zero() - b.to_complex()) / a.to_complex());
        *self = Self::from_iter([S::T::one()]);
    }

    pub fn quadratic_roots_inner(&mut self, roots_out: &mut Vec<C<S>>)
    where
        S::T: NonIntegerScalar,
    {
        if self.degree_raw_inner() != 2 {
            return;
        }
        let a = self.storage.borrow()[2].to_complex();
        let b = self.storage.borrow()[1].to_complex();
        let c = self.storage.borrow()[0].to_complex();
        let (x1, x2) = C::<S>::quadratic_formula(a, b, c);
        roots_out.push(x1);
        roots_out.push(x2);
        *self = Self::from_iter([S::T::one()]);
    }
}

#[cfg(test)]
mod unit {
    use crate::aliases::RealPoly32;

    #[test]
    fn len_inner() {
        let poly = RealPoly32::from_iter([1.0, 2.0, 3.0]);
        assert_eq!(poly.len_inner(), 3);
    }

    #[test]
    fn degree_raw_inner() {
        let poly = RealPoly32::from_iter([1.0, 2.0, 3.0]);
        assert_eq!(poly.degree_raw_inner(), 2);
        let one = RealPoly32::from_iter(vec![1.0]);
        assert_eq!(one.degree_raw_inner(), 0);
        let zero = RealPoly32::from_iter(vec![]);
        assert_eq!(zero.degree_raw_inner(), 0);
    }

    #[test]
    fn degree_inner() {
        let poly = RealPoly32::from_iter([1.0, 2.0, 3.0]);
        assert_eq!(poly.degree_inner(), 2);
        let one = RealPoly32::from_iter(vec![1.0]);
        assert_eq!(one.degree_inner(), 0);
        let zero = RealPoly32::from_iter(vec![]);
        assert_eq!(zero.degree_inner(), -1);
    }

    #[test]
    fn normalize_inner() {
        let mut poly = unsafe { RealPoly32::from_parts(vec![1.0, 2.0, 0.0, 0.0], 1) };
        assert_eq!(poly.len_inner(), 4);
        poly.normalize_inner();
        assert_eq!(poly.len_inner(), 2);
    }

    #[test]
    fn scale_inner() {
        let mut poly = RealPoly32::from_iter(vec![1.0, 2.0, 3.0]);
        poly.scale_inner(&1.5);
        assert_eq!(poly, RealPoly32::from_iter(vec![1.5, 3.0, 4.5]));
    }
}
