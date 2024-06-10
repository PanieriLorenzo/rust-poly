#![allow(clippy::op_ref)]

// Implementation of traits related to numeric operations, operators and number theory

use itertools::Itertools;
use num::{traits::CheckedRem, CheckedDiv, Complex, One, Zero};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

extern crate nalgebra as na;

use crate::{Poly, Scalar, ScalarOps, __util::linalg::convolve_1d};

impl<T: Scalar> Zero for Poly<T> {
    fn zero() -> Self {
        Poly::from_real_slice(&[T::zero()])
    }

    fn is_zero(&self) -> bool {
        debug_assert!(self.is_normalized());
        if self.len_raw() != 1 {
            return false;
        }
        self.0[0].is_zero()
    }
}

impl<T: ScalarOps> Poly<T> {
    /// Calculate the quotient and remainder using long division. More efficient than
    /// calculating them separately.
    ///
    /// # Panics
    /// Panics if a division by zero is attempted
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{Poly, poly};
    /// use num::{Complex, One};
    ///
    /// let c1 = poly![1.0, 2.0, 3.0];
    /// let c2 = poly![3.0, 2.0, 1.0];
    /// let expected1 = (poly![3.0], poly![-8.0, -4.0]);
    /// assert_eq!(c1.clone().div_rem(&c2).unwrap(), expected1);
    /// ```
    #[allow(clippy::cast_sign_loss)]
    #[allow(clippy::cast_possible_wrap)]
    #[must_use]
    pub fn div_rem(self, other: &Self) -> Option<(Self, Self)> {
        debug_assert!(self.is_normalized());
        debug_assert!(other.is_normalized());

        if other.is_zero() {
            // bail!("Attempted to divide a polynomial by zero");
            return None;
        }

        let expected_degree = self.degree_raw() - other.degree_raw();

        let den_c = other.as_slice().last().unwrap();
        let den_k = other.degree_raw();
        let mut this = self;
        let mut res = Self::zero();
        'for_else: {
            for _ in 0..u32::MAX {
                if this.degree_raw() < other.degree_raw() {
                    break 'for_else;
                }
                let num_c = this.as_slice().last().unwrap();
                let num_k = this.degree_raw();
                let c = num_c / den_c;
                let k = num_k - den_k;
                let new_term = Self::term(c, k as u32);
                this = this.clone() - new_term.clone() * other;
                res = res + new_term;
            }
            // else: did not converge
            return None;
        }
        let rem = this;

        // sanity check: result has the expected degree
        debug_assert_eq!(res.degree_raw(), expected_degree);
        Some((res, rem))
    }
}

impl<T: Scalar> One for Poly<T> {
    fn one() -> Self {
        Self(na::DVector::from_vec(vec![Complex::<T>::one()]))
    }
}

impl<T: Scalar> Add<Self> for Poly<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        let (mut longest, shortest) = if self.len_raw() >= rhs.len_raw() {
            (self.0, rhs.0)
        } else {
            (rhs.0, self.0)
        };
        longest
            .as_mut_slice()
            .iter_mut()
            .zip_longest(shortest.iter())
            .for_each(|p| {
                if let itertools::EitherOrBoth::Both(l, r) = p {
                    *l = *l + r;
                }
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Add<&Self> for Poly<T> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        self + rhs.clone()
    }
}

impl<T: Scalar> Add<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, rhs: Poly<T>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<T: Scalar> Add<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, rhs: &Poly<T>) -> Self::Output {
        self.clone() + rhs.clone()
    }
}

impl<T: Scalar> Mul<Self> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        if self.is_one() {
            return rhs;
        }
        if rhs.is_one() {
            return self;
        }

        let ret = convolve_1d(&self.0, &rhs.0);
        Self(ret).normalize()
    }
}

impl<T: Scalar> Mul<&Self> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        self * rhs.clone()
    }
}

impl<T: Scalar> Mul<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: Poly<T>) -> Self::Output {
        self.clone() * rhs
    }
}

impl<T: Scalar> Mul<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: &Poly<T>) -> Self::Output {
        self.clone() * rhs.clone()
    }
}

impl<T: Scalar> Mul<Complex<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: Complex<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: Scalar> Mul<&Complex<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        let mut lhs = self;
        lhs.0.apply(|c| *c = *c * rhs);
        lhs.normalize()
    }
}

impl<T: Scalar> Mul<Complex<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: Complex<T>) -> Self::Output {
        self.clone() * rhs
    }
}

impl<T: Scalar> Mul<&Complex<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        self.clone() * rhs
    }
}

impl<T: Scalar> Sub<Self> for Poly<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        let (mut longest, shortest) = if self.len_raw() >= rhs.len_raw() {
            (self.0, rhs.0)
        } else {
            (rhs.0, self.0)
        };
        longest
            .as_mut_slice()
            .iter_mut()
            .zip_longest(shortest.iter())
            .for_each(|p| {
                if let itertools::EitherOrBoth::Both(l, r) = p {
                    *l = *l - r;
                }
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Sub<&Self> for Poly<T> {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        self - rhs.clone()
    }
}

impl<T: Scalar> Sub<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, rhs: Poly<T>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<T: Scalar> Sub<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, rhs: &Poly<T>) -> Self::Output {
        self.clone() - rhs.clone()
    }
}

impl<T: ScalarOps> CheckedDiv for Poly<T> {
    fn checked_div(&self, rhs: &Self) -> Option<Self> {
        self.clone().checked_div_impl(rhs)
    }
}

impl<T: ScalarOps> CheckedRem for Poly<T> {
    fn checked_rem(&self, rhs: &Self) -> Option<Self> {
        self.clone().checked_rem_impl(rhs)
    }
}

impl<T: ScalarOps> Div<Self> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self / &rhs
    }
}

impl<T: ScalarOps> Div<&Self> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: &Self) -> Self::Output {
        self.checked_div_impl(rhs).expect("Division by zero")
    }
}

impl<T: ScalarOps> Div<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: Poly<T>) -> Self::Output {
        self.clone() / &rhs
    }
}

impl<T: ScalarOps> Div<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: &Poly<T>) -> Self::Output {
        self.clone() / rhs
    }
}

impl<T: Scalar> Div<Complex<T>> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: Complex<T>) -> Self::Output {
        self / &rhs
    }
}

impl<T: Scalar> Div<&Complex<T>> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: &Complex<T>) -> Self::Output {
        let mut lhs = self;
        lhs.0.apply(|c| *c = *c / rhs);
        lhs.normalize()
    }
}

impl<T: Scalar> Div<Complex<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: Complex<T>) -> Self::Output {
        self.clone() / rhs
    }
}

impl<T: Scalar> Div<&Complex<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: &Complex<T>) -> Self::Output {
        self.clone() / rhs
    }
}

impl<T: ScalarOps> Rem<Self> for Poly<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        self % &rhs
    }
}

impl<T: ScalarOps> Rem<&Self> for Poly<T> {
    type Output = Self;

    fn rem(self, rhs: &Self) -> Self::Output {
        self.checked_rem_impl(rhs).expect("Division by zero")
    }
}

impl<T: ScalarOps> Rem<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn rem(self, rhs: Poly<T>) -> Self::Output {
        self.clone() % &rhs
    }
}

impl<T: ScalarOps> Rem<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn rem(self, rhs: &Poly<T>) -> Self::Output {
        self.clone() % rhs
    }
}

impl<T: Scalar> Neg for Poly<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::zero() - self
    }
}

impl<T: Scalar> Neg for &Poly<T> {
    type Output = Poly<T>;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl<T: Scalar> std::iter::Sum for Poly<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x).normalize()
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

#[cfg(test)]
mod test {
    #[test]
    fn div() {
        let dividend = poly![-4.0, 0.0, -2.0, 1.0];
        let divisor = poly![-3.0, 1.0];
        let (q, r) = dividend.div_rem(&divisor).unwrap();
        assert_eq!(q, poly![3.0, 1.0, 1.0]);
        assert_eq!(r, poly![5.0]);
    }
}
