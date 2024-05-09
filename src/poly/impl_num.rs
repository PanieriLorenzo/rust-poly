// Implementation of traits related to numeric operations, operators and number theory

use itertools::Itertools;
use num::{traits::CheckedRem, CheckedDiv, Complex, One, Zero};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

extern crate nalgebra as na;

use crate::{util::linalg::convolve_1d, Poly, Scalar, ScalarOps};

impl<T: ScalarOps> Poly<T> {
    /// Calculate the quotient and remainder uwing long division. More efficient than
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
    pub fn div_rem(self, rhs: &Self) -> Option<(Self, Self)> {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        // pre-condition: don't divide by zero
        if rhs.is_zero() {
            // bail!("Attempted to divide a polynomial by zero");
            return None;
        }

        let lhs_len = self.len_raw();
        let rhs_len = self.len_raw();
        if lhs_len < rhs_len {
            return Some((Self::zero(), self));
        }
        if rhs_len == 1 {
            return Some((
                // TODO: should use checked operations
                Self(self.0 / rhs.0[rhs.len_raw() - 1].clone()),
                Self::zero(),
            ));
        }
        let len_delta = lhs_len - rhs_len;
        let scale = rhs.0[rhs.len_raw() - 1].clone();
        let rhs: na::DVector<_> = rhs
            .0
            .view_range(0..rhs.len_raw() - 1, 0..1)
            // HACK: this shouldn't be necessary, but nalgebra turns DVector into
            //       DMatrix when making a view, and needs to be politely reminded
            //       that this is a column vector.
            .column(0)
            .into();
        // TODO: useless clone of scale, it should be borrowed, but dvector does
        //       not implement Div<&_>
        // TODO: should use checked operations
        let rhs: na::DVector<_> = rhs / scale.clone();
        let mut lhs: na::DVector<_> = self.0.clone();
        let mut i = len_delta as isize;
        let mut j = (lhs_len - 1) as isize;
        while i >= 0 {
            lhs.view_range_mut(i as usize..j as usize, 0..1)
                .iter_mut()
                .zip((rhs.clone() * self.0[j as usize].clone()).iter())
                .for_each(|p| *p.0 -= p.1);
            i -= 1;
            j -= 1;
        }
        Some((
            Self(
                (lhs.view_range(j as usize + 1..lhs.len(), 0..1) / scale)
                    .column(0)
                    .into(),
            )
            .normalize(),
            Self(lhs.view_range(..(j + 1) as usize, 0..1).column(0).into()).normalize(),
        ))
    }
}

impl<T: Scalar> One for Poly<T> {
    fn one() -> Self {
        Self(na::DVector::from_vec(vec![Complex::<T>::one()]))
    }
}

impl<T: Scalar> Zero for Poly<T> {
    fn zero() -> Self {
        Self(na::DVector::from_vec(vec![]))
    }

    fn is_zero(&self) -> bool {
        self.is_empty()
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
                    *l = l.clone() + r;
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
        lhs.0.apply(|c| *c = c.clone() * rhs);
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
                    *l = l.clone() - r;
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
        lhs.0.apply(|c| *c = c.clone() / rhs);
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
