// Implementation of traits related to numeric operations, operators and number theory

use itertools::Itertools;
use num::{traits::CheckedRem, CheckedDiv, Complex, One, Zero};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

extern crate nalgebra as na;

use crate::{linalg_util::convolve_1d, Poly, Scalar};

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

impl<T: Scalar> Add<&Self> for Poly<T> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        let (mut longest, shortest) = if self.len_raw() >= rhs.len_raw() {
            (self.0.clone(), &rhs.0)
        } else {
            (rhs.0.clone(), &self.0)
        };
        longest
            .as_mut_slice()
            .iter_mut()
            .zip_longest(shortest.iter())
            .for_each(|p| {
                if let itertools::EitherOrBoth::Both(l, r) = p {
                    *l += r;
                }
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Add for Poly<T> {
    type Output = Self;

    /// Add two polynomials
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    ///
    /// let p1 = poly![1.0, 2.0, 3.0];
    /// let p2 = poly![4.0, 5.0];
    /// assert_eq!(p1 + p2, poly![5.0, 7.0, 3.0]);
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<T: Scalar> Add<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, rhs: Poly<T>) -> Self::Output {
        rhs + self
    }
}

impl<T: Scalar> Add<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, rhs: &Poly<T>) -> Self::Output {
        self + rhs.clone()
    }
}

impl<T: Scalar> Mul<&Self> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        if self.is_one() {
            return rhs.clone();
        }
        if rhs.is_one() {
            return self;
        }

        let ret = convolve_1d(&self.0, &rhs.0);
        Self(ret).normalize()
    }
}

impl<T: Scalar> Mul for Poly<T> {
    type Output = Self;

    /// Multiply two polynomials.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// let p1 = poly![1.0, 2.0, 3.0];
    /// let p2 = poly![3.0, 2.0, 1.0];
    /// assert_eq!(p1 * p2, poly![3.0, 8.0, 14.0, 8.0, 3.0]);
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: Scalar> Mul<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: Poly<T>) -> Self::Output {
        rhs * self
    }
}

impl<T: Scalar> Mul<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, rhs: &Poly<T>) -> Self::Output {
        self * rhs.clone()
    }
}

impl<T: Scalar> Mul<&Complex<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        Self(self.0.map(|e| e * rhs)).normalize()
    }
}

impl<T: Scalar> Mul<Complex<T>> for Poly<T> {
    type Output = Self;

    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p * Complex::from(2.0), poly![2.0, 4.0, 6.0]);
    /// ```
    fn mul(self, rhs: Complex<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: Scalar> Sub<&Self> for Poly<T> {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        let (mut longest, shortest) = if self.len_raw() >= rhs.len_raw() {
            (self.0.clone(), &rhs.0)
        } else {
            (rhs.0.clone(), &self.0)
        };
        longest
            .as_mut_slice()
            .iter_mut()
            .zip_longest(shortest.iter())
            .for_each(|p| {
                if let itertools::EitherOrBoth::Both(l, r) = p {
                    *l -= r;
                }
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Sub<Self> for Poly<T> {
    type Output = Self;

    /// Subtract two polynomials.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// let c1 = poly![1.0, 2.0, 3.0];
    /// let c2 = poly![3.0, 2.0, 1.0];
    /// assert_eq!(c1.clone() - c2.clone(), poly![-2.0, 0.0, 2.0]);
    /// assert_eq!(c2 - c1, poly![2.0, 0.0, -2.0]);
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: Scalar> Sub<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, rhs: Poly<T>) -> Self::Output {
        rhs - self
    }
}

impl<T: Scalar> Sub<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, rhs: &Poly<T>) -> Self::Output {
        self - rhs.clone()
    }
}

impl<T: Scalar> CheckedDiv for Poly<T> {
    fn checked_div(&self, rhs: &Self) -> Option<Self> {
        self.clone().checked_div_impl(rhs)
    }
}

impl<T: Scalar> CheckedRem for Poly<T> {
    fn checked_rem(&self, rhs: &Self) -> Option<Self> {
        self.clone().checked_rem_impl(rhs)
    }
}

impl<T: Scalar> Div<&Self> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: &Self) -> Self::Output {
        self.checked_div(rhs).unwrap()
    }
}

impl<T: Scalar> Div<Self> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self / &rhs
    }
}

impl<T: Scalar> Div<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: Poly<T>) -> Self::Output {
        rhs / self
    }
}

impl<T: Scalar> Div<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn div(self, rhs: &Poly<T>) -> Self::Output {
        self / rhs.clone()
    }
}

impl<T: Scalar> Div<&Complex<T>> for Poly<T> {
    type Output = Self;

    fn div(self, rhs: &Complex<T>) -> Self::Output {
        Self(self.0.map(|e| e / rhs)).normalize()
    }
}

impl<T: Scalar> Div<Complex<T>> for Poly<T> {
    type Output = Self;

    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// let p = poly![2.0, 4.0, 6.0];
    /// assert_eq!(p / Complex::from(2.0), poly![1.0, 2.0, 3.0]);
    /// ```
    fn div(self, rhs: Complex<T>) -> Self::Output {
        self / &rhs
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
        self / rhs.clone()
    }
}

impl<T: Scalar> Rem<&Self> for Poly<T> {
    type Output = Self;

    fn rem(self, rhs: &Self) -> Self::Output {
        self.checked_rem(rhs).unwrap()
    }
}

impl<T: Scalar> Rem<Self> for Poly<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        self % &rhs
    }
}

impl<T: Scalar> Rem<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn rem(self, rhs: Poly<T>) -> Self::Output {
        rhs % self
    }
}

impl<T: Scalar> Rem<&Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn rem(self, rhs: &Poly<T>) -> Self::Output {
        self % rhs.clone()
    }
}

impl<T: Scalar> Neg for Poly<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
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
