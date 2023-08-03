// Implementation of traits related to numeric operations, operators and number theory

use itertools::Itertools;
use na::DVector;
use num_complex::Complex;
use num_traits::{One, Zero};
use std::{
    borrow::BorrowMut,
    ops::{Add, Div, Mul, Neg, Sub},
};

extern crate nalgebra as na;

use crate::{linalg_util::convolve_1d, Poly, Scalar};

impl<T: Scalar> One for Poly<T> {
    fn one() -> Self {
        Self(na::DVector::from_vec(vec![Complex::<T>::one()]))
    }
}

impl<T: Scalar> Zero for Poly<T> {
    fn zero() -> Self {
        Self(na::DVector::from_vec(vec![Complex::<T>::zero()]))
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
            .for_each(|p| match p {
                itertools::EitherOrBoth::Both(l, r) => *l += r,
                _ => (),
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Add for Poly<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<T: Scalar> Mul<&Self> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        let ret = convolve_1d(self.0, &rhs.0);
        Self(ret).normalize()
    }
}

impl<T: Scalar> Mul for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: Scalar> Mul<&Complex<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        Self(self.0.map(|e| e * rhs))
    }
}

impl<T: Scalar> Mul<Complex<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, rhs: Complex<T>) -> Self::Output {
        self.mul(&rhs)
    }
}

impl<T: Scalar> Sub<&Poly<T>> for Poly<T> {
    type Output = Self;

    fn sub(self, rhs: &Poly<T>) -> Self::Output {
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
            .for_each(|p| match p {
                itertools::EitherOrBoth::Both(l, r) => *l -= r,
                _ => (),
            });
        Self(longest).normalize()
    }
}

impl<T: Scalar> Sub<Poly<T>> for Poly<T> {
    type Output = Self;

    fn sub(self, rhs: Poly<T>) -> Self::Output {
        self - &rhs
    }
}

impl<T: Scalar> Neg for Poly<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<T: Scalar> std::iter::Sum for Poly<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x).normalize()
    }
}
