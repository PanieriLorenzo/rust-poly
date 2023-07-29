// Implementation of traits related to numeric operations, operators and number theory

use num_complex::Complex;
use num_traits::{One, Zero};
use std::ops::{Add, Mul};

extern crate nalgebra as na;

use crate::{Poly, Scalar};

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
        Self(self.0 + rhs.0.clone())
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
        Self(self.0 * rhs.0.clone())
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

impl<T: Scalar> std::iter::Sum for Poly<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}
