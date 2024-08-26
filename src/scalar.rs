//! Traits for the coefficients of polynomials

use num::{
    complex::ComplexFloat,
    traits::{bounds::UpperBounded, Float, MulAdd, MulAddAssign},
    Complex, FromPrimitive, Num, ToPrimitive,
};
use std::{
    fmt::Display,
    ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign},
};

/// A rational number, less restrictive than num::Float
pub trait Rational: std::ops::Div<Output = Self> + num::One + ToPrimitive + FromPrimitive {
    fn recip(self) -> Self {
        Self::one() / self
    }

    fn upper_bound() -> Self;
}

impl Rational for f32 {
    fn recip(self) -> Self {
        f32::recip(self)
    }

    fn upper_bound() -> Self {
        <f32 as Float>::max_value()
    }
}

impl Rational for f64 {
    fn recip(self) -> Self {
        f64::recip(self)
    }

    fn upper_bound() -> Self {
        <f64 as Float>::max_value()
    }
}

impl Rational for f128 {
    fn recip(self) -> Self {
        Float::recip(self)
    }

    fn upper_bound() -> Self {
        <f128 as Float>::max_value()
    }
}

// pub trait ComplexScalar: Scalar {
//     type ComponentScalar;
// }

// impl ComplexScalar for Complex64 {
//     type ComponentScalar = f64;
// }
// impl ComplexScalar for Complex32 {
//     type ComponentScalar = f32;
// }

/// The trait bounds necessary to provide the basic functionality of this crate.
#[allow(clippy::module_name_repetitions)]
pub trait RealScalar:
    Clone
    + PartialEq
    + std::fmt::Debug
    + Num
    + FromPrimitive
    + ToPrimitive
    + SafeConstants
    + std::ops::Neg<Output = Self>
    + 'static
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
    //+ MulAdd<Output = Self>
    //+ MulAddAssign
    + Display
    + PartialEq
    + PartialOrd
    //+ UpperBounded
    + Rational
{
}
impl<
        T: Clone
            + PartialEq
            + PartialOrd
            + std::fmt::Debug
            + Num
            + FromPrimitive
            + ToPrimitive
            + SafeConstants
            + std::ops::Neg<Output = Self>
            + AddAssign
            + SubAssign
            + MulAssign
            + DivAssign
            + RemAssign
            //+ MulAdd<Output = Self>
            //+ MulAddAssign
            + 'static
            + PartialEq
            + PartialOrd
            //+ UpperBounded
            + Rational
            + Display,
    > RealScalar for T
{
}

/// A number that has a smallest positive _safe_ value for denominators and a largest
/// positive _safe_ value for numerators.
///
/// These are arbitrary ranges, but they need to fulfill a few requirements:
/// - `big_safe` should be very large, enough that practically most numbers will be smaller
/// - `small_safe` should be very small, enough that practically most numbers will be bigger
/// - `big_safe / small_safe <= MAX` where MAX is the largerst finite number representible
///
/// Knowledge of these numbers make it possible to check for overflow in division,
/// if `numberator <= big_safe` and `denominator >= small_safe` the division will
/// always be finite.
///
/// If either numerator or denominator is outside of the range, then the division
/// is not guaranteed to be finite, but is not guaranteed to overflow either.
///
/// The trait additionally provides a minimum safe value for reciprocals, i.e.
/// if `numerator <= 1.0`.
pub trait SafeConstants {
    /// Smallest value that is safe for reciprocals
    fn tiny_safe() -> Self;

    /// Smallest value that is safe for use as a denominator
    fn small_safe() -> Self;

    /// Biggest value that is safe for use as a numerator
    fn large_safe() -> Self;

    /// If the number is smaller than [`SafeConstants::small_safe`], i.e. it is not safe for division
    fn is_small(&self) -> bool;

    /// If the number is smaller than [`SafeConstants::tiny_safe`], i.e. it is not safe for reciprocals
    fn is_tiny(&self) -> bool;

    /// If the number is larger than [`SafeConstants::large_safe`], i.e. it is not safe for division
    fn is_large(&self) -> bool;
}

macro_rules! safe_constants_float_impl {
    ($t:ty) => {
        impl SafeConstants for $t {
            fn tiny_safe() -> Self {
                Self::MIN_POSITIVE
            }

            fn small_safe() -> Self {
                Self::MIN_POSITIVE.sqrt() / Self::EPSILON
            }

            fn large_safe() -> Self {
                Self::MAX * Self::small_safe()
            }

            fn is_tiny(&self) -> bool {
                self.abs() < Self::tiny_safe()
            }

            fn is_small(&self) -> bool {
                self.abs() < Self::small_safe()
            }

            fn is_large(&self) -> bool {
                self.abs() > Self::large_safe()
            }
        }
    };
}

safe_constants_float_impl!(f32);
safe_constants_float_impl!(f64);

impl<T: SafeConstants> SafeConstants for Complex<T> {
    fn tiny_safe() -> Self {
        Self {
            re: T::tiny_safe(),
            im: T::tiny_safe(),
        }
    }

    fn small_safe() -> Self {
        Self {
            re: T::small_safe(),
            im: T::small_safe(),
        }
    }

    fn large_safe() -> Self {
        Self {
            re: T::large_safe(),
            im: T::large_safe(),
        }
    }

    fn is_small(&self) -> bool {
        self.re.is_small() && self.im.is_small()
    }

    fn is_tiny(&self) -> bool {
        self.re.is_tiny() && self.im.is_tiny()
    }

    fn is_large(&self) -> bool {
        self.re.is_large() || self.im.is_large()
    }
}
