//! Traits for the coefficients of polynomials

use na::ComplexField;
use num::{
    complex::{Complex32, Complex64},
    traits::{MulAdd, MulAddAssign},
    Complex, FromPrimitive, Num,
};
use std::ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign};

/// The trait bounds necessary to provide the basic functionality of this crate.
pub trait Scalar:
    SafeConstants
    + Clone
    + PartialEq
    + std::fmt::Debug
    + Num
    + FromPrimitive
    + std::ops::Neg<Output = Self>
    + 'static
{
}
impl<
        T: Clone
            + PartialEq
            + std::fmt::Debug
            + Num
            + FromPrimitive
            + SafeConstants
            + std::ops::Neg<Output = Self>
            + 'static,
    > Scalar for T
{
}

pub trait ComplexScalar: Scalar {
    type ComponentScalar;
}

impl ComplexScalar for Complex64 {
    type ComponentScalar = f64;
}
impl ComplexScalar for Complex32 {
    type ComponentScalar = f32;
}

pub trait RealScalar: Scalar {}

impl RealScalar for f32 {}
impl RealScalar for f64 {}

// TODO: these are required by nalgebra for things that shouldn't require them.
//       perhaps in the future they can be dropped?
/// Trait bounds necessary to provide more advanced mathematical features.
#[allow(clippy::module_name_repetitions)]
pub trait ScalarOps:
    Scalar
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
    + MulAdd<Output = Self>
    + MulAddAssign
{
}
impl<
        T: Scalar
            + AddAssign
            + SubAssign
            + MulAssign
            + DivAssign
            + RemAssign
            + MulAdd<Output = Self>
            + MulAddAssign,
    > ScalarOps for T
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
        Complex {
            re: T::tiny_safe(),
            im: T::tiny_safe(),
        }
    }

    fn small_safe() -> Self {
        Complex {
            re: T::small_safe(),
            im: T::small_safe(),
        }
    }

    fn large_safe() -> Self {
        Complex {
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
