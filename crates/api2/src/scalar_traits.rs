use std::{
    fmt::Debug,
    ops::{Add, Div, Mul, Sub},
    process::Output,
};

use num_complex::{Complex32, Complex64};
use num_traits::{FromPrimitive, One, Zero};

use crate::errors::{CAST_OVERFLOW, INFALLIBLE_CONVERSION};

pub trait BasicScalar:
    Clone
    + Add<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Zero
    + One
    + PartialEq
    + Debug
    + FromPrimitive
{
    /// The type of the real part of a number. If the number is already real,
    /// this should just be `Self`.
    type RealPartScalar: RealScalar;

    /// Return the real part of a number. If the number is already real, this
    /// just returns itself.
    fn real_part(&self) -> &Self::RealPartScalar;

    fn taxicab_norm(&self) -> Self::RealPartScalar;
}

impl BasicScalar for f32 {
    type RealPartScalar = f32;

    fn real_part(&self) -> &Self::RealPartScalar {
        self
    }

    fn taxicab_norm(&self) -> Self::RealPartScalar {
        *self
    }
}

impl BasicScalar for f64 {
    type RealPartScalar = f64;

    fn real_part(&self) -> &Self::RealPartScalar {
        self
    }

    fn taxicab_norm(&self) -> Self::RealPartScalar {
        *self
    }
}

impl BasicScalar for Complex32 {
    type RealPartScalar = f32;

    fn real_part(&self) -> &Self::RealPartScalar {
        &self.re
    }

    fn taxicab_norm(&self) -> Self::RealPartScalar {
        self.norm_sqr()
    }
}

impl BasicScalar for Complex64 {
    type RealPartScalar = f64;

    fn real_part(&self) -> &Self::RealPartScalar {
        &self.re
    }

    fn taxicab_norm(&self) -> Self::RealPartScalar {
        self.norm_sqr()
    }
}

/// A scalar that is not an integer, but may be either a complex or a real.
///
/// This is used for all float-like operations that do not require [`PartialOrd`].
pub trait NonIntegerScalar: BasicScalar {
    /// The type of a complex number whose real part is this type. For a type
    /// that is already complex this is just `Self`.
    type ToComplexScalar: ComplexScalar;

    fn sqrt(&self) -> Self;

    fn to_complex(&self) -> Self::ToComplexScalar;

    /// Override to provide a more accurate or fast specialization for the quadratic formula.
    fn quadratic_formula(a: Self, b: Self, c: Self) -> (Self, Self) {
        let four = Self::from_u8(4).expect(CAST_OVERFLOW);
        let two = Self::from_u8(2).expect(CAST_OVERFLOW);

        // TODO: switch to different formula when b^2 and 4c are very close due
        //       to loss of precision
        let plus_minus_term = (b.clone() * b.clone() - four * a.clone() * c.clone()).sqrt();
        let x1 = (plus_minus_term.clone() - b.clone()) / (two.clone() * a.clone());
        let x2 = ((Self::zero() - b.clone()) - plus_minus_term) / (two * a.clone());
        (x1, x2)
    }
}

impl NonIntegerScalar for f32 {
    type ToComplexScalar = Complex32;

    fn sqrt(&self) -> Self {
        f32::sqrt(*self)
    }

    fn to_complex(&self) -> Self::ToComplexScalar {
        Complex32::from_f32(*self).expect(INFALLIBLE_CONVERSION)
    }
}

impl NonIntegerScalar for f64 {
    type ToComplexScalar = Complex64;

    fn sqrt(&self) -> Self {
        f64::sqrt(*self)
    }

    fn to_complex(&self) -> Self::ToComplexScalar {
        Complex64::from_f64(*self).expect(INFALLIBLE_CONVERSION)
    }
}

impl NonIntegerScalar for Complex32 {
    type ToComplexScalar = Self;

    fn sqrt(&self) -> Self {
        Complex32::sqrt(*self)
    }

    fn to_complex(&self) -> Self::ToComplexScalar {
        *self
    }
}

impl NonIntegerScalar for Complex64 {
    type ToComplexScalar = Self;

    fn sqrt(&self) -> Self {
        Complex64::sqrt(*self)
    }

    fn to_complex(&self) -> Self::ToComplexScalar {
        *self
    }
}

pub trait RealScalar: NonIntegerScalar + PartialOrd {
    /// Smallest number that can be safely used in reciprocals without causing
    /// a division by zero error, NaN, infinite or similar.
    const TINY: Self;

    /// Is smaller than or equal to [`Self::TINY`]
    fn is_tiny(&self) -> bool {
        self <= &Self::TINY
    }
}

impl RealScalar for f32 {
    const TINY: Self = Self::MIN_POSITIVE;
}

impl RealScalar for f64 {
    const TINY: Self = Self::MIN_POSITIVE;
}

pub trait ComplexScalar: NonIntegerScalar {
    fn new(re: Self::RealPartScalar, im: Self::RealPartScalar) -> Self;
}

impl ComplexScalar for Complex32 {
    fn new(re: Self::RealPartScalar, im: Self::RealPartScalar) -> Self {
        Complex32::new(re, im)
    }
}

impl ComplexScalar for Complex64 {
    fn new(re: Self::RealPartScalar, im: Self::RealPartScalar) -> Self {
        Complex64::new(re, im)
    }
}
