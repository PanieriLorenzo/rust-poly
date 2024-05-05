use duplicate::duplicate_item;
use nalgebra::RealField;

use num::{FromPrimitive, Num};
use std::ops::{AddAssign, DivAssign, MulAssign, Neg, RemAssign, SubAssign};

use numeric_constant_traits::{Eight, Four, Three, Two};

/// The trait bounds necessary to provide the basic functionality of this crate.
pub trait Scalar: Clone + PartialEq + std::fmt::Debug + Num + 'static {}
impl<T: Clone + PartialEq + std::fmt::Debug + Num + 'static> Scalar for T {}

// TODO: these are required by nalgebra for things that shouldn't require them.
//       perhaps in the future they can be dropped?
/// Trait bounds necessary to provide more advanced mathematical features.
pub trait ScalarOps: Scalar + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign {}
impl<T: Scalar + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign> ScalarOps for T {}

#[deprecated]
pub trait FloatScalar:
    Scalar
    + Neg
    + Num
    + Clone
    + DivAssign
    + RemAssign
    + SubAssign
    + MulAssign
    + AddAssign
    + core::fmt::Debug
    + 'static
    + RealField
    + PartialOrd
    + Two
    + Three
    + Four
    + Eight
    + FromPrimitive
{
}

// TODO: commented out types requires nalgebra to relax its requirements
#[duplicate_item(
    scalar_type;
    // [ i8 ];
    // [ i16 ];
    // [ i32 ];
    // [ i64 ];
    // [ i128 ];
    // [ isize ];
    [ f32 ];
    [ f64 ];
    // [ Complex<i8> ];
    // [ Complex<i16> ];
    // [ Complex<i32> ];
    // [ Complex<i64> ];
    // [ Complex<i128> ];
    // [ Complex<isize> ];
    // [ Complex<u8> ];
    // [ Complex<u16> ];
    // [ Complex<u32> ];
    // [ Complex<u64> ];
    // [ Complex<u128> ];
    // [ Complex<usize> ];
    // [ Complex<f32> ];
    // [ Complex<f64> ];
    // [ Ratio<i8> ];
    // [ Ratio<i16> ];
    // [ Ratio<i32> ];
    // [ Ratio<i64> ];
    // [ Ratio<i128> ];
    // [ Ratio<isize> ];
)]
impl FloatScalar for scalar_type {}
