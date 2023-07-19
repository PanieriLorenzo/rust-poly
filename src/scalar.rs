use duplicate::duplicate_item;
use ndarray::ScalarOperand;
use num_bigint::{BigInt, BigUint};
use num_complex::Complex;
use num_rational::{BigRational, Ratio};
use num_traits::{Float, Num};

pub trait Scalar: Num + Clone /*+ ScalarOperand*/ + core::fmt::Debug {}

#[duplicate_item(
    scalar_type;
    [ i8 ];
    [ i16 ];
    [ i32 ];
    [ i64 ];
    [ i128 ];
    [ isize ];
    [ u8 ];
    [ u16 ];
    [ u32 ];
    [ u64 ];
    [ u128 ];
    [ usize ];
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
    [ Ratio<i8> ];
    [ Ratio<i16> ];
    [ Ratio<i32> ];
    [ Ratio<i64> ];
    [ Ratio<i128> ];
    [ Ratio<isize> ];
    [ Ratio<u8> ];
    [ Ratio<u16> ];
    [ Ratio<u32> ];
    [ Ratio<u64> ];
    [ Ratio<u128> ];
    [ Ratio<usize> ];
)]
impl Scalar for scalar_type {}
