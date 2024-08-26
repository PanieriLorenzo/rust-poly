use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

use anyhow::bail;
use ethnum::u256;
use num::{
    traits::{bounds::UpperBounded, MulAdd, MulAddAssign},
    Float, FromPrimitive, Num, NumCast, One, ToPrimitive, Zero,
};
use softfloat::F64;

use crate::{
    scalar::{Rational, SafeConstants},
    RealScalar,
};

/// TODO: document that this is taken from softfloat crate
pub mod add;
pub mod cmp;
pub mod conv;
pub mod div;
pub mod mul;
pub mod pow;

#[derive(Default, Copy, Clone, Debug)]
#[repr(transparent)]
struct Bits128(u128);

#[derive(Default, Copy, Clone, Debug)]
#[repr(transparent)]
pub struct F128(Bits128);

type F = F128;
type FInt = u128;

/// Private
impl F128 {
    const BITS: u32 = 128;
    const SIGNIFICAND_BITS: u32 = 112;
    const EXPONENT_BITS: u32 = Self::BITS - Self::SIGNIFICAND_BITS - 1;
    const EXPONENT_MAX: u32 = (1 << Self::EXPONENT_BITS) - 1;
    const EXPONENT_BIAS: u32 = Self::EXPONENT_MAX >> 1;
    const SIGN_MASK: FInt = 1 << (Self::BITS - 1);
    const SIGNIFICAND_MASK: FInt = (1 << Self::SIGNIFICAND_BITS) - 1;
    const IMPLICIT_BIT: FInt = 1 << Self::SIGNIFICAND_BITS;
    const EXPONENT_MASK: FInt = !(Self::SIGN_MASK | Self::SIGNIFICAND_MASK);
    const ZERO: F128 = F128::const_from_f64(0.0);
    const ONE: F128 = F128::const_from_f64(1.0);
    const MIN_POSITIVE: F128 = F128::from_repr(0x0001_0000_0000_0000_0000_0000_0000_0000);
    const MAX: F128 = F128::from_repr(0x7ffe_ffff_ffff_ffff_ffff_ffff_ffff_ffff);

    const fn repr(self) -> FInt {
        self.0 .0
    }

    const fn from_repr(a: FInt) -> Self {
        Self(Bits128(a))
    }

    const fn signed_repr(self) -> i128 {
        self.repr() as i128
    }

    const fn normalize(significand: FInt) -> (i32, FInt) {
        let shift = significand
            .leading_zeros()
            .wrapping_sub((1u128 << Self::SIGNIFICAND_BITS).leading_zeros());
        (
            1i32.wrapping_sub(shift as i32),
            significand << shift as FInt,
        )
    }
}

/// Public
impl F128 {
    pub const fn const_from_f64(a: f64) -> Self {
        conv::extend(F64::from_native_f64(a))
    }

    pub const fn const_to_f64(self) -> f64 {
        conv::trunc(self).to_native_f64()
    }

    pub const fn const_add(self, rhs: Self) -> Self {
        add::add(self, rhs)
    }

    pub const fn const_neg(self) -> Self {
        Self::from_repr(self.repr() ^ Self::SIGN_MASK)
    }

    pub const fn const_sub(self, rhs: Self) -> Self {
        self.const_add(rhs.const_neg())
    }

    pub const fn const_cmp(self, rhs: Self) -> Option<std::cmp::Ordering> {
        cmp::cmp(self, rhs)
    }

    pub fn epsilon() -> Self {
        F128::const_from_f64(2.0).powi(-(Self::EXPONENT_BITS as i32)) / F128::const_from_f64(2.0)
    }
}

impl Add for F128 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.const_add(rhs)
    }
}

impl AddAssign for F128 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Mul for F128 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        mul::mul(self, rhs)
    }
}

impl MulAssign for F128 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl MulAdd for F128 {
    type Output = Self;

    fn mul_add(self, a: Self, b: Self) -> Self::Output {
        self * a + b
    }
}

impl MulAddAssign for F128 {
    fn mul_add_assign(&mut self, a: Self, b: Self) {
        *self = *self * a + b;
    }
}

impl Div for F128 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        div::div(self, rhs)
    }
}

impl DivAssign for F128 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Neg for F128 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.const_neg()
    }
}

impl Sub for F128 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.const_sub(rhs)
    }
}

impl SubAssign for F128 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Rem for F128 {
    type Output = Self;

    fn rem(self, _rhs: Self) -> Self::Output {
        unimplemented!()
    }
}

impl RemAssign for F128 {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs;
    }
}

impl PartialEq for F128 {
    fn eq(&self, other: &Self) -> bool {
        match self.const_cmp(*other) {
            Some(ordering) => ordering.is_eq(),
            None => false,
        }
    }
}

impl PartialOrd for F128 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.const_cmp(*other)
    }
}

impl ToPrimitive for F128 {
    fn to_i64(&self) -> Option<i64> {
        self.const_to_f64().to_i64()
    }

    fn to_u64(&self) -> Option<u64> {
        self.const_to_f64().to_u64()
    }
}

impl NumCast for F128 {
    fn from<T: num::ToPrimitive>(n: T) -> Option<Self> {
        Some(Self::const_from_f64(n.to_f64()?))
    }
}

impl Zero for F128 {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl One for F128 {
    fn one() -> Self {
        Self::ONE
    }
}

impl Num for F128 {
    type FromStrRadixErr = anyhow::Error;

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        bail!("unimplemented")
    }
}

impl Float for F128 {
    fn nan() -> Self {
        todo!()
    }

    fn infinity() -> Self {
        todo!()
    }

    fn neg_infinity() -> Self {
        todo!()
    }

    fn neg_zero() -> Self {
        todo!()
    }

    fn min_value() -> Self {
        todo!()
    }

    fn min_positive_value() -> Self {
        todo!()
    }

    fn max_value() -> Self {
        todo!()
    }

    fn is_nan(self) -> bool {
        todo!()
    }

    fn is_infinite(self) -> bool {
        todo!()
    }

    fn is_finite(self) -> bool {
        todo!()
    }

    fn is_normal(self) -> bool {
        todo!()
    }

    fn classify(self) -> std::num::FpCategory {
        todo!()
    }

    fn floor(self) -> Self {
        todo!()
    }

    fn ceil(self) -> Self {
        todo!()
    }

    fn round(self) -> Self {
        todo!()
    }

    fn trunc(self) -> Self {
        todo!()
    }

    fn fract(self) -> Self {
        todo!()
    }

    fn abs(self) -> Self {
        todo!()
    }

    fn signum(self) -> Self {
        todo!()
    }

    fn is_sign_positive(self) -> bool {
        todo!()
    }

    fn is_sign_negative(self) -> bool {
        todo!()
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        todo!()
    }

    fn recip(self) -> Self {
        Self::ONE / self
    }

    fn powi(self, n: i32) -> Self {
        pow::pow(self, n)
    }

    fn powf(self, n: Self) -> Self {
        todo!()
    }

    fn sqrt(self) -> Self {
        todo!()
    }

    fn exp(self) -> Self {
        todo!()
    }

    fn exp2(self) -> Self {
        todo!()
    }

    fn ln(self) -> Self {
        todo!()
    }

    fn log(self, base: Self) -> Self {
        todo!()
    }

    fn log2(self) -> Self {
        todo!()
    }

    fn log10(self) -> Self {
        todo!()
    }

    fn max(self, other: Self) -> Self {
        todo!()
    }

    fn min(self, other: Self) -> Self {
        todo!()
    }

    fn abs_sub(self, other: Self) -> Self {
        todo!()
    }

    fn cbrt(self) -> Self {
        todo!()
    }

    fn hypot(self, other: Self) -> Self {
        todo!()
    }

    fn sin(self) -> Self {
        todo!()
    }

    fn cos(self) -> Self {
        todo!()
    }

    fn tan(self) -> Self {
        todo!()
    }

    fn asin(self) -> Self {
        todo!()
    }

    fn acos(self) -> Self {
        todo!()
    }

    fn atan(self) -> Self {
        todo!()
    }

    fn atan2(self, other: Self) -> Self {
        todo!()
    }

    fn sin_cos(self) -> (Self, Self) {
        todo!()
    }

    fn exp_m1(self) -> Self {
        todo!()
    }

    fn ln_1p(self) -> Self {
        todo!()
    }

    fn sinh(self) -> Self {
        todo!()
    }

    fn cosh(self) -> Self {
        todo!()
    }

    fn tanh(self) -> Self {
        todo!()
    }

    fn asinh(self) -> Self {
        todo!()
    }

    fn acosh(self) -> Self {
        todo!()
    }

    fn atanh(self) -> Self {
        todo!()
    }

    fn integer_decode(self) -> (u64, i16, i8) {
        todo!()
    }
}

impl FromPrimitive for F128 {
    fn from_i64(n: i64) -> Option<Self> {
        Some(F128::const_from_f64(f64::from_i64(n)?))
    }

    fn from_u64(n: u64) -> Option<Self> {
        Some(F128::const_from_f64(f64::from_u64(n)?))
    }
}

impl SafeConstants for F128 {
    fn tiny_safe() -> Self {
        Self::MIN_POSITIVE
    }

    fn small_safe() -> Self {
        Self::MIN_POSITIVE.sqrt() / Self::epsilon()
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

impl UpperBounded for F128 {
    fn max_value() -> Self {
        Self::MAX
    }
}

impl Rational for F128 {}

fn u128_widen_mul(a: u128, b: u128) -> (u128, u128) {
    // TODO: open a PR to make wrappping_mul const
    let x = u256::wrapping_mul(u256::from(a), u256::from(b));
    (*x.low(), *x.high())
}

impl std::fmt::Display for F128 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use super::F128;

    #[test]
    fn convert() {
        let before = -1234.5678f64;
        let after = F128::const_from_f64(before).const_to_f64();
        assert_eq!(before, after);
    }
}
