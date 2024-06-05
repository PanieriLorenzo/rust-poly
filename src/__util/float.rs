//! floating point utilities

use softfloat::{F32, F64};

pub const F64_SMALL_NUM: f64 = f64_const_sqrt(f64::MIN_POSITIVE) / f64::EPSILON;
pub const F64_BIG_NUM: f64 = 1.0 / F64_SMALL_NUM;

pub const F32_SMALL_NUM: f32 = f32_const_sqrt(f32::MIN_POSITIVE) / f32::EPSILON;
pub const F32_BIG_NUM: f32 = 1.0 / F32_SMALL_NUM;

pub const fn f32_const_sqrt(x: f32) -> f32 {
    F32::from_native_f32(x).to_native_f32()
}

pub const fn f64_const_sqrt(x: f64) -> f64 {
    F64::from_native_f64(x).to_native_f64()
}

/// Makes a degenerate float normal again by either clamping it or replacing
/// NaN with zero.
pub fn f64_make_normal(x: f64) -> f64 {
    if x.is_nan() {
        return 0.0;
    }

    if x.is_infinite() && x.is_sign_positive() {
        return f64::MAX;
    }

    if x.is_infinite() && x.is_sign_negative() {
        return f64::MIN;
    }

    if x.is_subnormal() {
        return 0.0;
    }

    debug_assert!(
        x.is_normal(),
        "all denormal cases are checked by this point"
    );
    x
}

/// If a float is subnormal or zero, bump it to the nearest normal number or
/// `MIN_POSITIVE` if it's zero.
pub fn f64_make_nonzero(x: f64) -> f64 {
    let x = f64_make_normal(x);
    if x.abs() < f64::MIN_POSITIVE {
        f64::MIN_POSITIVE.copysign(x)
    } else {
        x
    }
}

/// If a float's absolute value is outside of the safe range
/// `[F64_SMALL_NUM, F64_BIG_NUM]`, bump it to the nearest safe value.
/// This also bumps zeros, as they are inherently unsafe, and NaNs are
/// treated as `F64_SMALL_NUM`.
pub fn f64_make_safe(x: f64) -> f64 {
    let x = f64_make_normal(x);
    if x.abs() < F64_SMALL_NUM {
        F64_SMALL_NUM.copysign(x)
    } else {
        x
    }
}
