//! floating point utilities

use crate::scalar::SafeConstants;

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
/// `[f64::small_safe(), f64::large_safe()]`, bump it to the nearest safe value.
/// This also bumps zeros, as they are inherently unsafe, and NaNs are
/// treated as `f64::small_safe()`.
pub fn f64_make_safe(x: f64) -> f64 {
    let x = f64_make_normal(x);
    if x.is_small() {
        f64::small_safe().copysign(x)
    } else if x.is_large() {
        f64::large_safe().copysign(x)
    } else {
        x
    }
}
