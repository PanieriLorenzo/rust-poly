use super::F128;
use softfloat::F64;

/// Dummy softfloat F64 for the constants
struct F64Consts;

impl F64Consts {
    const BITS: u32 = 64;
    const SIGNIFICAND_BITS: u32 = 52;
    const EXPONENT_BITS: u32 = Self::BITS - Self::SIGNIFICAND_BITS - 1;
    const EXPONENT_MAX: u32 = (1 << Self::EXPONENT_BITS) - 1;
    const EXPONENT_BIAS: u32 = Self::EXPONENT_MAX >> 1;
    const SIGN_MASK: u64 = 1 << (Self::BITS - 1);
    const SIGNIFICAND_MASK: u64 = (1 << Self::SIGNIFICAND_BITS) - 1;
    const IMPLICIT_BIT: u64 = 1 << Self::SIGNIFICAND_BITS;
    const EXPONENT_MASK: u64 = !(Self::SIGN_MASK | Self::SIGNIFICAND_MASK);
}

// Source: https://github.com/rust-lang/compiler-builtins/blob/3dea633a80d32da75e923a940d16ce98cce74822/src/float/extend.rs#L4
pub const fn extend(a: F64) -> F128 {
    let src_zero = 0;
    let src_one = 1;
    let src_bits = F64Consts::BITS;
    let src_sign_bits = F64Consts::SIGNIFICAND_BITS;
    let src_exp_bias = F64Consts::EXPONENT_BIAS;
    let src_min_normal = F64Consts::IMPLICIT_BIT;
    let src_infinity = F64Consts::EXPONENT_MASK;
    let src_sign_mask = F64Consts::SIGN_MASK;
    let src_abs_mask = src_sign_mask - src_one;
    let src_qnan = F64Consts::SIGNIFICAND_MASK;
    let src_nan_code = src_qnan - src_one;

    let dst_bits = F128::BITS;
    let dst_sign_bits = F128::SIGNIFICAND_BITS;
    let dst_inf_exp = F128::EXPONENT_MAX as u64;
    let dst_exp_bias = F128::EXPONENT_BIAS;
    let dst_min_normal = F128::IMPLICIT_BIT;

    let sign_bits_delta = dst_sign_bits - src_sign_bits;
    let exp_bias_delta = dst_exp_bias - src_exp_bias;
    let a_abs = a.to_bits() & src_abs_mask;
    let mut abs_result: u128 = 0;

    if a_abs.wrapping_sub(src_min_normal) < src_infinity.wrapping_sub(src_min_normal) {
        // a is a normal number.
        // Extend to the destination type by shifting the significand and
        // exponent into the proper position and rebiasing the exponent.
        let abs_dst = a_abs as u128;
        let bias_dst = exp_bias_delta as u128;
        abs_result = abs_dst.wrapping_shl(sign_bits_delta);
        abs_result += bias_dst.wrapping_shl(dst_sign_bits);
    } else if a_abs >= src_infinity {
        // a is NaN or infinity.
        // Conjure the result by beginning with infinity, then setting the qNaN
        // bit (if needed) and right-aligning the rest of the trailing NaN
        // payload field.
        let qnan_dst = (a_abs & src_qnan) as u128;
        let nan_code_dst = (a_abs & src_nan_code) as u128;
        let inf_exp_dst = dst_inf_exp as u128;
        abs_result = inf_exp_dst.wrapping_shl(dst_sign_bits);
        abs_result |= qnan_dst.wrapping_shl(sign_bits_delta);
        abs_result |= nan_code_dst.wrapping_shl(sign_bits_delta);
    } else if a_abs != src_zero {
        // a is denormal.
        // Renormalize the significand and clear the leading bit, then insert
        // the correct adjusted exponent in the destination type.
        let scale = a_abs.leading_zeros() - src_min_normal.leading_zeros();
        let abs_dst = a_abs as u128;
        let bias_dst = (exp_bias_delta - scale + 1) as u128;
        abs_result = abs_dst.wrapping_shl(sign_bits_delta + scale);
        abs_result = (abs_result ^ dst_min_normal) | (bias_dst.wrapping_shl(dst_sign_bits));
    }

    let sign_result = (a.to_bits() & src_sign_mask) as u128;
    F128::from_repr(abs_result | (sign_result.wrapping_shl(dst_bits - src_bits)))
}

// Source: https://github.com/rust-lang/compiler-builtins/blob/3dea633a80d32da75e923a940d16ce98cce74822/src/float/trunc.rs#L4
pub const fn trunc(a: F128) -> F64 {
    let src_zero = 0_u128;
    let src_one = 1_u128;
    let src_bits = F128::BITS;
    let src_exp_bias = F128::EXPONENT_BIAS as u64;

    let src_min_normal = F128::IMPLICIT_BIT;
    let src_significand_mask = F128::SIGNIFICAND_MASK;
    let src_infinity = F128::EXPONENT_MASK;
    let src_sign_mask = F128::SIGN_MASK;
    let src_abs_mask = src_sign_mask - src_one;
    let round_mask = (src_one << (F128::SIGNIFICAND_BITS - F64Consts::SIGNIFICAND_BITS)) - src_one;
    let halfway = src_one << (F128::SIGNIFICAND_BITS - F64Consts::SIGNIFICAND_BITS - 1);
    let src_qnan = src_one << (F128::SIGNIFICAND_BITS - 1);
    let src_nan_code = src_qnan - src_one;

    let dst_zero = 0_u64;
    let dst_one = 1_u64;
    let dst_bits = F64Consts::BITS;
    let dst_inf_exp = F64Consts::EXPONENT_MAX as u64;
    let dst_exp_bias = F64Consts::EXPONENT_BIAS as u64;

    let underflow_exponent: u128 = (src_exp_bias + 1 - dst_exp_bias) as u128;
    let overflow_exponent: u128 = (src_exp_bias + dst_inf_exp as u64 - dst_exp_bias) as u128;
    let underflow: u128 = underflow_exponent << F128::SIGNIFICAND_BITS;
    let overflow: u128 = overflow_exponent << F128::SIGNIFICAND_BITS;

    let dst_qnan = 1_u64 << (F64Consts::SIGNIFICAND_BITS - 1);
    let dst_nan_code = dst_qnan - dst_one;

    let sign_bits_delta = F128::SIGNIFICAND_BITS - F64Consts::SIGNIFICAND_BITS;
    // Break a into a sign and representation of the absolute value.
    let a_abs = a.repr() & src_abs_mask;
    let sign = a.repr() & src_sign_mask;
    let mut abs_result: u64;

    if a_abs.wrapping_sub(underflow) < a_abs.wrapping_sub(overflow) {
        // The exponent of a is within the range of normal numbers in the
        // destination format.  We can convert by simply right-shifting with
        // rounding and adjusting the exponent.
        abs_result = (a_abs >> sign_bits_delta) as u64;
        let tmp = src_exp_bias.wrapping_sub(dst_exp_bias) << F64Consts::SIGNIFICAND_BITS;
        abs_result = abs_result.wrapping_sub(tmp as u64);

        let round_bits = a_abs & round_mask;
        if round_bits > halfway {
            // Round to nearest.
            abs_result += dst_one;
        } else if round_bits == halfway {
            // Tie to even.
            abs_result += abs_result & dst_one;
        };
    } else if a_abs > src_infinity {
        // a is NaN.
        // Conjure the result by beginning with infinity, setting the qNaN
        // bit and inserting the (truncated) trailing NaN field.
        abs_result = (dst_inf_exp << F64Consts::SIGNIFICAND_BITS) as u64;
        abs_result |= dst_qnan;
        abs_result |= dst_nan_code
            & ((a_abs & src_nan_code) >> (F128::SIGNIFICAND_BITS - F64Consts::SIGNIFICAND_BITS))
                as u64;
    } else if a_abs >= overflow {
        // a overflows to infinity.
        abs_result = (dst_inf_exp << F64Consts::SIGNIFICAND_BITS) as u64;
    } else {
        // a underflows on conversion to the destination type or is an exact
        // zero.  The result may be a denormal or zero.  Extract the exponent
        // to get the shift amount for the denormalization.
        let a_exp: u32 = (a_abs >> F128::SIGNIFICAND_BITS) as u32;
        let shift = src_exp_bias - dst_exp_bias - a_exp as u64 + 1;

        let significand = (a.repr() & src_significand_mask) | src_min_normal;

        // Right shift by the denormalization amount with sticky.
        if shift > F64Consts::SIGNIFICAND_BITS as u64 {
            abs_result = dst_zero;
        } else {
            let sticky = if (significand << (src_bits as u64 - shift)) != src_zero {
                src_one
            } else {
                src_zero
            };
            let denormalized_significand: u128 = significand >> shift | sticky;
            abs_result = (denormalized_significand
                >> (F128::SIGNIFICAND_BITS - F64Consts::SIGNIFICAND_BITS))
                as u64;
            let round_bits = denormalized_significand & round_mask;
            // Round to nearest
            if round_bits > halfway {
                abs_result += dst_one;
            }
            // Ties to even
            else if round_bits == halfway {
                abs_result += abs_result & dst_one;
            };
        }
    }

    // Apply the signbit to the absolute value.
    F64::from_bits(abs_result | sign.wrapping_shr(src_bits - dst_bits) as u64)
}

#[cfg(test)]
mod test {
    use core::f64;

    use super::{extend, trunc};
    use softfloat::F64;

    #[test]
    fn back_and_forth() {
        let before = -1234.5678;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);

        let before = f64::INFINITY;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);

        let before = -f64::INFINITY;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);

        let before = f64::EPSILON;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);

        let before = 0.0;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);

        let before = -0.0;
        let after = trunc(extend(F64::from_native_f64(before))).to_native_f64();
        assert_eq!(before, after);
    }
}
