use crate::{
    casting_util::{usize_to_scalar, usize_to_u32},
    util::luts::factorial_lut,
    Poly, Scalar, ScalarOps,
};
use num::{BigUint, FromPrimitive, Num, Zero};

use crate::casting_util::{usize_to_f64, usize_to_i32};

impl<T: Scalar + FromPrimitive> Poly<T> {
    #[deprecated(note = "use cheby1 instead")]
    #[inline]
    #[must_use]
    pub fn cheby(n: usize) -> Self {
        Self::cheby1(n)
    }

    /// Get the nth [Chebyshev polynomial](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
    /// of the first kind.
    ///
    /// ```
    /// use rust_poly::{poly, Poly};
    ///
    /// assert_eq!(Poly::cheby(2), poly![-1.0, 0.0, 2.0]);
    /// assert_eq!(Poly::cheby(3), poly![0.0, -3.0, 0.0, 4.0]);
    /// assert_eq!(Poly::cheby(4), poly![1.0, 0.0, -8.0, 0.0, 8.0])
    /// ```
    #[must_use]
    pub fn cheby1(n: usize) -> Self {
        // TODO: make the first 32-ish explicit for performance
        match n {
            0 => poly![T::one()],
            1 => poly![T::zero(), T::one()],
            2 => poly![T::zero() - T::one(), T::zero(), usize_to_scalar(2)],
            3 => poly![
                T::zero(),
                T::zero() - usize_to_scalar(3),
                T::zero(),
                usize_to_scalar(4)
            ],
            4 => poly![
                T::one(),
                T::zero(),
                T::zero() - usize_to_scalar(8),
                T::zero(),
                usize_to_scalar(8)
            ],
            _ => poly![T::zero(), usize_to_scalar(2)] * Self::cheby1(n - 1) - Self::cheby1(n - 2),
        }
    }

    /// Get the nth [Bessel polynomial](https://en.wikipedia.org/wiki/Bessel_polynomials)
    #[must_use]
    pub fn bessel(n: usize) -> Option<Self> {
        let mut poly = poly![];
        for k in 0..=n {
            let c = T::from_f64(coeff(n, k))?;
            let term = Self::term(complex!(c), usize_to_u32(k));
            dbg!(&term);
            poly = poly + term;
        }
        Some(poly)
    }

    #[must_use]
    pub fn reverse_bessel(n: usize) -> Option<Self> {
        let p = Self::bessel(n)?;
        let v: Vec<_> = p.iter().cloned().rev().collect();
        Some(Self::from_complex_vec(v))
    }
}

fn factorial(n: usize) -> BigUint {
    if n <= 128 {
        return factorial_lut(n);
    }
    BigUint::from(n) * factorial(n - 1)
}

#[allow(clippy::cast_precision_loss)]
fn biguint_to_f64(x: &BigUint) -> f64 {
    let mut x_f: f64 = 0.0;
    for (i, d) in x.to_u64_digits().iter().enumerate() {
        x_f += (*d as f64) * 2.0f64.powi(usize_to_i32(i) * 64);
    }
    x_f
}

/// The coefficient for the k-th term of the n-th bessel polynomial
pub fn coeff(n: usize, k: usize) -> f64 {
    // NOTE: in theory f64 can do factorials up to 170!, but if we do it with
    //       BigUint, as we divide by factorial(n - k), we get a smaller number
    //       out, so this way we can squeeze a few more terms before f64
    //       goes to infinity than if we computed the factorials with f64 directly
    let aux_a = factorial(n + k) / factorial(n - k) / factorial(k);
    let aux_b = 1.0 / (usize_to_f64(k)).exp2();
    biguint_to_f64(&aux_a) * aux_b
}

#[cfg(test)]
#[allow(clippy::pedantic)]
mod test {
    use num::{BigUint, Num};

    use crate::poly::special_funcs::biguint_to_f64;

    use super::factorial;

    #[test]
    fn factorial_129() {
        let x = factorial(129);
        let x_expected = BigUint::from_str_radix("49745042224772874403902341504126809639656611137138843145968864022652168932196355119328515747917449637889876686464600208839390308261862352651828829226610077151044469167497022952331930501120000000000000000000000000000000", 10).unwrap();
        assert_eq!(x, x_expected);
    }

    #[test]
    fn test_biguint_to_f64() {
        let i = 1u128 << 90;
        let x = BigUint::from(i);
        assert_eq!(i as f64, biguint_to_f64(&x));
    }
}
