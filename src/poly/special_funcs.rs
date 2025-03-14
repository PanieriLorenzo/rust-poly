use crate::{
    poly2::OwnedUniPoly,
    util::{
        casting::{usize_to_f64, usize_to_i32},
        doc_macros::panic_absurd_size,
        luts::factorial_lut,
    },
    Poly, RealScalar,
};
use num::{BigUint, FromPrimitive};

impl<T: RealScalar + FromPrimitive> Poly<T> {
    #[deprecated(note = "use cheby1 instead")]
    #[inline]
    #[must_use]
    pub fn cheby(n: usize) -> Self {
        Self::cheby1(n)
    }

    /// Get the nth [Bessel polynomial](https://en.wikipedia.org/wiki/Bessel_polynomials)
    #[doc = panic_absurd_size!()]
    #[must_use]
    pub fn bessel(n: usize) -> Option<Self> {
        let mut poly = poly![];
        for k in 0..=n {
            let c = T::from_f64(coeff(n, k))?;
            let term = Self::term(complex!(c), k.try_into().expect("degree is too big"));
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

impl<T: RealScalar> Poly<T> {
    // TODO: technically it can panic in some extreme cases, would need to
    //       do some boundary testing to write a proper doc comment
    #[allow(clippy::missing_panics_doc)]
    #[must_use]
    pub fn legendre(n: usize) -> Self {
        match n {
            0 => return poly![T::one()],
            1 => return poly![T::zero(), T::one()],
            _ => {}
        }

        // this is the memoized form of the recursive recurrence relation definition
        let mut memo = Vec::with_capacity(n - 1);
        memo.push(poly![T::one()]);
        memo.push(poly![T::zero(), T::one()]);

        for i in 2..=n {
            let p1 = &memo[i - 1];
            let p2 = &memo[i - 2];
            let ns = T::from_usize(i).expect("overflow");
            memo.push(
                poly![
                    T::zero(),
                    T::from_usize(2 * i - 1).expect("overflow") / ns.clone()
                ] * p1
                    + poly![(T::one() - ns.clone()) / ns] * p2,
            );
        }
        memo.last().expect("infallible").clone()
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

    use crate::{poly::special_funcs::biguint_to_f64, Poly, Poly64};

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

    #[test]
    #[cfg_attr(miri, ignore)] // takes way too long on MIRI
    fn bessel_big() {
        // largest computable bessel polynomial
        let _ = Poly64::bessel(134).unwrap();
    }

    #[test]
    #[cfg_attr(miri, ignore)] // takes way too long on MIRI
    fn reverse_bessel_big() {
        // largest computable reverse bessel
        let _ = Poly64::reverse_bessel(134).unwrap();
    }

    #[test]
    fn legendre() {
        let p = Poly::<f32>::legendre(3);
        assert_eq!(p, poly![0.0, -1.5, 0.0, 2.5]);
    }
}
