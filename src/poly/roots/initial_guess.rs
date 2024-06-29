use itertools::Itertools;
use num::FromPrimitive;

use crate::{
    num::{Complex, Float, Zero},
    Poly, Scalar, ScalarOps,
    __util::complex::{c_min, c_neg},
};

impl<T: ScalarOps + Float> Poly<T> {
    /// [ref](https://doi.org/10.1007/BF01933524)
    pub(crate) fn initial_guess_smallest(&self) -> Complex<T> {
        debug_assert!(self.is_normalized());
        debug_assert!(self.len_raw() >= 2);

        let small = Float::recip(T::from_u16(1_000).expect("overflow"));
        let p_diff = self.clone().diff();
        let mut pz = self.eval(Complex::zero());
        let mut pdz = p_diff.eval(Complex::zero());

        // avoid divide by zero
        if pdz.norm() < small {
            pz += small;
            pdz += small;
        }

        let theta = (c_neg(pz) / pdz).arg();
        let mut iter_coeffs = self.0.iter();
        let a0 = iter_coeffs.next().expect("infallible");

        let mut guess = iter_coeffs
            .zip(1..)
            .map(|(ak, k)| {
                Complex::i()
                    .scale(theta)
                    .exp()
                    .scale((a0 / ak).norm())
                    .powf(T::one() / T::from_usize(k).expect("overflow"))
            })
            .reduce(c_min)
            .expect("infallible")
            .scale(Float::recip(T::from_u8(2).expect("overflow")));

        if guess.im.is_zero() {
            // add a small constant because some methods can't converge to
            // complex roots if the initial guess is real
            guess += Complex::i().scale(Float::recip(T::from_u16(1_000).expect("overflow")));
        }
        guess
    }

    pub(crate) fn initial_guess_lower_bound(&self) -> T {
        todo!()
    }

    pub(crate) fn initial_guess_uppwer_bound(&self) -> T {
        todo!()
    }
}

pub fn initial_guesses_random<T: Scalar>(poly: &Poly<T>, seed: u64, out: &mut [Complex<T>]) {
    let mut rng = fastrand::Rng::with_seed(seed);
    let low = lower_bound(poly).to_f64().expect("overflow");
    let high = upper_bound(poly).to_f64().expect("overflow");
    let span = high - low;
    for y in out {
        let radius = T::from_f64(rng.f64() * span + low).expect("overflow");
        let angle = T::from_f64(rng.f64() * std::f64::consts::TAU).expect("overflow");
        *y = Complex::from_polar(radius, angle);
    }
}

/// The radius of a disk containing all the roots
///
/// Uses Deutsch's simple formula \[[McNamee 2005](https://www.researchgate.net/publication/228745231_A_comparison_of_a_priori_bounds_on_real_or_complex_roots_of_polynomials)\]
fn upper_bound<T: Scalar>(poly: &Poly<T>) -> T {
    debug_assert!(
        poly.degree_raw() >= 1,
        "there are no bounds for a polynomial with no roots"
    );
    debug_assert!(
        poly.is_monic(),
        "Deuthsch's formula requires the polynomial to be monic"
    );

    let n = poly.len_raw();

    let next_last = poly.0[poly.len_raw() - 2];
    let coeffs_iter = poly.0.iter().take(n - 2);
    let coeffs_iter_shifted = poly.0.iter().skip(1).take(n - 2);
    let max_term = coeffs_iter
        .zip(coeffs_iter_shifted)
        .map(|(num, denom)| num / denom)
        .map(|z| z.norm())
        .reduce(|acc, z| if z > acc { z } else { acc })
        .expect("infallible");
    next_last.norm() + max_term
}

/// The radius of a disk containing none of the roots
fn lower_bound<T: Scalar>(poly: &Poly<T>) -> T {
    let mut this = Poly::from_complex_vec(poly.0.iter().cloned().rev().collect_vec());
    this.make_monic();
    upper_bound(&this).recip()
}
