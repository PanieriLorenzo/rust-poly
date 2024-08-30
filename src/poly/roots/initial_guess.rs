use std::cmp::Ordering;

use itertools::Itertools;

use crate::{
    num::{Complex, Zero},
    scalar::Rational,
    util::{
        casting::usize_to_f64,
        complex::{c_arg, c_exp, c_from_f64, c_min, c_neg, c_powf},
        doc_macros::{panic_t_from_f64, panic_t_from_int},
    },
    Poly, RealScalar,
};

// TODO: initial guesses Bini (see allow(unused) functions below)

/// Guess close to the root with the smallest magnitude ([Madsen 1973](https://doi.org/10.1007/BF01933524))
///
/// # Panics
#[doc = panic_t_from_int!(r"usize")]
#[allow(clippy::module_name_repetitions)]
pub fn initial_guess_smallest<T: RealScalar>(poly: &Poly<T>) -> Complex<T> {
    debug_assert!(poly.is_normalized());
    debug_assert!(poly.len_raw() >= 2);

    let small = Rational::recip(T::from_u16(1_000).expect("overflow"));
    let p_diff = poly.clone().diff();
    let mut pz = poly.eval(Complex::zero());
    let mut pdz = p_diff.eval(Complex::zero());

    // avoid divide by zero
    if pdz.norm_sqr() < small {
        pz += small.clone();
        pdz += small;
    }

    let theta = c_arg(c_neg(pz) / pdz);
    let mut iter_coeffs = poly.0.iter();
    let a0 = iter_coeffs.next().expect("infallible");

    let mut guess = iter_coeffs
        .zip(1..)
        .map(|(ak, k)| {
            c_powf(
                c_exp(Complex::i().scale(theta.clone())).scale((a0 / ak).norm_sqr()),
                T::one() / T::from_usize(k).expect("overflow"),
            )
        })
        .reduce(c_min)
        .expect("infallible")
        .scale(Rational::recip(T::from_u8(2).expect("overflow")));

    if guess.im.is_zero() {
        // add a small constant because some methods can't converge to
        // complex roots if the initial guess is real
        guess += Complex::i().scale(Rational::recip(T::from_u16(1_000).expect("overflow")));
    }
    guess
}

/// TODO: doc
///
/// # Panics
#[doc = panic_t_from_f64!()]
pub fn initial_guesses_random<T: RealScalar>(mut poly: Poly<T>, seed: u64, out: &mut [Complex<T>]) {
    poly.make_monic();
    let mut rng = fastrand::Rng::with_seed(seed);
    let low = lower_bound(&poly).to_f64().expect("overflow");
    let high = upper_bound(&poly).to_f64().expect("overflow");
    let span = high - low;
    for y in out {
        let radius = rng.f64().mul_add(span, low);
        let angle = rng.f64() * std::f64::consts::TAU;
        *y = c_from_f64(Complex::from_polar(radius, angle));
    }
}

/// Equidistant points around a circle.
///
/// The bias parameter controls the radius of the circle. A bias of 0 means the
/// lower bound is used, whereas a bias of 1 means the upper bound is used. This
/// parameter can be extrapolated.
///
/// The perturbation parameter controls the amount of randomization. At 0 no
/// randomization is applied. At 1, the radius is picked uniformly between the
/// lower bound and upper bound, and the angle has a uniform random offset of
/// `pi / n_odd` where `n_odd` is the degree of the polynomial, rounded up to
/// the next odd number. In essence, the annulus that contains all the roots is
/// partitioned into equal slices and then a guess is picked at random in each
/// of these slices. This parameter can be extrapolated above 1.
///
/// # Panics
#[doc = panic_t_from_f64!()]
pub fn initial_guesses_circle<T: RealScalar>(
    poly: &Poly<T>,
    bias: T,
    seed: u64,
    perturbation: T,
    out: &mut [Complex<T>],
) {
    let n = out.len();

    // ensuring n is odd makes the points always asymmetrical, even with even
    // number of roots. This is important as some methods may get stuck if the
    // guesses are symmetrical.
    let n_odd = if n % 2 == 0 { n + 1 } else { n };

    let mut rng = fastrand::Rng::with_seed(seed);
    let angle_increment = std::f64::consts::TAU / (usize_to_f64(n_odd));
    let low = lower_bound(poly);
    let high = upper_bound(poly);
    let span = high.clone() - low.clone();
    let radius = high * bias.clone() + low.clone() * (T::one() - bias);
    let mut angle_accumulator = 0.0;
    for y in out {
        let angle = T::from_f64(angle_accumulator).expect("overflow")
            + T::from_f64(rng.f64().mul_add(angle_increment, -(angle_increment / 2.0)))
                .expect("overflow")
                * perturbation.clone();
        let radius = radius.clone() * (T::one() - perturbation.clone())
            + (T::from_f64(rng.f64()).expect("overflow") * span.clone() + low.clone())
                * perturbation.clone();
        *y = c_from_f64(Complex::from_polar(
            radius.to_f64().expect("overflow"),
            angle.to_f64().expect("overflow"),
        ));
        angle_accumulator += angle_increment;
    }
}

/// The radius of a disk containing all the roots
///
/// Uses Deutsch's simple formula \[[McNamee 2005](https://www.researchgate.net/publication/228745231_A_comparison_of_a_priori_bounds_on_real_or_complex_roots_of_polynomials)\]
fn upper_bound<T: RealScalar>(poly: &Poly<T>) -> T {
    debug_assert!(
        poly.degree_raw() >= 2,
        "upper bound of small degree polynomials is not supported, use explicit solver"
    );
    debug_assert!(
        poly.is_monic(),
        "Deuthsch's formula requires the polynomial to be monic"
    );

    let n = poly.len_raw();

    let next_last = poly.0[poly.len_raw() - 2].clone();
    let coeffs_iter = poly.0.iter().take(n - 2);
    let coeffs_iter_shifted = poly.0.iter().skip(1).take(n - 2);
    let max_term = coeffs_iter
        .zip(coeffs_iter_shifted)
        .map(|(num, denom)| num / denom)
        .map(|z| Complex::norm_sqr(&z))
        .reduce(|acc, z| if z > acc { z } else { acc })
        .expect("infallible");
    next_last.norm_sqr() + max_term
}

/// The radius of a disk containing none of the roots
fn lower_bound<T: RealScalar>(poly: &Poly<T>) -> T {
    let mut this = Poly::from_complex_vec(poly.0.iter().cloned().rev().collect_vec());
    this.make_monic();
    upper_bound(&this).recip()
}

/// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
/// Returns a positive value, if OAB makes a counter-clockwise turn,
/// negative for clockwise turn, and zero if the points are collinear.
///
/// [From Wiki Books](https://web.archive.org/web/20240617105108/https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python)
#[allow(unused)]
fn cross_2d<T: RealScalar>(o: (T, T), a: (T, T), b: (T, T)) -> T {
    (a.0 - o.0.clone()) * (b.1 - o.1.clone()) - (a.1 - o.1) * (b.0 - o.0)
}

/// Extract upper envelope of the convex hull of a set of points
///
/// [From Wiki Books](https://web.archive.org/web/20240617105108/https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python)
#[allow(unused)]
// TODO: should return an option instead of panic
fn upper_convex_envelope<T: RealScalar>(points: &mut [(T, T)]) -> Vec<(T, T)> {
    points.sort_by(
        |a, b| match (a.0.partial_cmp(&b.0), a.1.partial_cmp(&b.1)) {
            (None, _) | (Some(Ordering::Equal), None) => panic!("cannot order NaNs"),
            (Some(Ordering::Equal), Some(ord)) | (Some(ord), _) => ord,
        },
    );

    let mut upper: Vec<(T, T)> = vec![];
    for p in points.iter_mut().rev() {
        while upper.len() >= 2
            && cross_2d(
                upper[upper.len() - 2].clone(),
                upper[upper.len() - 1].clone(),
                p.clone(),
            ) <= T::zero()
        {
            upper.pop();
        }
        upper.push(p.clone());
    }
    upper
}
