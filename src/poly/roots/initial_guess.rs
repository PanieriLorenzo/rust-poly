use std::cmp::Ordering;

use itertools::Itertools;

use crate::{
    num::{Complex, Float, Zero},
    Poly, RealScalar,
    __util::complex::{c_min, c_neg},
};

/// [ref](https://doi.org/10.1007/BF01933524)
pub fn initial_guess_smallest<T: RealScalar>(poly: &Poly<T>) -> Complex<T> {
    debug_assert!(poly.is_normalized());
    debug_assert!(poly.len_raw() >= 2);

    let small = Float::recip(T::from_u16(1_000).expect("overflow"));
    let p_diff = poly.clone().diff();
    let mut pz = poly.eval(Complex::zero());
    let mut pdz = p_diff.eval(Complex::zero());

    // avoid divide by zero
    if pdz.norm() < small {
        pz += small;
        pdz += small;
    }

    let theta = (c_neg(pz) / pdz).arg();
    let mut iter_coeffs = poly.0.iter();
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

pub fn initial_guesses_random<T: RealScalar>(mut poly: Poly<T>, seed: u64, out: &mut [Complex<T>]) {
    poly.make_monic();
    let mut rng = fastrand::Rng::with_seed(seed);
    let low = lower_bound(&poly).to_f64().expect("overflow");
    let high = upper_bound(&poly).to_f64().expect("overflow");
    let span = high - low;
    for y in out {
        let radius = T::from_f64(rng.f64() * span + low).expect("overflow");
        let angle = T::from_f64(rng.f64() * std::f64::consts::TAU).expect("overflow");
        *y = Complex::from_polar(radius, angle);
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
    let angle_increment = std::f64::consts::TAU / (n_odd as f64);
    let low = lower_bound(poly);
    let high = upper_bound(poly);
    let span = high - low;
    let radius = high * bias + low * (T::one() - bias);
    let mut angle_accumulator = 0.0;
    for y in out {
        let angle = T::from_f64(angle_accumulator).expect("overflow")
            + T::from_f64(rng.f64() * angle_increment - angle_increment / 2.0).expect("overflow")
                * perturbation;
        let radius = radius * (T::one() - perturbation)
            + (T::from_f64(rng.f64()).expect("overflow") * span + low) * perturbation;
        *y = Complex::from_polar(radius, angle);
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
fn cross_2d<T: RealScalar>(o: (T, T), a: (T, T), b: (T, T)) -> T {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

/// Extract upper envelope of the convex hull of a set of points
///
/// [From Wiki Books](https://web.archive.org/web/20240617105108/https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python)
fn upper_convex_envelope<T: RealScalar>(points: &mut Vec<(T, T)>) -> Vec<(T, T)> {
    points.sort_by(
        |a, b| match (a.0.partial_cmp(&b.0), a.1.partial_cmp(&b.1)) {
            (None, _) => panic!("cannot order NaNs"),
            (Some(Ordering::Equal), None) => panic!("cannot order NaNs"),
            (Some(Ordering::Equal), Some(ord)) => ord,
            (Some(ord), _) => ord,
        },
    );

    let mut upper = vec![];
    for p in points.into_iter().rev() {
        while upper.len() >= 2
            && cross_2d(upper[upper.len() - 2], upper[upper.len() - 1], *p) <= T::zero()
        {
            upper.pop();
        }
        upper.push(*p);
    }
    upper
}
