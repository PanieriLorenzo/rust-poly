use core::f64;

use crate::{
    aliases::{C, R},
    base::{BasePoly, UnivariateMarker},
    errors::CAST_OVERFLOW,
    scalar_traits::{BasicScalar, ComplexScalar, NonIntegerScalar},
    storage_traits::BasicStorage,
};

use fastrand::Rng;
use num_traits::{FromPrimitive, One, cast};

pub fn initial_guess_uniform<T: ComplexScalar>(
    seed: u64,
    min_re: T::RealPartScalar,
    max_re: T::RealPartScalar,
    min_im: T::RealPartScalar,
    max_im: T::RealPartScalar,
    out: &mut Vec<T>,
    num: usize,
) {
    let mut rng = fastrand::Rng::with_seed(seed);
    let re_span = max_re.clone() - min_re.clone();
    let im_span = max_im.clone() - min_im.clone();
    for _ in 0..num {
        let re = T::RealPartScalar::from_f64(rng.f64()).expect(CAST_OVERFLOW) * re_span.clone()
            + min_re.clone();
        let im = T::RealPartScalar::from_f64(rng.f64()).expect(CAST_OVERFLOW) * im_span.clone()
            + min_im.clone();
        out.push(T::cartesian(re, im));
    }
}

pub fn initial_guess_annulus<S: BasicStorage>(
    seed: u64,
    poly: &BasePoly<S, UnivariateMarker>,
    bias: R<S>,
    perturbation: R<S>,
    out: &mut Vec<C<S>>,
    num: usize,
) where
    S::T: NonIntegerScalar,
{
    // ensure that generated guesses have a "gap" so that guesses are always radially
    // asymmetric. This is important as some methods may get stuck if the guesses
    // are symmetrical.
    let n_asym = if num % 2 == 0 { num + 1 } else { num + 2 };
    let mut rng = Rng::with_seed(seed);
    let angle_increment = std::f64::consts::TAU / cast(n_asym).unwrap_or(f64::MAX);
    let low: R<S> = poly.lower_bound_inner();
    let high: R<S> = poly.upper_bound_inner();
    let span: R<S> = high.clone() - low.clone();
    let radius: R<S> = high * bias.clone() + low.clone() * (R::<S>::one() - bias);
    let mut angle_accumulator = 0.0;
    for _ in 0..num {
        let angle: R<S> = R::<S>::from_f64(angle_accumulator).expect(CAST_OVERFLOW)
            + R::<S>::from_f64(rng.f64().mul_add(angle_increment, -(angle_increment / 2.0)))
                .expect(CAST_OVERFLOW)
                * perturbation.clone();
        let radius = radius.clone() * (R::<S>::one() - perturbation.clone())
            + (R::<S>::from_f64(rng.f64()).expect(CAST_OVERFLOW) * span.clone() + low.clone())
                * perturbation.clone();
        out.push(C::<S>::polar(radius, angle));
        angle_accumulator += angle_increment;
    }
}
