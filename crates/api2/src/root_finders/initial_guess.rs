use crate::{
    errors::CAST_OVERFLOW,
    scalar_traits::{BasicScalar, ComplexScalar, NonIntegerScalar},
};

use num_traits::FromPrimitive;

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
        out.push(T::new(re, im));
    }
}
