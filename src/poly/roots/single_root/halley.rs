use super::{
    line_search_accelerate, line_search_decelerate, multiplicity_lagouanelle, LazyDerivatives,
};
use crate::{
    num::{Complex, One, Zero},
    poly::roots,
    roots::initial_guess::initial_guess_smallest,
    util::{
        self,
        complex::c_from_f64,
        doc_macros::{errors_no_converge, panic_t_from_f64, panic_t_from_int},
    },
    Poly, Poly2, RealScalar,
};

/// Find a single root
///
/// # Returns
/// - vector of roots (usually 1)
/// - number of evaluations
///
/// # Errors
#[doc = errors_no_converge!()]
///
/// # Panics
#[doc = panic_t_from_f64!()]
#[doc = panic_t_from_int!(r"usize")]
#[allow(clippy::items_after_statements)]
// TODO: remove this when there's a better alternative
#[allow(clippy::type_complexity)]
#[allow(clippy::similar_names)]
pub fn halley<T: RealScalar>(
    poly: &Poly<T>,
    epsilon: T,
    max_iter: Option<usize>,
    initial_guess: Option<Complex<T>>,
) -> std::result::Result<(Vec<Complex<T>>, u128), roots::Error<Vec<Complex<T>>>> {
    let mut eval_counter = 0;
    let mut guess = initial_guess.unwrap_or_else(|| initial_guess_smallest(poly));
    let mut guess_old = guess.clone();
    let mut guess_old_old = guess.clone();
    let mut guess_delta_old = Complex::one();

    // number of iterations without improvements after which we assume we're in
    // a cycle set it to a prime number to avoid meta-cycles forming
    const CYCLE_COUNT_THRESHOLD: usize = 17;
    let mut cycle_counter = 0;
    let mut best_guess = guess.clone();
    let mut best_px_norm = guess.clone().norm_sqr();

    let mut diffs = LazyDerivatives::new(poly);

    // until convergence
    for i in util::iterator::saturating_counter() {
        let px = poly.eval(guess.clone());
        log::trace!("{{current_guess: {guess}, error: {}}}", px.norm_sqr());

        // stopping criterion 1: converged
        if px.norm_sqr() <= epsilon {
            return Ok((vec![guess], eval_counter));
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3
            && super::stopping_criterion_garwick(guess.clone(), guess_old.clone(), guess_old_old)
        {
            return Ok((vec![guess], eval_counter));
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        // check for cycles
        if px.norm_sqr() >= best_px_norm {
            cycle_counter += 1;
            if cycle_counter > CYCLE_COUNT_THRESHOLD {
                // arbitrary constants
                const ROTATION_RADIANS: f64 = 0.925_024_5;
                const SCALE: f64 = 5.0;

                cycle_counter = 0;
                log::trace!("cycle detected, backing off {{current_guess: {guess}, best_guess: {best_guess}}}");

                // TODO: when const trait methods are supported, this should be
                //       made fully const.
                let backoff = c_from_f64(&Complex::from_polar(SCALE, ROTATION_RADIANS));
                // reverting to older base guess, but offset
                guess = best_guess.clone() - guess_delta_old.clone() * backoff;
            }
        } else {
            cycle_counter = 0;
            best_guess = guess.clone();
            best_px_norm = px.norm_sqr();
        }

        let pdx = diffs.get_nth_derivative(1).eval(guess.clone());
        let pddx = diffs.get_nth_derivative(2).eval(guess.clone());
        eval_counter += 2;
        let denom = (pdx.clone() * pdx.clone()).scale(T::from_u8(2).expect("overflow"))
            - px.clone() * pddx.clone();

        let guess_delta = if denom.is_zero() || pdx.is_zero() {
            // these are arbitrary, originally chosen by Madsen.
            const ROTATION_RADIANS: f64 = 0.925_024_5;
            const SCALE: f64 = 5.0;

            log::trace!("local minimum, backing off");

            // TODO: when const trait methods are supported, this should be
            //       made fully const.
            let backoff = c_from_f64(&Complex::from_polar(SCALE, ROTATION_RADIANS));
            <Complex<T> as std::ops::Mul>::mul(guess_delta_old.clone(), backoff)
        } else {
            let m = multiplicity_lagouanelle(px.clone(), pdx.clone(), pddx);
            (m + Complex::one()) * (px.clone() * pdx) / denom
        };

        const EXPLODE_THRESHOLD: f64 = 5.0;
        let guess_delta = if guess_delta.norm_sqr()
            > guess_delta_old.norm_sqr() * T::from_f64(EXPLODE_THRESHOLD).expect("overflow")
        {
            // these are arbitrary, originally chosen by Madsen.
            const ROTATION_RADIANS: f64 = 0.925_024_5;
            const SCALE: f64 = 5.0;

            log::trace!("exploding gradient, backing off");

            // TODO: when const trait methods are supported, this should be
            //       made fully const.
            let backoff = c_from_f64(&Complex::from_polar(SCALE, ROTATION_RADIANS))
                .scale(guess_delta_old.norm_sqr() / guess_delta.norm_sqr());
            guess_delta * backoff
        } else {
            guess_delta
        };

        eval_counter += 1;
        let guess_new =
            if poly.eval(guess.clone() - guess_delta.clone()).norm_sqr() >= px.norm_sqr() {
                log::trace!("overshooting, shortening step");
                let res = line_search_decelerate(poly, guess.clone(), guess_delta.clone());
                eval_counter += res.1;
                res.0
            } else {
                log::trace!("undershooting, lengthening step");
                let res = line_search_accelerate(poly, guess.clone(), guess_delta.clone());
                eval_counter += res.1;
                res.0
            };

        guess_delta_old = guess_delta;
        guess_old_old = guess_old;
        guess_old = guess;
        guess = guess_new;
    }
    unreachable!()
}

#[cfg(test)]
mod test {
    use crate::{
        num::One, roots::halley_deflate, util::__testing::check_roots, OwnedPoly, Poly2, Poly64,
    };

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-7));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-7));
    }

    #[test]
    fn degree_5_multiplicity_3() {
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(2.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-3),
            "{roots:?}"
        );
    }

    #[test]
    fn degree_5_2_zeros() {
        let roots_expected = vec![
            complex!(0.0),
            complex!(0.0),
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = halley_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-7));
    }
}
