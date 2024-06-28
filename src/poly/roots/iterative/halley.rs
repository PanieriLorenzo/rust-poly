use crate::{
    __util,
    num::{Complex, One, Zero},
    poly::roots,
    roots::{
        line_search_accelerate, line_search_decelerate, multiplicity_lagouanelle, LazyDerivatives,
    },
    Poly, ScalarOps,
};
use na::RealField;

/// Find all roots of a polynomial using a modified Halley's method.
///
/// This implementation is based on [Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
pub fn halley<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> roots::Result<T> {
    let mut eval_counter = 0;
    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let mut roots = vec![];
    let mut initial_guesses = initial_guesses.iter().copied();

    // until we've found all roots
    loop {
        let trivial_roots = poly.trivial_roots(epsilon);
        eval_counter += trivial_roots.1;
        roots.extend(trivial_roots.0.iter());

        debug_assert!(poly.is_normalized());
        if poly.degree_raw() == 0 {
            log::debug!("{{evaluations: {eval_counter}}}");
            return Ok(roots);
        }

        let next_guess = initial_guesses.next();
        let (next_roots, num_evals) = next_root(poly, epsilon, max_iter, next_guess)?;
        let root = next_roots[0];
        eval_counter += num_evals;
        roots.push(root);
        // TODO: deflate_composite should borrow instead
        *poly = poly.clone().deflate_composite(root);
    }
}

/// Find a single root
///
/// # Returns
/// - vector of roots (usually 1)
/// - number of evaluations
fn next_root<T: ScalarOps + RealField>(
    poly: &Poly<T>,
    epsilon: T,
    max_iter: Option<usize>,
    initial_guess: Option<Complex<T>>,
) -> std::result::Result<(Vec<Complex<T>>, u128), roots::Error<Vec<Complex<T>>>> {
    let mut eval_counter = 0;
    let mut guess = initial_guess.unwrap_or_else(|| poly.initial_guess_smallest());
    let mut guess_old = guess;
    let mut guess_old_old = guess;
    let mut guess_delta_old = Complex::one();

    // number of iterations without improvements after which we assume we're in
    // a cycle set it to a prime number to avoid meta-cycles forming
    const CYCLE_COUNT_THRESHOLD: usize = 17;
    let mut cycle_counter = 0;
    let mut best_guess = guess;
    let mut best_px_norm = guess.norm();

    let mut diffs = LazyDerivatives::new(poly);

    // until convergence
    for i in __util::iterator::saturating_counter() {
        let px = poly.eval_point(guess);
        log::trace!("{{current_guess: {guess}, error: {}}}", px.norm());

        // stopping criterion 1: converged
        if px.norm() <= epsilon {
            return Ok((vec![guess], eval_counter));
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3 && super::stopping_criterion_garwick(guess, guess_old, guess_old_old) {
            return Ok((vec![guess], eval_counter));
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        // check for cycles
        if px.norm() >= best_px_norm {
            cycle_counter += 1;
            if cycle_counter > CYCLE_COUNT_THRESHOLD {
                cycle_counter = 0;
                log::trace!("cycle detected, backing off {{current_guess: {guess}, best_guess: {best_guess}}}");
                // arbitrary constants
                const ROTATION_RADIANS: f64 = 0.925_024_5;
                const SCALE: f64 = 5.0;
                // TODO: when const trait methods are supported, this should be
                //       made fully const.
                let backoff = Complex::from_polar(
                    T::from_f64(SCALE).expect("overflow"),
                    T::from_f64(ROTATION_RADIANS).expect("overflow"),
                );
                // reverting to older base guess, but offset
                guess = best_guess - guess_delta_old * backoff;
            }
        } else {
            cycle_counter = 0;
            best_guess = guess;
            best_px_norm = px.norm();
        }

        let pdx = diffs.get_nth_derivative(1).eval_point(guess);
        let pddx = diffs.get_nth_derivative(2).eval_point(guess);
        eval_counter += 2;
        let denom = (pdx * pdx).scale(T::from_u8(2).expect("overflow")) - px * pddx;

        let guess_delta = if denom.is_zero() || pdx.is_zero() {
            log::trace!("local minimum, backing off");

            // these are arbitrary, originally chosen by Madsen.
            const ROTATION_RADIANS: f64 = 0.925_024_5;
            const SCALE: f64 = 5.0;
            // TODO: when const trait methods are supported, this should be
            //       made fully const.
            let backoff = Complex::from_polar(
                T::from_f64(SCALE).expect("overflow"),
                T::from_f64(ROTATION_RADIANS).expect("overflow"),
            );
            <Complex<T> as std::ops::Mul>::mul(guess_delta_old, backoff)
        } else {
            let m = multiplicity_lagouanelle(px, pdx, pddx);
            (m + Complex::one()) * (px * pdx) / denom
        };

        const EXPLODE_THRESHOLD: f64 = 5.0;
        let guess_delta = if guess_delta.norm()
            > guess_delta_old.norm() * T::from_f64(EXPLODE_THRESHOLD).expect("overflow")
        {
            log::trace!("exploding gradient, backing off");

            // these are arbitrary, originally chosen by Madsen.
            const ROTATION_RADIANS: f64 = 0.925_024_5;
            const SCALE: f64 = 5.0;
            // TODO: when const trait methods are supported, this should be
            //       made fully const.
            let backoff = Complex::from_polar(
                T::from_f64(SCALE).expect("overflow"),
                T::from_f64(ROTATION_RADIANS).expect("overflow"),
            )
            .scale(guess_delta_old.norm() / guess_delta.norm());
            guess_delta * backoff
        } else {
            guess_delta
        };

        eval_counter += 1;
        let guess_new = if poly.eval_point(guess - guess_delta).norm() >= px.norm() {
            log::trace!("overshooting, shortening step");
            let res = line_search_decelerate(poly, guess, guess_delta);
            eval_counter += res.1;
            res.0
        } else {
            log::trace!("undershooting, lengthening step");
            let res = line_search_accelerate(poly, guess, guess_delta);
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
    use super::halley;
    use crate::{__util::testing::check_roots, num::One, Poly64};

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
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
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-4),
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
        let roots = super::halley(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }
}
