use crate::{
    __util,
    num::{Complex, Zero},
    poly::roots,
    roots::{line_search_accelerate, line_search_decelerate, LazyDerivatives},
    scalar::SafeConstants,
    Poly, Scalar, ScalarOps,
};
use na::RealField;
use num::One;

/// Modified Newton's method, originally conceived by [Kaj Madsen 1973](https://doi.org/10.1007/BF01933524).
///
/// This method is a much more robust than traditional naive Newton iteration.
/// It can detect when convergence slows down and increase the step size. It can
/// also detect when it gets stuck in a local minimum or maximum and unstuck
/// itself.
///
/// This implementation was based on [Henrik Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
pub fn newton<T: ScalarOps + RealField>(
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
    let mut guess_delta = Complex::one();

    // number of iterations without improvements after which we assume we're in
    // a cycle set it to a prime number to avoid meta-cycles forming
    const CYCLE_COUNT_THRESHOLD: usize = 17;
    let mut cycle_counter = 0;
    let mut best_guess = guess;
    let mut best_px_norm = guess.norm();

    let mut diffs = LazyDerivatives::new(poly);

    // until convergence
    for i in __util::iterator::saturating_counter() {
        let px = poly.eval(guess);
        eval_counter += 1;

        // stopping criterion 1: converged
        if px.norm() <= epsilon {
            return Ok((vec![best_guess], eval_counter));
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3 && super::stopping_criterion_garwick(guess, guess_old, guess_old_old) {
            return Ok((vec![best_guess], eval_counter));
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            log::trace!("did not converge {{best_guess: {best_guess}, poly: {poly}}}");
            return Err(roots::Error::NoConverge(vec![best_guess]));
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
                guess = best_guess - guess_delta * backoff;
            }
        } else {
            cycle_counter = 0;
            best_guess = guess;
            best_px_norm = px.norm();
        }

        let pdx = diffs.get_nth_derivative(1).eval(guess);
        eval_counter += 1;

        guess_delta = compute_delta(px, pdx, guess_delta);

        // naive newton step, we're gonna improve this.
        let mut guess_new = guess - guess_delta;
        let px_new = poly.eval(guess_new);
        let pdx_new = diffs.get_nth_derivative(1).eval(guess_new);

        if !check_will_converge(guess_new, guess, px_new, pdx_new, pdx) {
            // if the current guess isn't captured, we adjust our guess more
            // aggressively until it is captured, i.e. it is close enough that
            // it "falls" towards the correct guess instead of jumping somewhere
            // else.
            if px_new.norm() > px.norm() {
                let res = line_search_decelerate(poly, guess, guess_delta);
                guess_new = res.0;
                eval_counter += res.1;
            } else {
                let res = line_search_accelerate(poly, guess, guess_delta);
                guess_new = res.0;
                eval_counter += res.1;
            }
        }

        guess_old_old = guess_old;
        guess_old = guess;
        guess = guess_new;
    }
    unreachable!()
}

fn compute_delta<T: Scalar>(px: Complex<T>, pdx: Complex<T>, delta_old: Complex<T>) -> Complex<T> {
    const EXPLODE_THRESHOLD: f64 = 5.0;

    if pdx.is_zero() {
        // these are arbitrary, originally chosen by Madsen.
        const ROTATION_RADIANS: f64 = 0.925_024_5;
        const SCALE: f64 = 5.0;
        // TODO: when const trait methods are supported, this should be
        //       made fully const.
        let backoff = Complex::from_polar(
            T::from_f64(SCALE).expect("overflow"),
            T::from_f64(ROTATION_RADIANS).expect("overflow"),
        );
        return delta_old * backoff;
    }

    let delta = px / pdx;

    if delta.norm() > delta_old.norm() * T::from_f64(EXPLODE_THRESHOLD).expect("overflow") {
        // these are arbitrary, originally chosen by Madsen.
        const ROTATION_RADIANS: f64 = 0.925_024_5;
        const SCALE: f64 = 5.0;
        // TODO: when const trait methods are supported, this should be
        //       made fully const.

        let backoff = Complex::from_polar(
            T::from_f64(SCALE).expect("overflow"),
            T::from_f64(ROTATION_RADIANS).expect("overflow"),
        )
        .scale(delta_old.norm() / delta.norm());

        return delta * backoff;
    }

    delta
}

/// Heuristic that tries to determine if the current guess is close enough to the
/// real root to be "captured", i.e. if the current guess is guaranteed to
/// converge to this root, rather than jumping off to a different root.
///
/// This condition is based on [Ostrowski 1966](https://doi.org/10.2307/2005025),
/// but uses an approximation by [Henrik Vestermark 2020](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
fn check_will_converge<T: Scalar>(
    guess: Complex<T>,
    guess_old: Complex<T>,
    px: Complex<T>,
    pdx: Complex<T>,
    pdx_old: Complex<T>,
) -> bool {
    !(px * pdx).is_small()
        && ((px / pdx).norm()
            * T::from_u8(2).expect("overflow")
            * ((pdx_old - pdx) / (guess_old - guess)).norm()
            <= pdx.norm())
}

#[cfg(test)]
mod test {
    use super::newton;
    use crate::{__util::testing::check_roots, num::One, Poly64};

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-6));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-8));
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
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
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
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-6));
    }
}
