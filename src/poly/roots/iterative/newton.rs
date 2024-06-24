use crate::{
    __util::{self, float::F64_PHI},
    num::{Complex, Float, Zero},
    roots::{line_search_accelerate, line_search_decelerate},
    scalar::SafeConstants,
    Poly,
};
use na::{ComplexField, RealField};
use num::One;

use crate::{
    poly::roots::{self, FinderConfig, FinderState, RootFinder},
    roots::FinderHistory,
    Scalar, ScalarOps,
};

use super::IterativeRootFinder;

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
    let mut initial_guesses = initial_guesses.iter().cloned();

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

    let p_diff = poly.clone().diff();

    // until convergence
    for i in __util::iterator::saturating_counter() {
        let px = poly.eval_point(guess);
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
                const ROTATION_RADIANS: f64 = 0.9250245;
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

        let pdx = p_diff.eval_point(guess);
        eval_counter += 1;

        guess_delta = compute_delta(px, pdx, guess_delta);

        // naive newton step, we're gonna improve this.
        let mut guess_new = guess - guess_delta;
        let px_new = poly.eval_point(guess_new);
        let pdx_new = p_diff.eval_point(guess_new);

        if !check_will_converge(guess_new, guess, px_new, pdx_new, pdx) {
            // if the current guess isn't captured, we adjust our guess more
            // aggressively until it is captured, i.e. it is close enough that
            // it "falls" towards the correct guess instead of jumping somewhere
            // else.
            if px_new.norm() > px.norm() {
                let res = line_search_decelerate(&poly, guess, guess_delta);
                guess_new = res.0;
                eval_counter += res.1;
            } else {
                let res = line_search_accelerate(&poly, guess, guess_delta);
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
        const ROTATION_RADIANS: f64 = 0.9250245;
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
        const ROTATION_RADIANS: f64 = 0.9250245;
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

/// Modified Newton's method, originally conceived by [Kaj Madsen 1973](https://doi.org/10.1007/BF01933524).
///
/// This method is a much more robust than traditional naive Newton iteration.
/// It can detect when convergence slows down and increase the step size. It can
/// also detect when it gets stuck in a local minimum or maximum and unstuck
/// itself.
///
/// This implementation was based on [Henrik Vestermark 2020](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
#[allow(clippy::module_name_repetitions)]
pub struct Newton<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    history: Option<FinderHistory<T>>,
}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for Newton<T> {
    fn from_poly(poly: crate::Poly<T>) -> Self {
        Self {
            state: FinderState::new(poly),
            config: FinderConfig::new(),
            history: None,
        }
    }

    fn run(&mut self) -> std::result::Result<(), roots::Error<()>> {
        debug_assert!(self.state.poly.is_normalized());
        self.next_n_roots(self.state.poly.degree_raw().try_into().expect("overflow"))
            .map(|_| ())
            .map_err(|e| e.map_no_converge(|_| ()))
    }

    fn state_mut(&mut self) -> &mut FinderState<T> {
        &mut self.state
    }

    fn state(&self) -> &FinderState<T> {
        &self.state
    }

    fn config(&mut self) -> &mut roots::FinderConfig<T> {
        &mut self.config
    }

    fn history(&mut self) -> &mut Option<roots::FinderHistory<T>> {
        &mut self.history
    }
}

impl<T: ScalarOps + Float + RealField> IterativeRootFinder<T> for Newton<T> {
    fn next_root(&mut self) -> roots::Result<T> {
        // these are tuned by hand
        const DZ_STUCK_ROTATION_DEGREES: f64 = 90.0 / F64_PHI;
        const DZ_STUCK_SCALE: f64 = 5.0;
        const DZ_EXPLODE_ROTATION_DEGREES: f64 = 90.0 / F64_PHI;
        const DZ_EXPLODE_THRESHOLD: f64 = 5.0;

        let dz_stuck_scale = T::from_f64(DZ_STUCK_SCALE).expect("overflow");
        let dz_stuck_rotation =
            T::from_f64(DZ_STUCK_ROTATION_DEGREES.to_radians()).expect("overflow");
        let dz_stuck_factor = Complex::from_polar(dz_stuck_scale, dz_stuck_rotation);
        let dz_explode_threshold = T::from_f64(DZ_EXPLODE_THRESHOLD).expect("overflow");
        let dz_explode_rotation =
            T::from_f64(DZ_EXPLODE_ROTATION_DEGREES.to_radians()).expect("overflow");
        let min_multiplicity = T::one();

        // handle trivial cases
        let epsilon = self.config().epsilon;
        let roots = self.state_mut().poly.trivial_roots(epsilon).0;
        if !roots.is_empty() {
            return Ok(roots);
        }

        // early return for degree zero
        if self.state.poly.degree_raw() == 0 {
            return Ok(vec![]);
        }

        let p_diff = self.state.poly.clone().diff();
        let p_diff2 = p_diff.clone().diff();

        let mut guess = self
            .state
            .dirty_roots
            .pop()
            .unwrap_or_else(|| self.state.poly.initial_guess_smallest());
        let mut guess_old = guess;
        let mut guess_delta = guess;
        let mut guess_delta_old = guess;
        let mut guess_old_old = guess;

        for i in 0..self.config.max_iter {
            let px = self.state.poly.eval_point(guess);

            // stopping criterion 1: reached requested epsilon
            if px.norm() <= self.config.epsilon {
                return Ok(vec![guess]);
            }

            // stopping criterion 2: can't improve guess after at least 3 iterations
            if i >= 3 && self.stop_iteration(guess, guess_old, guess_old_old) {
                return Ok(vec![guess]);
            }

            let pdx = p_diff.eval_point(guess);
            guess_old_old = guess_old;
            guess_old = guess;
            guess_delta_old = guess_delta;

            if pdx.is_zero() {
                // if stuck in local minimum, backoff and rotate instead of converging
                guess_delta *= dz_stuck_factor;
                guess -= guess_delta;
                if let Some(history_handle) = &mut self.history {
                    let mut roots_history = history_handle.roots_history.pop().unwrap_or(vec![]);
                    roots_history.push(guess);
                    history_handle.roots_history.push(roots_history);
                }
                continue;
            } else {
                // normal Newton step
                guess_delta = px / pdx;
            }

            let guess_delta_norm = guess_delta.norm();
            let guess_delta_old_norm = guess_delta_old.norm();
            if guess_delta_norm > dz_explode_threshold * guess_delta_old_norm {
                // delta is exploding, limit it
                guess_delta *= Complex::from_polar(
                    dz_explode_threshold * guess_delta_old_norm / guess_delta_norm,
                    dz_explode_rotation,
                );
                guess -= guess_delta;
                if let Some(stats_handle) = &mut self.history {
                    let mut roots_history = stats_handle.roots_history.pop().unwrap_or(vec![]);
                    roots_history.push(guess);
                    stats_handle.roots_history.push(roots_history);
                }
                continue;
            }

            // [Schroeder 1870](https://doi.org/10.1007/BF01444024) proposes this
            // method for estimating the multiplicity efficiently
            let pddx = p_diff2.eval_point(guess);
            let pdx_2 = pdx.powu(2);
            let multiplicity = (pdx_2 / (pdx_2 - px * pddx)).norm();
            guess_delta = guess_delta.scale(Float::max(multiplicity, min_multiplicity));

            guess -= guess_delta;

            if (guess - guess_old).norm() < self.config.epsilon {
                // the solver got stuck
                guess_delta *= dz_stuck_factor;
                guess -= guess_delta;
                if let Some(stats_handle) = &mut self.history {
                    let mut roots_history = stats_handle.roots_history.pop().unwrap_or(vec![]);
                    roots_history.push(guess);
                    stats_handle.roots_history.push(roots_history);
                }
                continue;
            }

            // collect stats
            if let Some(stats_handle) = &mut self.history {
                let mut roots_history = stats_handle.roots_history.pop().unwrap_or(vec![]);
                roots_history.push(guess);
                stats_handle.roots_history.push(roots_history);
            }
        }
        self.state.dirty_roots.push(guess);
        Err(roots::Error::NoConverge(vec![guess]))
    }
}

#[cfg(test)]
mod test {
    use crate::__util::testing::check_roots;

    use crate::num::One;

    use crate::Poly64;

    use super::newton;

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
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
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
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
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
        let roots = newton(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }
}
