use crate::{
    __util::float::F64_PHI,
    num::{Complex, Float, Zero},
};
use na::{ComplexField, RealField};

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
/// This implementation was based on [Henrik Vestermark 2020](http://www.hvks.com/Numerical/Downloads/HVE%20Practical%20Implementation%20of%20Polynomial%20root%20finders%20vs%207.pdf).
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
        let roots = self.state_mut().poly.trivial_roots(epsilon);
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
    use crate::{poly::roots::RootFinder, Poly, __util::testing::check_roots};

    use super::Newton;

    #[test]
    fn newton_degree_3() {
        // easy: only real roots
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let p = Poly::from_roots(&roots_expected);
        let roots = Newton::from_poly(p)
            .with_epsilon(1E-14)
            .with_max_iter(100)
            .roots()
            .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn newton_degree_5() {
        // medium: some conjugate roots
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
            complex!(0.0, -1.0),
            complex!(0.0, 1.0),
        ];
        let p = Poly::from_roots(&roots_expected);
        let roots = Newton::from_poly(p)
            .with_epsilon(1E-14)
            .with_max_iter(100)
            .roots()
            .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-1));
    }
}
