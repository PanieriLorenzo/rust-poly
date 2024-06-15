use crate::{
    __util::float::F64_PHI,
    num::{Complex, Float, Zero},
    roots::Newton,
};
use na::{ComplexField, RealField};

use crate::{
    poly::roots::{self, FinderConfig, FinderState, RootFinder},
    roots::FinderHistory,
    Scalar, ScalarOps,
};

use super::IterativeRootFinder;

/// Naive Newton's method.
///
/// It's strongly recommended to use a different root finder, as naive Newton
/// iteration is not robust at all. Use this only if dealing with small well-behaved
/// polynomials were the simplicity gives a performance advantage, or as a baseline
/// for benchmarking purposes.
///
/// Unlike [`roots::Madsen`] this has none of the tweaks that improve convergence
/// on pathological roots. This will often get stuck on any slightly challenging
/// root.
#[allow(clippy::module_name_repetitions)]
pub struct Naive<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    history: Option<FinderHistory<T>>,
}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for Naive<T> {
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

    fn config(&mut self) -> &mut FinderConfig<T> {
        &mut self.config
    }

    fn history(&mut self) -> &mut Option<FinderHistory<T>> {
        &mut self.history
    }
}

impl<T: ScalarOps + Float + RealField> IterativeRootFinder<T> for Naive<T> {
    fn next_root(&mut self) -> roots::Result<T> {
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

        // TODO: move this to the next_n_roots method
        // TODO: it should retry if more initial guesses are available
        let mut guess = self
            .state
            .dirty_roots
            .pop()
            .unwrap_or_else(|| self.state.poly.initial_guess_smallest());
        let mut guess_old = guess;
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

            if pdx.is_zero() {
                // got stuck at local minimum
                break;
            }

            let guess_delta = px / pdx;

            guess_old_old = guess_old;
            guess_old = guess;
            guess -= guess_delta;

            // collect history
            if let Some(history_handle) = &mut self.history {
                let mut roots_history = history_handle.roots_history.pop().unwrap_or(vec![]);
                roots_history.push(guess);
                history_handle.roots_history.push(roots_history);
            }
        }
        self.state.dirty_roots.push(guess);
        Err(roots::Error::NoConverge(vec![guess]))
    }
}
