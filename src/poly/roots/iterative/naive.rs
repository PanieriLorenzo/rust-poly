use crate::{
    __util::{self, float::F64_PHI},
    num::{Complex, Float, Zero},
    roots::{History, Newton},
};
use na::{ComplexField, RealField};

use crate::{
    poly::roots::{self, FinderConfig, FinderState, RootFinder},
    roots::FinderHistory,
    Poly, Scalar, ScalarOps,
};

/// A builder for the [`naive`] root finder. See [`naive`] for more info.
pub struct NaiveA;

/// Find all roots of a polynomial using a naive Newton-Raphson approach.
///
/// You are likely looking for [`roots::newton`] instead, which is a more
/// robust implementation of the same algorithm.
///
/// This algorithm is not very stable, and will often fail to find fairly ordinary
/// roots. So only use this approach for specific purposes, namely:
/// - You know its not going to get stuck and you want a slight performance improvement
/// - You are benchmarking your custom root finder against a classical Newton-Raphson approach.
pub fn naive<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
) -> roots::Result<T> {
    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let mut roots = vec![];

    // until we've found all roots
    loop {
        let trivial_roots = poly.trivial_roots(epsilon);
        roots.extend(trivial_roots.iter());

        debug_assert!(poly.is_normalized());
        if poly.degree_raw() == 0 {
            return Ok(roots);
        }

        let root = next_root(poly, epsilon, max_iter)?[0];
        roots.push(root);
        // TODO: deflate_composite should borrow instead
        *poly = poly.clone().deflate_composite(root);
    }
}

fn next_root<T: ScalarOps + RealField>(
    poly: &Poly<T>,
    epsilon: T,
    max_iter: Option<usize>,
) -> roots::Result<T> {
    let mut guess = poly.initial_guess_smallest();
    let mut guess_old = guess;
    let mut guess_old_old = guess;
    let p_diff = poly.clone().diff();

    // until convergence
    for i in __util::iterator::saturating_counter() {
        let px = poly.eval_point(guess);

        // stopping criterion 1: converged
        if px.norm() <= epsilon {
            return Ok(vec![guess]);
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3 && super::stopping_criterion_garwick(guess, guess_old, guess_old_old) {
            return Ok(vec![guess]);
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        let pdx = p_diff.eval_point(guess);

        // got stuck at local minimum
        if pdx.is_zero() {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        let guess_delta = px / pdx;

        guess_old_old = guess_old;
        guess_old = guess;
        guess = guess - guess_delta;
    }
    unreachable!()
}

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
pub struct NaiveOld<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    history: Option<FinderHistory<T>>,
}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for NaiveOld<T> {
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

impl<T: ScalarOps + Float + RealField> IterativeRootFinder<T> for NaiveOld<T> {
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

#[cfg(test)]
mod test {
    use crate::{__util::testing::check_roots, roots::History};

    #[test]
    fn test_naive() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
        panic!();
    }
}
