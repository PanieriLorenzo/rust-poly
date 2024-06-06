use na::RealField;
use num::Float;

use crate::{
    poly::roots::{self, FinderConfig, FinderState, RootFinder},
    roots::FinderStatistics,
    Poly, Scalar, ScalarOps,
};

use super::IterativeRootFinder;

pub struct NewtonFinder<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    statistics: Option<FinderStatistics<T>>,
    poly_original: Poly<T>,
}

impl<T: ScalarOps + Float + RealField> NewtonFinder<T> {}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for NewtonFinder<T> {
    fn from_poly(poly: crate::Poly<T>) -> Self {
        Self {
            state: FinderState::new(poly.clone()),
            config: FinderConfig::new(),
            statistics: None,
            poly_original: poly,
        }
    }

    fn roots(&mut self) -> roots::Result<T> {
        debug_assert!(self.state.poly.is_normalized());
        self.next_n_roots(self.state.poly.degree_raw() as usize)
    }

    fn state(&mut self) -> &mut FinderState<T> {
        &mut self.state
    }

    fn config(&mut self) -> &mut roots::FinderConfig<T> {
        &mut self.config
    }

    fn statistics(&mut self) -> &mut Option<roots::FinderStatistics<T>> {
        &mut self.statistics
    }
}

impl<T: ScalarOps + Float + RealField> IterativeRootFinder<T> for NewtonFinder<T> {
    fn next_root(&mut self) -> roots::Result<T> {
        //self.state.poly.make_monic();
        let mut guess = self
            .state
            .dirty_roots
            .pop()
            .unwrap_or(self.state.poly.initial_guess_smallest());
        let mut old_guess = guess;
        let p_diff = self.state.poly.clone().diff();
        for _ in 0..self.config.max_iter {
            let px = self.state.poly.eval_point(guess);
            if px.norm() <= self.config.epsilon {
                return Ok(vec![guess]);
            }
            let pdx = p_diff.eval_point(guess);
            guess -= px / pdx;
            if guess == old_guess {
                // the solver got stuck, unstuck it early
                break;
            }
            old_guess = guess;
            // collect stats
            if let Some(stats_handle) = &mut self.statistics {
                let mut roots_history = stats_handle.roots_history.pop().unwrap_or(vec![]);
                let mut roots_err_sq_history =
                    stats_handle.roots_err_sq_history.pop().unwrap_or(vec![]);
                roots_history.push(guess);
                let err = self.poly_original.eval_point(guess).norm_sqr();
                roots_err_sq_history.push(err);
                stats_handle.roots_history.push(roots_history);
                stats_handle.roots_err_sq_history.push(roots_err_sq_history);
            }
        }
        self.state.dirty_roots.push(guess);
        Err(roots::Error::NoConverge(vec![guess]))
    }
}

#[cfg(test)]
mod test {
    use crate::{poly::roots::RootFinder, Poly, __util::testing::check_roots};

    use super::NewtonFinder;

    #[test]
    fn newton_degree_3() {
        // easy: only real roots
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let p = Poly::from_roots(&roots_expected);
        let roots = NewtonFinder::from_poly(p)
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
        let roots = NewtonFinder::from_poly(p)
            .with_epsilon(1E-14)
            .with_max_iter(100)
            .roots()
            .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-1));
    }

    #[test]
    fn newton_degree_15() {
        // hard: multiplicity and non-conjugate roots
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
            complex!(4.0),
            complex!(5.0),
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
            complex!(-1.0),
            complex!(-2.0),
            complex!(-3.0),
            complex!(0.0, -1.0),
            complex!(0.0, 1.0),
            complex!(0.5, -2.0),
            complex!(1.0, -3.0),
        ];
        let p = Poly::from_roots(&roots_expected);
        let roots = NewtonFinder::from_poly(p)
            .with_epsilon(1E-8)
            .with_max_iter(100)
            .roots()
            .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-1));
    }
}
