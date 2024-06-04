use fastrand::Rng;
use na::RealField;
use num::{Complex, Float};

use crate::{
    poly::roots::{self, FinderConfig, FinderState, RootFinder},
    Poly, Scalar,
};

use super::IterativeRootFinder;

pub struct NewtonFinder<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
}

impl<T: Scalar + Float + RealField> NewtonFinder<T> {}

impl<T: Scalar + Float + RealField> RootFinder<T> for NewtonFinder<T> {
    fn from_poly(poly: crate::Poly<T>) -> Self {
        Self {
            state: FinderState::new(poly.clone()),
            config: FinderConfig::new(),
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
}

impl<T: Scalar + Float + RealField> IterativeRootFinder<T> for NewtonFinder<T> {
    fn next_root(&mut self) -> roots::Result<T> {
        let mut guess = self
            .state
            .dirty_roots
            .pop()
            .unwrap_or(self.state.poly.initial_guess_smallest());
        let p_diff = self.state.poly.clone().diff();
        for _ in 0..self.config.max_iter {
            let px = self.state.poly.eval_point(guess);
            if px.norm() <= self.config.epsilon {
                return Ok(vec![guess]);
            }
            let pdx = p_diff.eval_point(guess);
            guess = guess - px / pdx;
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
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
            complex!(0.0, -1.0),
            complex!(0.0, 1.0),
        ];
        let p = Poly::from_roots(&roots_expected);
        let mut roots = NewtonFinder::from_poly(p)
            .with_epsilon(1E-14)
            .with_max_iter(100)
            .roots()
            .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-1));
    }
}
