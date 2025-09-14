use std::ops::ControlFlow::{self, Break, Continue};

use crate::{
    aliases::{C, R},
    base::{BasePoly, UnivariateMarker},
    scalar_traits::{BasicScalar, NonIntegerScalar},
    storage_traits::{BasicStorage, MutStorage, OwnedStorage},
};

use num_traits::{One, Zero};

use super::errors::{BreakReason, ContinueReason};

pub fn aberth_ehrlich<S: MutStorage>(
    poly: &mut BasePoly<S, UnivariateMarker>,
    guesses: &mut Vec<C<S>>,
    min_iter: usize,
    max_iter: Option<usize>,
    stop_epsilon: R<S>,
) -> ControlFlow<BreakReason, ContinueReason>
where
    S::T: NonIntegerScalar,
{
    debug_assert!(poly.is_normalized_inner());
    debug_assert!(guesses.len() == poly.degree_raw_inner());
    debug_assert!(poly.degree_raw_inner() > 2);
    let d = poly.degree_raw_inner();
    let n = poly.len_inner();
    for i in 0..d {
        for j in (i + 1)..d {
            if !((guesses[i].clone() - guesses[j].clone()).taxicab_norm() > R::<S>::zero()) {
                return Break(BreakReason::DuplicateGuesses);
            }
        }
    }

    poly.make_monic_inner();
    let mut alphas = vec![C::<S>::zero(); d];
    let mut betas = vec![C::<S>::zero(); d];
    let mut i = 0usize;
    loop {
        if max_iter.is_some_and(|max| i >= max) {
            return Continue(ContinueReason::MaxIter);
        }

        // update alphas
        let mut p_diff = poly.into_owned_poly_inner();
        p_diff.diff_inner();
        poly.eval_multiple_complex_inner(guesses, &mut alphas);
        for (y, x) in alphas.iter_mut().zip(guesses.iter()) {
            *y = y.clone() / p_diff.eval_complex_inner(x.clone());
        }

        // update betas
        betas.fill(C::<S>::zero());
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                betas[i] =
                    betas[i].clone() + C::<S>::one() / (guesses[i].clone() - guesses[j].clone());
            }
        }

        // alphas become deltas in-place
        for (a, b) in alphas.iter_mut().zip(betas.iter()) {
            *a = a.clone() / (C::<S>::one() - a.clone() * b.clone());
        }
        let deltas = &mut alphas;

        for (y, d) in guesses.iter_mut().zip(deltas.iter()) {
            *y = y.clone() - d.clone();
        }

        // stopping criteria
        if i >= min_iter && deltas.iter().all(|d| d.taxicab_norm() <= stop_epsilon) {
            return Continue(ContinueReason::Ok);
        }
        i = i.saturating_add(1);
    }
}
