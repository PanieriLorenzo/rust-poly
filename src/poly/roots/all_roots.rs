//! Root finders that find all roots at once.

pub mod aberth_ehrlich;
pub use aberth_ehrlich::aberth_ehrlich;
use na::{Complex, RealField};

use crate::{Poly, ScalarOps};

use super::{halley, naive, newton, single_root::NextRootFun};

/// Combinator that makes an all-roots root finder from a single-root root finder.
pub fn deflate<T: ScalarOps + RealField>(
    next_root_fun: NextRootFun<T>,
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
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
        let (next_roots, num_evals) = next_root_fun(poly, epsilon, max_iter, next_guess)?;
        let root = next_roots[0];
        eval_counter += num_evals;
        roots.push(root);
        // TODO: deflate_composite should borrow instead
        *poly = poly.clone().deflate_composite(root);
    }
}

/// Find all roots of a polynomial using a modified Halley's method.
///
/// This implementation is based on [Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
#[inline]
pub fn halley_deflate<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    deflate(halley, poly, epsilon, max_iter, initial_guesses)
}

/// Find all roots of a polynomial using a naive Newton-Raphson approach.
///
/// You are likely looking for [`roots::newton`] instead, which is a more
/// robust implementation of the same algorithm.
///
/// This algorithm is not very stable, and will often fail to find fairly ordinary
/// roots. So only use this approach for specific purposes, namely:
/// - You know its not going to get stuck and you want a slight performance improvement
/// - You are benchmarking your custom root finder against a classical Newton-Raphson approach.
#[inline]
pub fn naive_deflate<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    deflate(naive, poly, epsilon, max_iter, initial_guesses)
}

/// Modified Newton's method, originally conceived by [Kaj Madsen 1973](https://doi.org/10.1007/BF01933524).
///
/// This method is a much more robust than traditional naive Newton iteration.
/// It can detect when convergence slows down and increase the step size. It can
/// also detect when it gets stuck in a local minimum or maximum and unstuck
/// itself.
///
/// This implementation was based on [Henrik Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
#[inline]
pub fn newton_deflate<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    deflate(newton, poly, epsilon, max_iter, initial_guesses)
}
