//! Root finders that find all roots at once.

pub mod aberth_ehrlich;
pub use aberth_ehrlich::aberth_ehrlich;
use num::Complex;

use crate::{util::doc_macros::errors_no_converge, Poly, RealScalar};

use super::{halley, naive, newton, single_root::NextRootFun};

/// Combinator that makes an all-roots root finder from a single-root root finder.
///
/// # Errors
/// Passes through any errors returned by the passed closure `next_root_fun`.
pub fn deflate<T: RealScalar>(
    next_root_fun: NextRootFun<T>,
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    debug_assert!(poly.is_normalized());
    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let mut roots = vec![];
    let mut initial_guesses = initial_guesses.iter().cloned();

    // until we've found all roots
    loop {
        let trivial_roots = {
            let mut roots = poly.zero_roots(&epsilon);
            match poly.degree_raw() {
                1 => roots.extend(poly.linear_roots()),
                2 => roots.extend(poly.quadratic_roots()),
                _ => {}
            }
            roots
        };
        roots.extend(trivial_roots.iter().cloned());

        debug_assert!(poly.is_normalized());
        if poly.degree_raw() == 0 {
            return Ok(roots);
        }

        let next_guess = initial_guesses.next();
        let (next_roots, _) = next_root_fun(poly, epsilon.clone(), max_iter, next_guess)?;
        let root = next_roots[0].clone();
        roots.push(root.clone());
        poly.deflate_composite(&root);
    }
}

/// Find all roots of a polynomial using a modified Halley's method.
///
/// This implementation is based on [Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728).
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn halley_deflate<T: RealScalar>(
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
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn naive_deflate<T: RealScalar>(
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
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn newton_deflate<T: RealScalar>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    deflate(newton, poly, epsilon, max_iter, initial_guesses)
}
