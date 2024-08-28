//! Methods that find multiple roots at the same time (but not all roots)

use super::{halley, naive, newton, single_root::NextRootFun};
use crate::{
    num::Complex,
    roots::initial_guesses_circle,
    util::doc_macros::{errors_no_converge, panic_t_from_f64},
    Poly, RealScalar,
};
use num::Zero;

/// Run multiple single-root root finders in parallel.
///
/// # Errors
/// Passes through whichever errors the passed closure `next_root_fun` returned.
///
/// # Panics
#[doc = panic_t_from_f64!()]
pub fn parallel<T: RealScalar>(
    next_root_fun: NextRootFun<T>,
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    poly.make_monic();

    let mut roots = vec![];
    for z in initial_guesses {
        roots.extend(
            next_root_fun(
                poly,
                epsilon.clone().unwrap_or(T::zero()),
                max_iter,
                Some(z.clone()),
            )?
            .0,
        );
    }
    Ok(roots)
}

/// Use Naive Newton's method in parallel for multiple initial guesses.
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn naive_parallel<T: RealScalar>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(naive, poly, epsilon, max_iter, initial_guesses)
}

/// Use Newton's method in parallel for multiple initial guesses.
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn newton_parallel<T: RealScalar>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(newton, poly, epsilon, max_iter, initial_guesses)
}

/// Use Halley's method in parallel for multiple initial guesses.
///
/// # Errors
#[doc = errors_no_converge!()]
#[inline]
pub fn halley_parallel<T: RealScalar>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(halley, poly, epsilon, max_iter, initial_guesses)
}
