//! Methods that find multiple roots at the same time (but not all roots)

use super::{halley, naive, newton, single_root::NextRootFun};
use crate::{num::Complex, roots::initial_guesses_circle, Poly, ScalarOps};
use na::RealField;
use num::Zero;

/// Run multiple single-root root finders in parallel.
pub fn parallel<T: ScalarOps + RealField>(
    next_root_fun: NextRootFun<T>,
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    num: usize,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    poly.make_monic();

    let initial_guesses = {
        let mut complete_initial_guesses = Vec::with_capacity(num);
        for z in initial_guesses {
            complete_initial_guesses.push(*z);
        }
        let remaining_guesses_delta = num - complete_initial_guesses.len();
        let mut remaining_guesses = vec![Complex::zero(); remaining_guesses_delta];
        initial_guesses_circle(
            poly,
            T::from_f64(0.5).expect("overflow"),
            1,
            T::from_f64(0.5).expect("overflow"),
            &mut remaining_guesses,
        );
        for z in remaining_guesses.drain(..) {
            complete_initial_guesses.push(z);
        }
        complete_initial_guesses
    };

    let mut roots = vec![];
    for z in initial_guesses {
        roots.extend(next_root_fun(poly, epsilon.unwrap_or(T::zero()), max_iter, Some(z))?.0);
    }
    Ok(roots)
}

#[inline]
pub fn naive_parallel<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    num: usize,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(naive, poly, epsilon, max_iter, num, initial_guesses)
}

#[inline]
pub fn newton_parallel<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    num: usize,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(newton, poly, epsilon, max_iter, num, initial_guesses)
}

#[inline]
pub fn halley_parallel<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    num: usize,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    parallel(halley, poly, epsilon, max_iter, num, initial_guesses)
}
