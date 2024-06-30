use super::{
    line_search_accelerate, line_search_decelerate, multiplicity_lagouanelle, LazyDerivatives,
};
use crate::{
    __util,
    num::{Complex, Zero},
    poly::roots,
    Poly, ScalarOps,
};
use na::RealField;

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
pub fn naive<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> roots::Result<T> {
    super::deflate(next_root, poly, epsilon, max_iter, initial_guesses)
}

/// Find a single root
///
/// # Returns
/// - vector of roots (usually 1)
/// - number of evaluations
fn next_root<T: ScalarOps + RealField>(
    poly: &Poly<T>,
    epsilon: T,
    max_iter: Option<usize>,
    initial_guess: Option<Complex<T>>,
) -> std::result::Result<(Vec<Complex<T>>, u128), roots::Error<Vec<Complex<T>>>> {
    let mut eval_counter = 0;
    let mut guess = initial_guess.unwrap_or_else(|| poly.initial_guess_smallest());
    let mut guess_old = guess;
    let mut guess_old_old = guess;
    let mut diffs = LazyDerivatives::new(poly);

    // until convergence
    for i in __util::iterator::saturating_counter() {
        let px = poly.eval(guess);

        // stopping criterion 1: converged
        if px.norm() <= epsilon {
            return Ok((vec![guess], eval_counter));
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3 && super::stopping_criterion_garwick(guess, guess_old, guess_old_old) {
            return Ok((vec![guess], eval_counter));
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        eval_counter += 1;
        let pdx = diffs.get_nth_derivative(1).eval(guess);

        // got stuck at local minimum
        if pdx.is_zero() {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

        let guess_delta = px / pdx;

        guess_old_old = guess_old;
        guess_old = guess;
        guess -= guess_delta;
    }
    unreachable!()
}

#[cfg(test)]
mod test {
    use super::naive;
    use crate::__util::testing::check_roots;
    use crate::num::One;
    use crate::Poly64;

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_5_multiplicity_3() {
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(2.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-4),
            "{roots:?}"
        );
    }

    #[test]
    fn degree_5_2_zeros() {
        let roots_expected = vec![
            complex!(0.0),
            complex!(0.0),
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }
}
