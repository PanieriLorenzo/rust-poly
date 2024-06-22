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

#[cfg(test)]
mod test {
    use crate::{__util::testing::check_roots, roots::History};

    use crate::num::One;

    use crate::Poly64;

    use super::naive;

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn test_degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn test_degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn test_degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn test_degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn test_degree_5_2_zeros() {
        let roots_expected = vec![
            complex!(0.0),
            complex!(0.0),
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::naive(&mut p, Some(1E-14), Some(100)).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }
}
