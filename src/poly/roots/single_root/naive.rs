use super::LazyDerivatives;
use crate::{
    num::{Complex, Zero},
    poly::roots,
    roots::initial_guess::initial_guess_smallest,
    util::{
        self,
        doc_macros::{errors_no_converge, panic_t_from_f64},
    },
    Poly, RealScalar,
};

/// Find a single root using naive Newton's method
///
/// # Returns
/// - vector of roots (usually 1)
/// - number of evaluations
///
/// # Errors
#[doc = errors_no_converge!()]
///
/// # Panics
#[doc = panic_t_from_f64!()]
#[allow(clippy::similar_names)]
#[allow(clippy::items_after_statements)]
#[allow(clippy::type_complexity)]
pub fn naive<T: RealScalar>(
    poly: &Poly<T>,
    epsilon: T,
    max_iter: Option<usize>,
    initial_guess: Option<Complex<T>>,
) -> std::result::Result<(Vec<Complex<T>>, u128), roots::Error<Vec<Complex<T>>>> {
    let mut guess = initial_guess.unwrap_or_else(|| initial_guess_smallest(poly));
    let mut guess_old = guess;
    let mut guess_old_old = guess;
    let mut diffs = LazyDerivatives::new(poly);

    // until convergence
    for i in util::iterator::saturating_counter() {
        let px = poly.eval(guess);

        // stopping criterion 1: converged
        if px.norm() <= epsilon {
            return Ok((vec![guess], i as u128));
        }

        // stopping criterion 2: no improvement predicted due to numeric precision
        if i > 3 && super::stopping_criterion_garwick(guess, guess_old, guess_old_old) {
            return Ok((vec![guess], i as u128));
        }

        // max iter exceeded
        if max_iter.is_some_and(|max| i >= max) {
            return Err(roots::Error::NoConverge(vec![guess]));
        }

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
    use crate::{num::One, roots::naive_deflate, util::__testing::check_roots, Poly64};

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
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
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
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
        let roots = naive_deflate(&mut p, Some(1E-14), Some(100), &[]).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }
}
