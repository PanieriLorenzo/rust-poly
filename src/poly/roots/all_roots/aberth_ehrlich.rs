use crate::{
    num::{Complex, One, Zero},
    roots::{self},
    util::{
        self,
        doc_macros::{errors_no_converge, panic_t_from_f64},
    },
    Poly, RealScalar,
};

/// Find all roots using Aberth Ehrlich method.
///
/// # Caveats
/// This method performs poorly around zero roots, so you should remove them
/// first (zero roots are trivial to factor out).
///
/// # Errors
#[doc = errors_no_converge!()]
///
/// # Panics
/// If the provided guesses are not unique, i.e. if two or more are the same.
///
/// If the number of provided guesses is wrong, there must be exactly one guess
/// per root, i.e. as many guesses as the degree of the polynomial.
///
#[doc = panic_t_from_f64!()]
pub fn aberth_ehrlich<T: RealScalar>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> roots::Result<T> {
    debug_assert!(poly.is_normalized());
    debug_assert!(initial_guesses.len() == poly.degree_raw());

    let n = poly.degree_raw();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            assert!(
                (initial_guesses[i].clone() - initial_guesses[j].clone()).norm_sqr() > T::zero(),
                "initial guesses must be distinct"
            );
        }
    }

    let epsilon = epsilon.unwrap_or(T::tiny_safe());

    if poly.degree_raw() == 0 {
        return Ok(vec![]);
    }

    poly.make_monic();

    let mut points = Vec::from(&initial_guesses[..n]);
    let mut alphas_buff = vec![Complex::<T>::zero(); n];
    let mut betas_buff = vec![Complex::<T>::zero(); n];

    for i in util::iterator::saturating_counter() {
        if max_iter.is_some_and(|max| i > max) {
            return Err(roots::Error::NoConverge(points));
        }

        alphas(poly, &points, &mut alphas_buff);
        betas(&points, &mut betas_buff);

        // alphas become deltas in-place
        for (a, b) in alphas_buff.iter_mut().zip(betas_buff.iter()) {
            *a /= Complex::<T>::one() - a.clone() * b;
        }
        let deltas_buff = &mut alphas_buff;

        for (y, d) in points.iter_mut().zip(deltas_buff.iter()) {
            *y -= d;
        }

        log::trace!("{points:?}");

        // stopping criteria
        if deltas_buff.iter().all(|d| d.norm_sqr() <= epsilon) {
            return Ok(points);
        }
    }
    unreachable!();
}

/// Alpha coefficients of the Aberth-Ehrlich method
///
/// Needs `points.len() == out.len()`.
fn alphas<T: RealScalar>(poly: &Poly<T>, points: &[Complex<T>], out: &mut [Complex<T>]) {
    debug_assert_eq!(points.len(), out.len());

    let p_diff = poly.clone().diff();

    // TODO: division by zero
    poly.eval_multiple(points, out);
    for (y, x) in out.iter_mut().zip(points) {
        *y /= p_diff.eval(x.clone());
    }
}

/// Beta coefficients of the Aberth-Ehrlich method
///
/// Needs `points.len() == out.len()`.
fn betas<T: RealScalar>(points: &[Complex<T>], out: &mut [Complex<T>]) {
    debug_assert_eq!(points.len(), out.len());

    let n = points.len();
    out.fill(Complex::zero());
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            out[i] += Complex::<T>::one() / (points[i].clone() - points[j].clone());
        }
    }
}

#[cfg(test)]
mod test {
    use num::{complex::Complex64, Zero};

    use super::aberth_ehrlich;
    use crate::{roots::initial_guess::initial_guesses_circle, util::__testing::check_roots};

    // #[test]
    // pub fn degree_0() {
    //     let mut p = Poly64::one();
    //     let mut guesses = [Complex64::zero(); 1];
    //     initial_guesses_random(&p, 1, &mut guesses);
    //     let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &mut guesses).unwrap();
    //     assert!(roots.is_empty());
    //     assert!(p.is_one());
    // }

    // #[test]
    // fn degree_1() {
    //     let roots_expected = vec![complex!(1.0)];
    //     let mut p = crate::Poly::from_roots(&roots_expected);
    //     let mut guesses = [Complex64::zero(); 1];
    //     initial_guesses_random(&p, 1, &mut guesses);
    //     let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &mut guesses).unwrap();
    //     assert!(check_roots(roots, roots_expected, 1E-12));
    // }

    // #[test]
    // fn degree_2() {
    //     let roots_expected = vec![complex!(1.0), complex!(2.0)];
    //     let mut p = crate::Poly::from_roots(&roots_expected);
    //     let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &[]).unwrap();
    //     assert!(check_roots(roots, roots_expected, 1E-12));
    // }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let mut guesses = [Complex64::zero(); 3];
        initial_guesses_circle(&p, 0.5, 1, 0.5, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &guesses).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let mut guesses = [Complex64::zero(); 3];
        initial_guesses_circle(&p, 0.5, 1, 0.5, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &guesses).unwrap();
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
        let mut guesses = [Complex64::zero(); 5];
        initial_guesses_circle(&p, 0.5, 1, 0.5, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-8), Some(100), &guesses).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-4),
            "{roots:?}"
        );
    }

    #[test]
    fn degree_10_multiplicity_3_2() {
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(2.0),
            complex!(2.0),
            complex!(3.0),
            complex!(0.0, 1.0),
            complex!(0.0, 2.0),
            complex!(0.0, 2.0),
            complex!(0.0, 3.0),
            complex!(-1.0, -1.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let mut guesses = [Complex64::zero(); 10];
        initial_guesses_circle(&p, 0.5, 1, 0.5, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-8), Some(100), &guesses).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-4),
            "{roots:?}"
        );
    }
}
