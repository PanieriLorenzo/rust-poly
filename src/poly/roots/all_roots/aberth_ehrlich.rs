use na::RealField;

use crate::{
    num::{Complex, One, Zero},
    roots::{self, initial_guess::initial_guesses_circle},
    util::{
        self,
        doc_macros::{errors_no_converge, panic_t_from_f64},
    },
    Poly, RealScalar,
};

/// TODO: document this
///
/// # Errors
#[doc = errors_no_converge!()]
///
/// # Panics
/// If the provided guesses are not unique, i.e. if two or more are the same.
///
#[doc = panic_t_from_f64!()]
pub fn aberth_ehrlich<T: RealScalar + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> roots::Result<T> {
    debug_assert!(poly.is_normalized());

    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let trivial_roots = poly.trivial_roots(epsilon).0;

    debug_assert!(poly.is_normalized());
    if poly.degree_raw() == 0 {
        return Ok(trivial_roots);
    }

    poly.make_monic();

    // fill remaining initial guesses
    // TODO: this should be factored out to a separate function
    let initial_guesses = {
        let mut complete_initial_guesses = Vec::with_capacity(poly.degree_raw());
        for z in initial_guesses {
            complete_initial_guesses.push(*z);
        }
        let remaining_guesses_delta = poly.degree_raw() - complete_initial_guesses.len();
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

    let n = poly.degree_raw();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            assert!(
                (initial_guesses[i] - initial_guesses[j]).norm() > T::zero(),
                "initial guesses must be distinct"
            );
        }
    }

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
            *a /= Complex::<T>::one() - *a * b;
        }
        let deltas_buff = &mut alphas_buff;

        for (y, d) in points.iter_mut().zip(deltas_buff.iter()) {
            *y -= d;
        }

        // stopping criteria
        if deltas_buff.iter().all(|d| d.norm() <= epsilon) {
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
        *y /= p_diff.eval(*x);
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
            out[i] += Complex::<T>::one() / (points[i] - points[j]);
        }
    }
}

#[cfg(test)]
mod test {
    use num::{complex::Complex64, Zero};

    use super::aberth_ehrlich;
    use crate::{
        roots::initial_guess::{initial_guesses_circle, initial_guesses_random},
        util::__testing::check_roots,
    };

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
        initial_guesses_random(p.clone(), 1, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-14), Some(100), &guesses).unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let mut guesses = [Complex64::zero(); 3];
        initial_guesses_random(p.clone(), 1, &mut guesses);
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
        initial_guesses_random(p.clone(), 1, &mut guesses);
        let roots = aberth_ehrlich(&mut p, Some(1E-5), Some(100), &guesses).unwrap();
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
        let roots = aberth_ehrlich(&mut p, Some(1E-5), Some(100), &guesses).unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-4),
            "{roots:?}"
        );
    }
}
