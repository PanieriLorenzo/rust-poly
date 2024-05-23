use na::{Complex, ComplexField, Normed, RealField, Scalar};
use num::{traits::float::FloatCore, Float, FromPrimitive, Num, One, Zero};

use crate::{
    Poly, ScalarOps,
    __util::{
        casting::usize_to_scalar,
        complex::{c_min, c_neg, complex_sort_mut},
    },
};

/// Polynomial root-finding algorithms
#[non_exhaustive]
pub enum OneRootAlgorithms {
    Newton,
    Halley,
    Laguerre,
}

#[non_exhaustive]
pub enum AllRootsAlgorithms {
    Schur,
}

// private
impl<T: Scalar + RealField + Float> Poly<T> {
    /// Ref: https://doi.org/10.1007/BF01933524
    fn initial_guess_smallest(&self) -> Complex<T> {
        debug_assert!(self.is_normalized());
        debug_assert!(self.len_raw() >= 2);

        let small = Float::recip(usize_to_scalar::<T>(1_000));
        let p_diff = self.clone().diff();
        let mut pz = self.eval_point(Complex::zero());
        let mut pdz = p_diff.eval_point(Complex::zero());

        // avoid divide by zero
        if pdz.norm() < small {
            pz = pz + small;
            pdz = pdz + small;
        }

        let theta = (c_neg(pz) / pdz).arg();
        let mut iter_coeffs = self.0.iter();
        let a0 = iter_coeffs.next().expect("infallible");

        let mut guess = iter_coeffs
            .zip(1..)
            .map(|(ak, k)| {
                Complex::i()
                    .scale(theta)
                    .exp()
                    .scale((a0 / ak).norm())
                    .powf(T::one() / usize_to_scalar(k))
            })
            .reduce(c_min)
            .expect("infallible")
            .scale(Float::recip(usize_to_scalar::<T>(2)));

        if guess.im.is_zero() {
            // add a small constant because some methods can't converge to
            // complex roots if the initial guess is real
            guess = guess + Complex::i().scale(Float::recip(usize_to_scalar::<T>(1_000)));
        }
        guess
    }

    fn one_root_newton(
        &self,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
    ) -> Result<Complex<T>, Complex<T>> {
        let p_diff = self.clone().diff();
        let mut x = initial_guess.unwrap_or(self.initial_guess_smallest());
        for _ in 0..max_iter {
            let px = self.eval_point(x);
            if px.norm() <= epsilon {
                return Ok(x);
            }
            let pdx = p_diff.eval_point(x);
            x = x - px / pdx;
        }
        Err(x)
    }

    fn one_root_halley(
        &self,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
    ) -> Result<Complex<T>, Complex<T>> {
        let p_diff = self.clone().diff();
        let p_diff2 = p_diff.clone().diff();

        let mut x = initial_guess.unwrap_or(self.initial_guess_smallest());
        for _ in 0..max_iter {
            let px = self.eval_point(x);
            if px.norm() <= epsilon {
                return Ok(x);
            }
            let pdx = p_diff.eval_point(x);
            let pddx = p_diff2.eval_point(x);
            let two = Complex::from_u32(2).expect("infallible");
            x = x - (px * pdx * two) / (pdx.powu(2) * two - px * pddx);
        }
        Err(x)
    }

    fn one_root_laguerre(
        &self,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
    ) -> Result<Complex<T>, Complex<T>> {
        let p_diff = self.clone().diff();
        let p_diff2 = p_diff.clone().diff();
        let degree = Complex::from_i32(self.degree_raw()).expect("degree too big");
        let mut x = initial_guess.unwrap_or(self.initial_guess_smallest());
        for _ in 0..max_iter {
            let px = self.eval_point(x);
            if px.norm() <= epsilon {
                return Ok(x);
            }
            let pdx = p_diff.eval_point(x);
            let pddx = p_diff2.eval_point(x);
            let g = pdx / px;
            let g2 = g.powu(2);
            let h = g2 - pddx / px;
            let degree_m_1 = <Complex<T> as std::ops::Sub>::sub(degree, Complex::one());
            let denom_pm = (degree_m_1 * (degree * h - g2)).sqrt();
            let denom1 = g + denom_pm;
            let denom2 = g - denom_pm;
            let a = if denom1.norm() > denom2.norm() {
                denom1
            } else {
                denom2
            };
            x = x - a;
        }
        Err(x)
    }

    fn one_root_ostrowski(&self, epsilon: T, max_iter: usize) -> Result<Complex<T>, Complex<T>> {
        todo!()
    }

    fn one_root_ostrowski_sq(&self, epsilon: T, max_iter: usize) -> Result<Complex<T>, Complex<T>> {
        todo!()
    }

    fn one_root_householder3(&self, epsilon: T, max_iter: usize) -> Result<Complex<T>, Complex<T>> {
        todo!()
    }

    fn roots_durand_kerner(
        &self,
        epsilon: T,
        max_iter: usize,
    ) -> Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        todo!()
    }

    fn roots_jenkins_traub(
        &self,
        epsilon: T,
        max_iter: usize,
    ) -> Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        todo!()
    }

    fn roots_schur(&self, epsilon: T, max_iter: usize) -> Result<Vec<Complex<T>>, ()> {
        // NOTE: this algorithm has a pathological case when the polynomial
        //       has only odd-degree or even-degree terms, where it does not
        //       converge. Should have some sort of pattern defeating pre-processing step
        let comp = self.companion();

        Ok(comp
            .try_schur(epsilon, max_iter)
            .ok_or(())?
            .eigenvalues()
            .expect("never fails on complex matrices")
            .as_slice()
            .into())
    }
}

impl<T: Scalar + Float + RealField> Poly<T> {
    /// Find only some of the roots of the polynomial.
    ///
    /// Note that for large `n`, using [`Poly::try_roots`] is probably faster.
    ///
    /// It utilizes an iterative method, so the precision gets progressively
    /// worse the more roots are found. For small `n` this is negligible.
    ///
    /// `Err` result contains the roots it was able to find, even if they are
    /// fewer than requested.
    ///
    /// Use [`Poly::try_n_roots_algo`] to specify which algorithm to use, if
    /// you already know which one will perform best.
    pub fn try_n_roots(
        &self,
        n: usize,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
        algorithm: Option<OneRootAlgorithms>,
    ) -> Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        debug_assert!(self.is_normalized());
        assert!(
            n as i32 <= self.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

        let algorithm = algorithm.unwrap_or(OneRootAlgorithms::Newton);

        let mut roots = vec![];
        let mut this = self.clone();
        for i in 0..n {
            let r = match algorithm {
                OneRootAlgorithms::Newton => this
                    .clone()
                    .one_root_newton(initial_guess, epsilon, max_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::Halley => this
                    .clone()
                    .one_root_halley(initial_guess, epsilon, max_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::Laguerre => this
                    .clone()
                    .one_root_laguerre(initial_guess, epsilon, max_iter)
                    .map_err(|_| roots.clone())?,
            };
            roots.push(r.clone());
            if i < (n - 1) {
                this = this / Poly::from_roots(&[r.clone()]);
            }
        }
        Ok(roots)
    }

    pub fn try_roots(
        &self,
        epsilon: T,
        max_iter: usize,
        max_tries: usize,
        max_recovery_iter: Option<usize>,
        algorithm: Option<AllRootsAlgorithms>,
        recovery_algorithm: Option<OneRootAlgorithms>,
    ) -> Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        debug_assert!(self.is_normalized());

        if self.len_raw() < 2 {
            return Ok(vec![]);
        }

        let max_recovery_iter = max_recovery_iter.unwrap_or(max_iter);
        let algorithm = algorithm.unwrap_or(AllRootsAlgorithms::Schur);
        let recovery_algorithm = recovery_algorithm.unwrap_or(OneRootAlgorithms::Newton);

        let all_roots = match algorithm {
            AllRootsAlgorithms::Schur => |this: Poly<T>| this.roots_schur(epsilon, max_iter),
            _ => unimplemented!(),
        };

        let mut roots = vec![];
        let mut this = self.clone();
        for _ in 0..max_tries {
            // TODO: automatically skip to recovery if odd-degree only or even-degree
            //       only as these are pathological cases for eigenvalue root finders
            let maybe_roots = all_roots(this.clone());
            if let Ok(found_roots) = maybe_roots {
                roots.extend(found_roots);
                // TODO: sort
                return Ok(roots);
            }

            // uses one iteration of single-root algorithm for recovery when
            // the multi-root algorithm gets stuck (for pathological cases like
            // Legendre polynomials), this shrinks the problem by 1 degree
            // and moves around the coefficients so they are not pathologic anymore
            let r = match recovery_algorithm {
                OneRootAlgorithms::Newton => this
                    .clone()
                    .one_root_newton(None, epsilon, max_recovery_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::Halley => this
                    .clone()
                    .one_root_halley(None, epsilon, max_recovery_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::Laguerre => this
                    .clone()
                    .one_root_laguerre(None, epsilon, max_recovery_iter)
                    .map_err(|_| roots.clone())?,
            };
            roots.push(r.clone());
            this = this / Poly::from_roots(&[r.clone()]);
        }

        Err(roots)
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use na::Complex;
    use num::complex::{Complex64, ComplexFloat};

    use crate::{poly::roots::OneRootAlgorithms, Poly, Poly64};

    #[test]
    fn initial_guess_smallest() {
        assert!(
            (poly![24.0, -14.0, -13.0, 2.0, 1.0].initial_guess_smallest()
                - Complex::new(0.68, 0.0))
            .norm()
                < 0.01
        );
    }

    #[test]
    fn roots_schur() {
        dbg!(poly![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
            .try_roots(0.01, 100, 10, None, None, None)
            .unwrap());
    }

    #[test]
    fn roots_newton() {
        let p = poly![1.0, 0.0, 1.0, 0.0, 1.0];

        // takes exactly 10 iterations
        let roots = p
            .try_n_roots(4, None, 1E-14, 10, Some(OneRootAlgorithms::Newton))
            .unwrap();
        assert!((Poly::from_roots(&roots) - p).almost_zero(&1E-14));
    }

    #[test]
    fn roots_halley() {
        let p = poly![1.0, 0.0, 1.0, 0.0, 1.0];

        // takes exactly 5 iterations
        let roots = p
            .try_n_roots(4, None, 1E-14, 5, Some(OneRootAlgorithms::Halley))
            .unwrap();
        assert!((Poly::from_roots(&roots) - p).almost_zero(&1E-14));
    }

    #[test]
    fn roots_laguerre() {
        let p = poly![1.0, 0.0, 1.0, 0.0, 1.0];

        // takes exactly 5 iterations
        let roots = p
            .try_n_roots(4, None, 1E-14, 100, Some(OneRootAlgorithms::Laguerre))
            .unwrap();
        assert!((Poly::from_roots(&roots) - p).almost_zero(&1E-14));
    }

    #[test]
    #[ignore]
    fn combinatorial_schur() {
        fn u16_to_poly(x: u16) -> Poly64 {
            let coeffs = (0..8)
                .map(|n| {
                    let re = ((x >> n) & 1) as f64;
                    let im = ((x >> (n + 1)) & 1) as f64;
                    Complex64::new(re, im)
                })
                .collect_vec();
            Poly64::from_complex_vec(coeffs)
        }
        for n in 0..u16::MAX {
            let poly = u16_to_poly(n);
            if poly.len() < 2 {
                continue;
            }
            match poly.clone().roots_schur(0.01, 30) {
                Err(_) => println!("{poly}"),
                _ => (),
            }
        }
    }

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn schur_roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.roots_schur(1E-14, 1000).unwrap();
        assert_eq!(roots[0].re(), -1.5);
        assert!((roots[0].im().abs() - 0.866) < 0.01);
        assert_eq!(roots[1].re(), -1.5);
        assert!((roots[1].im().abs() - 0.866) < 0.01);
    }

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn newton_roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly
            .try_n_roots(2, None, 1E-14, 1000, Some(OneRootAlgorithms::Newton))
            .unwrap();
        assert_eq!(roots[0].re(), -1.5);
        assert!((roots[0].im().abs() - 0.866) < 0.01);
        assert_eq!(roots[1].re(), -1.5);
        assert!((roots[1].im().abs() - 0.866) < 0.01);
    }
}
