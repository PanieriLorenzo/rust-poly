use na::{Complex, ComplexField, Normed, RealField, Scalar};
use num::{traits::float::FloatCore, Float, One, Zero};

use crate::{
    util::{
        casting::usize_to_scalar,
        complex::{c_min, c_neg},
    },
    Poly, ScalarOps,
};

/// Polynomial root-finding algorithms
#[non_exhaustive]
pub enum OneRootAlgorithms {
    Newton,
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
            if px.norm_sqr() <= epsilon {
                return Ok(x);
            }
            let pdx = p_diff.eval_point(x);
            x = x - px / pdx;
        }
        Err(x)
    }

    fn one_root_halley(&self, epsilon: T, max_iter: usize) -> Result<Complex<T>, Complex<T>> {
        todo!()
    }

    fn one_root_laguerre(&self, epsilon: T, max_iter: usize) -> Result<Complex<T>, Complex<T>> {
        todo!()
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
        let one_root = match algorithm {
            OneRootAlgorithms::Newton => {
                |this: Poly<T>| this.one_root_newton(initial_guess, epsilon, max_iter)
            }
            _ => unimplemented!(),
        };

        let mut roots = vec![];
        let mut this = self.clone();
        for i in 0..n {
            let r = one_root(this.clone()).map_err(|_| roots.clone())?;
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
        let one_root = match recovery_algorithm {
            OneRootAlgorithms::Newton => {
                |this: Poly<T>| this.one_root_newton(None, epsilon, max_recovery_iter)
            }
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
                return Ok(roots);
            }

            // uses one iteration of single-root algorithm for recovery when
            // the multi-root algorithm gets stuck (for pathological cases like
            // Legendre polynomials)
            let r = one_root(this.clone()).map_err(|_| roots.clone())?;
            roots.push(r.clone());
            this = this / Poly::from_roots(&[r.clone()]);
        }

        Err(roots)
    }
}

#[cfg(test)]
mod test {
    use na::{Complex, ComplexField};

    use crate::Poly;

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
}
