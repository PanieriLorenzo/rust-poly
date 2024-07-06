use crate::{Poly, Scalar, ScalarOps, __util::complex::c_neg};
use na::{Complex, ComplexField, RealField};
use num::{traits::real::Real, Float, FromPrimitive, One, Zero};

mod single_root;
pub use single_root::{halley, naive, newton};
mod all_roots;
pub use all_roots::{aberth_ehrlich, deflate, halley_deflate, naive_deflate, newton_deflate};
mod many_roots;
pub use many_roots::{halley_parallel, naive_parallel, newton_parallel, parallel};

mod initial_guess;
pub use initial_guess::{initial_guess_smallest, initial_guesses_circle, initial_guesses_random};

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error<T> {
    #[error("root finder did not converge within the given constraints")]
    NoConverge(T),

    #[error("unexpected error while running root finder")]
    Other(#[from] anyhow::Error),
}

impl<T> Error<T> {
    pub(crate) fn map_no_converge<U>(self, mut f: impl FnMut(T) -> U) -> Error<U> {
        match self {
            Self::NoConverge(t) => Error::NoConverge(f(t)),
            Self::Other(o) => Error::Other(o),
        }
    }
}

pub type Result<T> = std::result::Result<Vec<Complex<T>>, Error<Vec<Complex<T>>>>;

impl<T: ScalarOps + RealField + Float> Poly<T> {
    /// A convenient way of finding roots, with a pre-configured root finder.
    /// Should work well for most real polynomials of low degree.
    ///
    /// Use a root finder if you need more control over performance or accuracy.
    ///
    /// # Errors
    /// - Solver did not converge within `max_iter` iterations
    /// - Some other edge-case was encountered which could not be handled (please
    ///   report this, as we can make this solver more robust!)
    pub fn roots(&self, epsilon: T, max_iter: usize) -> Result<T> {
        debug_assert!(self.is_normalized());
        let roots = aberth_ehrlich(&mut self.clone(), Some(epsilon), Some(max_iter), &[])?;

        // don't bother refining small polys
        if roots.len() <= 2 {
            return Ok(roots);
        }

        // further polishing of roots
        newton_parallel(
            &mut self.clone(),
            Some(epsilon),
            Some(max_iter),
            roots.len(),
            &roots,
        )
    }
}

// private
impl<T: ScalarOps + RealField + Float> Poly<T> {
    fn trivial_roots(&mut self, epsilon: T) -> (Vec<Complex<T>>, u128) {
        let mut eval_counter = 0;
        debug_assert!(self.is_normalized());

        let mut roots = vec![];
        for _ in 0..self.degree_raw() {
            if self.eval(Complex::zero()).norm() < epsilon {
                eval_counter += 1;
                roots.push(Complex::zero());
                // TODO: deflating zero roots can be accomplished simply by shifting
                *self = self.clone().deflate_composite(Complex::zero());
            } else {
                break;
            }
        }

        match self.degree_raw() {
            1 => roots.extend(self.linear()),
            2 => roots.extend(self.quadratic()),
            _ => {}
        }

        // post-condition: polynomials of degree 1 or 2 have been reduced
        debug_assert!(self.degree_raw() != 1 && self.degree_raw() != 2);
        (roots, eval_counter)
    }

    fn linear(&mut self) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());
        debug_assert_eq!(self.degree_raw(), 1);

        self.trim();
        if self.degree_raw() < 1 {
            return vec![];
        }

        let a = self.0[1];
        let b = self.0[0];

        // we found all the roots
        *self = Self::one();

        vec![-b / a]
    }

    /// Quadratic formula
    fn quadratic(&mut self) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());
        debug_assert_eq!(self.degree_raw(), 2);

        // trimming trailing almost zeros to avoid overflow
        self.trim();
        if self.degree_raw() == 1 {
            return self.linear();
        }
        if self.degree_raw() == 0 {
            return vec![];
        }

        let a = self.0[2];
        let b = self.0[1];
        let c = self.0[0];
        let four = Complex::<T>::from_u8(4).expect("overflow");
        let two = Complex::<T>::from_u8(2).expect("overflow");

        // TODO: switch to different formula when b^2 and 4c are very close due
        //       to loss of precision
        let plus_minus_term = (b * b - four * a * c).sqrt();
        let x1 = (plus_minus_term - b) / (two * a);
        let x2 = (c_neg(b) - plus_minus_term) / (two * a);

        // we found all the roots
        *self = Self::one();

        vec![x1, x2]
    }
}

#[cfg(test)]
mod test {
    use num::complex::ComplexFloat;

    use crate::Poly64;

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.roots(1E-10, 1000).unwrap();
        assert!((roots[0].re() - -1.5).abs() < 0.01);
        assert!((roots[0].im().abs() - 0.866).abs() < 0.01);
        assert!((roots[1].re() - -1.5).abs() < 0.01);
        assert!((roots[1].im().abs() - 0.866).abs() < 0.01);
    }
}
