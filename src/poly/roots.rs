use na::{Complex, ComplexField, RealField};
use num::{Float, FromPrimitive, One, Zero};

use crate::{
    Poly, Scalar, ScalarOps,
    __util::{
        self,
        complex::{c_min, c_neg},
    },
};

mod eigenvalue;
pub use eigenvalue::{EigenvalueRootFinder, FrancisQR};
mod iterative;
pub use iterative::{IterativeRootFinder, Naive, Newton};
mod initial_guess;
mod multiroot;

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
            Error::NoConverge(t) => Error::NoConverge(f(t)),
            Error::Other(o) => Error::Other(o),
        }
    }
}

pub type Result<T> = std::result::Result<Vec<Complex<T>>, Error<Vec<Complex<T>>>>;

/// Helper struct for implementing stateful root finders and converting between them
pub struct FinderState<T: Scalar> {
    pub poly: Poly<T>,

    /// Store roots that are within the given epsilon and are considered "done"
    pub clean_roots: Vec<Complex<T>>,

    /// Store "best guess" roots that are not within the epsilon bounds, but
    /// nevertheless can be useful as initial guesses for further refinement.
    pub dirty_roots: Vec<Complex<T>>,
}

/// Helper struct for implementing root finders and converting between them
pub struct FinderConfig<T: Scalar> {
    pub epsilon: T,
    pub max_iter: usize,
}

/// Statistics about the session (for benchmarking)
#[derive(Debug)]
pub struct FinderHistory<T: Scalar> {
    pub roots_history: Vec<Vec<Complex<T>>>,
}

impl<T: Scalar> FinderState<T> {
    fn new(poly: Poly<T>) -> Self {
        Self {
            poly,
            clean_roots: vec![],
            dirty_roots: vec![],
        }
    }
}

impl<T: Scalar> FinderConfig<T> {
    fn new() -> Self {
        Self {
            epsilon: T::small_safe(),
            max_iter: 0,
        }
    }
}

impl<T: Scalar> FinderHistory<T> {
    const fn new() -> Self {
        Self {
            roots_history: vec![],
        }
    }

    #[must_use]
    pub fn total_iter(&self) -> usize {
        self.roots_history.iter().map(std::vec::Vec::len).sum()
    }
}

/// Complex root finder for complex-valued polynomials
pub trait RootFinder<T: Scalar>: Sized {
    fn from_poly(poly: Poly<T>) -> Self;

    /// The maximum error within which a root is considered to be found.
    #[must_use]
    fn with_epsilon(mut self, epsilon: T) -> Self {
        self.config().epsilon = epsilon;
        self
    }

    #[must_use]
    fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.config().max_iter = max_iter;
        self
    }

    /// Initialize the state with some previously obtained guesses for roots.
    ///
    /// These will be used first before generating synthetic initial guesses
    /// and allow for recycling unsuccesful attempts from other root finders
    /// or if the user has some prior knowledge of the roots.
    ///
    /// Note that some finders may ignore these, but they will be carried over
    /// if multiple finders are used for the same polynomial.
    #[must_use]
    fn with_initial_guesses(mut self, guesses: Vec<Complex<T>>) -> Self {
        self.state_mut().dirty_roots.extend(guesses);
        self
    }

    #[must_use]
    fn collect_history(mut self) -> Self {
        if self.history().is_none() {
            *self.history() = Some(FinderHistory::new());
        }
        self
    }

    /// Find roots and store them in the state.
    fn run(&mut self) -> std::result::Result<(), Error<()>>;

    /// Tries to find as many roots as possible with the given configuration,
    /// returning a `std::result::Result` containing either the roots or the best guess so far
    /// in case of no convergence.
    ///
    /// Consecutive calls to `try_roots` will attempt to resume root finding
    /// from where it left off if possible. Use [`RootFinder::reset`] to
    /// create a new clean session.
    ///
    /// # Errors
    /// - [`Error::NoConverge`] - solver did not converge within `max_iter` iterations
    /// - [`Error::Other`] - solver encountered an unhandled edge-case
    fn roots(&mut self) -> Result<T> {
        self.run()
            .map(|()| self.state().clean_roots.clone())
            .map_err(|e| e.map_no_converge(|()| self.state().dirty_roots.clone()))
    }

    /// Get a mutable reference to the current finder state
    fn state_mut(&mut self) -> &mut FinderState<T>;

    /// Get a reference to the current finder state
    fn state(&self) -> &FinderState<T>;

    /// Get a mutable reference to the current finder configuration
    fn config(&mut self) -> &mut FinderConfig<T>;

    /// Get a mutable reference to the history
    fn history(&mut self) -> &mut Option<FinderHistory<T>>;
}

/// Polynomial root-finding algorithms
#[non_exhaustive]
pub enum OneRootAlgorithms {
    Newton,
    Halley,
    JenkinsTraub,
}

#[non_exhaustive]
pub enum AllRootsAlgorithms {
    #[deprecated]
    Schur,

    FrancisQR,
}

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
        FrancisQR::from_poly(self.clone())
            .with_epsilon(epsilon)
            .with_max_iter(max_iter)
            .with_companion_matrix_type(eigenvalue::CompoanionMatrixType::Schmeisser)
            .roots()
    }
}

// private
impl<T: ScalarOps + RealField + Float> Poly<T> {
    fn trivial_roots(&mut self, epsilon: T) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());

        let mut roots = vec![];
        for _ in 0..self.degree_raw() {
            if self.eval_point(Complex::zero()).norm() < epsilon {
                roots.push(Complex::zero());
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

        roots
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
        *self = Poly::one();

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
        *self = Poly::one();

        vec![x1, x2]
    }
}

#[cfg(test)]
mod test {

    use na::Complex;
    use num::complex::{Complex64, ComplexFloat};

    use crate::__util::testing::{binary_coeffs, check_roots};
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

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn schur_roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.roots(1E-14, 1000).unwrap();
        assert!((roots[0].re() - -1.5).abs() < 0.01);
        assert!(
            (roots[0].im().abs() - 0.866).abs() < 0.01,
            "{}",
            roots[0].im()
        );
        assert!((roots[1].re() - -1.5).abs() < 0.01);
        assert!(
            (roots[1].im().abs() - 0.866).abs() < 0.01,
            "{}",
            roots[1].im()
        );
    }

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn newton_roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.roots(1E-14, 1000).unwrap();
        assert!((roots[0].re() - -1.5).abs() < 0.01);
        assert!((roots[0].im().abs() - 0.866).abs() < 0.01);
        assert!((roots[1].re() - -1.5).abs() < 0.01);
        assert!((roots[1].im().abs() - 0.866).abs() < 0.01);
    }
}
