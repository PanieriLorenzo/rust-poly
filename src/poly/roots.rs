use na::{Complex, ComplexField, Normed, RealField};
use num::{
    traits::{float::FloatCore, MulAdd},
    Float, FromPrimitive, Num, One, Zero,
};

use crate::{
    Poly, Scalar, ScalarOps,
    __util::{
        self,
        casting::usize_to_scalar,
        complex::{c_min, c_neg},
    },
};

mod eigenvalue;
mod iterative;
pub use iterative::{newton::NewtonFinder, IterativeRootFinder};
mod multiroot;

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error<T> {
    #[error("root finder did not converge within the given constraints")]
    NoConverge(Vec<Complex<T>>),

    #[error("unexpected error while running root finder")]
    Other,
}

pub type Result<T> = std::result::Result<Vec<Complex<T>>, Error<T>>;

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

/// Complex root finder for complex-valued polynomials
pub trait RootFinder<T: Scalar>: Sized {
    fn from_poly(poly: Poly<T>) -> Self;

    /// The maximum error within which a root is considered to be found.
    fn with_epsilon(mut self, epsilon: T) -> Self {
        self.config().epsilon = epsilon;
        self
    }

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
    fn with_initial_guesses(mut self, guesses: Vec<Complex<T>>) -> Self {
        self.state().dirty_roots.extend(guesses);
        self
    }

    /// Tries to find as many roots as possible with the given configuration,
    /// returning a `std::result::Result` containing either the roots or the best guess so far
    /// in case of no convergence.
    ///
    /// Consecutive calls to `try_roots` will attempt to resume root finding
    /// from where it left off if possible. Use [`RootFinder::reset`] to
    /// create a new clean session.
    fn roots(&mut self) -> Result<T>;

    /// Get a mutable reference to the current finder state
    fn state(&mut self) -> &mut FinderState<T>;

    /// Get a mutable reference to the current finder configuration
    fn config(&mut self) -> &mut FinderConfig<T>;
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

// private
impl<T: ScalarOps + RealField + Float> Poly<T> {
    fn linear(mut self) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());
        debug_assert_eq!(self.degree_raw(), 1);

        self.trim();
        if self.degree_raw() < 1 {
            return vec![];
        }

        let a = self.0[1];
        let b = self.0[0];

        vec![-b / a]
    }

    /// Quadratic formula
    fn quadratic(mut self) -> Vec<Complex<T>> {
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
        let four = Complex::<T>::from_u8(4).expect("should always fit for small ints");
        let two = Complex::<T>::from_u8(2).expect("should always fit for small ints");

        // TODO: switch to different formula when b^2 and 4c are very close due
        //       to loss of precision
        let plus_minus_term = (b * b - four * a * c).sqrt();
        let x1 = (plus_minus_term - b) / (two * a);
        let x2 = (c_neg(b) - plus_minus_term) / (two * a);
        vec![x1, x2]
    }

    /// Ref: https://doi.org/10.1007/BF01933524
    pub(crate) fn initial_guess_smallest(&self) -> Complex<T> {
        debug_assert!(self.is_normalized());
        debug_assert!(self.len_raw() >= 2);

        let small = Float::recip(usize_to_scalar::<T>(1_000));
        let p_diff = self.clone().diff();
        let mut pz = self.eval_point(Complex::zero());
        let mut pdz = p_diff.eval_point(Complex::zero());

        // avoid divide by zero
        if pdz.norm() < small {
            pz += small;
            pdz += small;
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
            guess += Complex::i().scale(Float::recip(usize_to_scalar::<T>(1_000)));
        }
        guess
    }

    fn one_root_newton(
        &self,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
    ) -> std::result::Result<Complex<T>, Complex<T>> {
        let p_diff = self.clone().diff();
        let mut x = initial_guess.unwrap_or(self.initial_guess_smallest());
        for _ in 0..max_iter {
            let px = self.eval_point(x);
            if px.norm() <= epsilon {
                return Ok(x);
            }
            let pdx = p_diff.eval_point(x);
            x -= px / pdx;
        }
        Err(x)
    }

    fn one_root_halley(
        &self,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
    ) -> std::result::Result<Complex<T>, Complex<T>> {
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
            x -= (px * pdx * two) / (pdx.powu(2) * two - px * pddx);
        }
        Err(x)
    }

    fn one_root_jenkins_traub(
        &mut self,
        epsilon: T,
        max_iter: usize,
    ) -> std::result::Result<Complex<T>, Complex<T>> {
        // TODO: tune these to the size of the polynomial with a lookup table
        const M: usize = 5;
        const L: usize = 100;

        self.make_monic();

        // stage one
        let mut h_poly = self.clone().diff();
        for _ in 0..M {
            let pz = self.eval_point(Complex::zero());

            // TODO: too many clones
            let hz = h_poly.clone().eval_point(Complex::zero());
            h_poly = h_poly - self.clone().scaled(hz / pz);

            // TODO: linear division can be done with synthetic division
            h_poly = h_poly / poly![T::zero(), T::one()];
        }

        // stage two
        todo!();

        // stage three
        todo!()
    }

    fn roots_francis_qr(
        &self,
        epsilon: T,
        max_iter: usize,
    ) -> std::result::Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        // TODO: this implementation uses an outdated Francis shift algorithm,
        //       the "state of the art" (from the 90s lmao), is to use a
        //       multishift QR algorithm. This is what LAPACK does.
        //       read the papers by Karen Braman, Ralph Byers and Roy Mathias
        debug_assert!(self.is_normalized());
        debug_assert!(
            self.degree_raw() > 2,
            "cannot use roots_francis_qr method on polynomials of degree 2 or smaller"
        );

        // TODO: tune max_iter_per_deflation to input

        // a small safe number, this is often done in LAPACK to avoid overflows
        let small_num = T::small_safe();

        // TODO: remove this clone
        let mut this = self.clone();

        let mut needs_unshifting = false;
        // leading zero coefficients can be problematic, so shifting to
        // avoid (the householder reflectors cannot be constructed if the
        // first element is zero)
        if this.0[0].abs() < small_num {
            needs_unshifting = true;
            this = this.translate(Complex::one(), Complex::zero());
        }

        let mut comp = this.companion();

        // rotating the matrix 180 degrees. This is equivalent to using similarity
        // transforms so it does not move the eigenvalues. But NumPy does it
        // and apparently it makes it more precise
        let n = comp.nrows();
        for i in 0..n / 2 {
            comp.swap_rows(i, n - i - 1);
            comp.swap_columns(i, n - i - 1);
        }

        let mut roots =
            __util::linalg::eigen_francis_shift(comp.as_view_mut(), epsilon, max_iter, max_iter)?;

        if needs_unshifting {
            for r in &mut roots {
                *r -= Complex::one();
            }
        }

        Ok(roots)
    }
}

impl<T: ScalarOps + Float + RealField> Poly<T> {
    /// Find only some of the roots of the polynomial.
    ///
    /// Note that for large `n`, using [`Poly::try_roots`] is probably faster.
    ///
    /// It utilizes an iterative method, so the precision gets progressively
    /// worse the more roots are found. For small `n` this is negligible.
    ///
    /// `Err` `std::result::Result` contains the roots it was able to find, even if they are
    /// fewer than requested.
    ///
    /// Use [`Poly::try_n_roots_algo`] to specify which algorithm to use, if
    /// you already know which one will perform best.
    #[deprecated = "use NewtonFinder instead"]
    pub fn try_n_roots(
        &self,
        n: usize,
        initial_guess: Option<Complex<T>>,
        epsilon: T,
        max_iter: usize,
        algorithm: Option<OneRootAlgorithms>,
    ) -> std::result::Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        debug_assert!(self.is_normalized());
        assert!(
            n as i32 <= self.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

        // TODO: if you use monic polynomials, there is a simpler algorithm
        //       for polynomial division that improves accuracy, because no
        //       divisions are performed on the coefficients

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
                OneRootAlgorithms::JenkinsTraub => unimplemented!(),
            };
            roots.push(r);
            if i < (n - 1) {
                this = this / Self::from_roots(&[r]);
            }
        }
        Ok(roots)
    }

    #[deprecated = "use NewtonFinder instead"]
    pub fn try_roots(
        &self,
        epsilon: T,
        max_iter: usize,
        max_tries: usize,
        max_recovery_iter: Option<usize>,
        algorithm: Option<AllRootsAlgorithms>,
        recovery_algorithm: Option<OneRootAlgorithms>,
    ) -> std::result::Result<Vec<Complex<T>>, Vec<Complex<T>>> {
        debug_assert!(self.is_normalized());

        let max_recovery_iter = max_recovery_iter.unwrap_or(max_iter);
        let algorithm = algorithm.unwrap_or(AllRootsAlgorithms::FrancisQR);
        let recovery_algorithm = recovery_algorithm.unwrap_or(OneRootAlgorithms::Newton);

        let mut roots = vec![];
        let mut this = self.clone();

        for _ in 0..max_tries {
            this.trim();

            match self.degree_raw() {
                ..=0 => return Ok(roots),
                1 => {
                    roots.extend(this.linear());
                    return Ok(roots);
                }
                2 => {
                    roots.extend(this.quadratic());
                    return Ok(roots);
                }
                _ => {}
            }

            let maybe_roots = match algorithm {
                AllRootsAlgorithms::Schur | AllRootsAlgorithms::FrancisQR => {
                    this.roots_francis_qr(epsilon, max_iter)
                }
            };

            if let Ok(found_roots) = maybe_roots {
                roots.extend(found_roots);
                // TODO: sort
                return Ok(roots);
            }
            let err = maybe_roots.expect_err("infallible");
            // TODO: use recovered eigens as initial guesses

            // uses one iteration of single-root algorithm for recovery when
            // the multi-root algorithm gets stuck (for pathological cases like
            // Legendre polynomials), this shrinks the problem by 1 degree
            // and moves around the coefficients so they are not pathological anymore
            let r = match recovery_algorithm {
                OneRootAlgorithms::Newton => this
                    .clone()
                    .one_root_newton(None, epsilon, max_recovery_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::Halley => this
                    .clone()
                    .one_root_halley(None, epsilon, max_recovery_iter)
                    .map_err(|_| roots.clone())?,
                OneRootAlgorithms::JenkinsTraub => unimplemented!(),
            };
            roots.push(r);
            this = this / Self::from_roots(&[r]);
        }

        Err(roots)
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use na::Complex;
    use num::complex::{Complex64, ComplexFloat};

    use crate::__util::{
        complex::complex_sort_mut,
        testing::{binary_coeffs, check_roots},
    };
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
        let mut roots_expected = vec![
            Complex64 { re: 0.0, im: 0.0 },
            Complex64 { re: 1.0, im: 0.0 },
            Complex64 { re: 2.5, im: 0.0 },
            // Complex64 { re: 1.0, im: 1.0 },
            // Complex64 { re: 1.0, im: -2.5 },
            // Complex64 { re: -1.0, im: 1.0 },
            // Complex64 { re: -1.0, im: -1.0 },
        ];
        let poly = Poly::from_roots(&roots_expected);
        let mut roots = poly.roots_francis_qr(1E-9, 1000).unwrap();

        assert!(
            check_roots(roots.clone(), roots_expected.clone(), 1E-4),
            "{:?} != {:?}",
            roots,
            roots_expected
        );
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

    /// This test is to find the minimum number of iterations at various degrees
    /// that achieves 90% success rate, using "binary coefficients" polynomials,
    /// which are difficult for the nalgebra implementation of the Schur
    /// algorithm, using a somewhat realistic epsilon of `1E-5`. This is a heuristic
    /// for deciding when schur did not converge and should switch to single-root
    /// strategy instead.
    ///
    /// If the [LAPACK implementation](https://netlib.org/lapack/explore-html/d5/d38/group__gees_ga59d7d13222ddc67c5ae0e9aa1a62cb61.html#ga59d7d13222ddc67c5ae0e9aa1a62cb61)
    /// was used, this would not be necessary, because the LAPACK implementation
    /// is more robust. But nalgebra has had this issue in a long time and they
    /// don't seem to be working on fixing it at the moment. See nalgebra issues
    /// [nalgebra-#1291](https://github.com/dimforge/nalgebra/issues/1291),
    /// [nalgebra-#764](https://github.com/dimforge/nalgebra/issues/764),
    /// [nalgebra-#611](https://github.com/dimforge/nalgebra/issues/611).
    #[test]
    #[ignore]
    fn schur_tuning() {
        fn scenario(deg: usize, iter: usize) {
            let mut ok = 0;
            let mut err = 0;
            for p in binary_coeffs(deg as i32, deg) {
                match p.try_roots(1E-5, iter, 1, None, None, None) {
                    Ok(_) => ok += 1,
                    Err(_) => err += 1,
                }
            }
            assert!(ok as f32 / (ok + err) as f32 >= 0.9, "{}/{}", ok, ok + err);
        }
        scenario(4, 100);
        // scenario(5, 49);
        // scenario(6, 59);
        // scenario(7, 36);
        // scenario(8, 23);
        // scenario(9, 23);
        // scenario(10, 25);
        //scenario(15, 25);
    }

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn schur_roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.try_roots(1E-14, 1000, 1, None, None, None).unwrap();
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
