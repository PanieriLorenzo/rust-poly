#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]
use duplicate::duplicate_item;
use ndarray::{array, s, Array, Array1, ArrayView1, ArrayViewMut1, Axis, Dimension};
use num_complex::Complex;
use num_rational::Ratio;
use num_traits::{Float, One, Zero};

mod scalar;
pub use scalar::Scalar;

// mod roots;
// pub use roots::Roots;

mod impl_num;

/// polynomial as a list of coefficients of terms of descending degree
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Poly<T: Scalar>(Array1<Complex<T>>);

pub(crate) type C<T> = Complex<T>;
pub(crate) type A<T> = Array1<T>;
pub(crate) type AC<T> = Array1<Complex<T>>;

impl<T: Scalar> Poly<T> {
    /// Create a new polynomial from a 1D array of complex coefficients
    #[must_use]
    pub fn new(coeffs: AC<T>) -> Self {
        Self(coeffs).trim_zeros()
    }

    /// Create a new polynomial from a 1D array of real coefficients
    #[must_use]
    pub fn from_reals(coeffs: A<T>) -> Self {
        Self(coeffs.map(|x| Complex::new(x.clone(), T::zero()))).trim_zeros()
    }

    /// Creates a polynomial with a single term of degree `n`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let t1 = Poly::term(Complex::from(1i32), 0);
    /// let t2 = Poly::term(Complex::from(2i32), 1);
    /// let t3 = Poly::term(Complex::from(3i32), 2);
    /// assert_eq!(t1, Poly::new(array![Complex::from(1)]));
    /// assert_eq!(t2, Poly::new(array![Complex::from(2), Complex::from(0)]));
    /// assert_eq!(t3, Poly::new(array![Complex::from(3), Complex::from(0), Complex::from(0)]));
    /// ```
    #[allow(clippy::missing_panics_doc)]
    pub fn term(coeff: C<T>, degree: usize) -> Self {
        let zeros: AC<T> = AC::<T>::zeros([degree]);
        let mut term: AC<T> = array![coeff];
        term.append(Axis(0), zeros.view())
            .unwrap_or_else(|_| unreachable!()); // TODO: proof
        Self(term)
    }

    /// Checks whether the polynomial has leading zeros
    fn is_normalized(&self) -> bool {
        if self.raw_len() == 0 {
            return true;
        }
        !self.0[0].is_zero()
    }

    /// Removes leading zero coefficients
    fn trim_zeros(&self) -> Self {
        if self.is_normalized() {
            return self.clone();
        }
        let mut first: usize = 0;
        for e in &self.0 {
            if !e.is_zero() {
                break;
            }
            first += 1;
        }
        Self(self.0.slice(s![first..]).to_owned())
    }

    // TODO: trim in-place for better performance

    /// Length of the polynomial
    ///
    /// Note that this does not include leading zeros, as polynomials are
    /// stored in their normalized form internally.

    // internal NOTE: strictly speaking, polynomials are only normalized
    //     when necessary, but this *should* be invisible to the user
    #[must_use]
    pub fn len(&self) -> usize {
        self.trim_zeros().raw_len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn degree(&self) -> i32 {
        (self.len() - 1) as i32
    }

    /// Length of the polynomial, without trimming zeros
    fn raw_len(&self) -> usize {
        self.0.len()
    }

    /// Evaluate a polynomial at a specific input value `x`. This may be an
    /// ndarray of any dimension
    ///
    /// # Examples
    ///
    /// Evaluate a real polynomial at real points
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// // x^2 + 2x + 1
    /// let p = Poly::new(array![Complex::from(1), Complex::from(2), Complex::from(1)]);
    /// let x = array![Complex::from(-1), Complex::from(0), Complex::from(1)];
    /// let y = p.eval(&x);
    /// assert_eq!(y, array![Complex::from(0), Complex::from(1), Complex::from(4)]);
    /// ```
    ///
    /// Evaluate a complex polynomial at complex points
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex64;
    ///
    /// // (2+i)x^2 + 2i
    /// let p = Poly::new(array![
    ///     Complex64::new(2.0, 1.0),
    ///     Complex64::new(0.0, 0.0),
    ///     Complex64::new(0.0, 2.0),
    /// ]);
    /// let x = array![Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0)];
    /// let y = p.eval(&x);
    /// assert_eq!(y, array![Complex64::new(2.0, 3.0), Complex64::new(-2.0, 1.0)]);
    /// ```
    pub fn eval<D: Dimension>(&self, x: &Array<C<T>, D>) -> Array<C<T>, D> {
        let mut y: Array<C<T>, D> = Array::<C<T>, D>::zeros(x.raw_dim());
        for pv in &self.0 {
            y = (y * x).map(|e| e.clone() + pv.clone());
        }
        y
    }

    #[must_use]
    pub fn as_ndarray(&self) -> ArrayView1<C<T>> {
        self.0.view()
    }

    pub fn as_ndarray_mut(&mut self) -> ArrayViewMut1<C<T>> {
        self.0.view_mut()
    }

    #[must_use]
    pub fn to_vec(&self) -> Vec<C<T>> {
        self.0.to_vec()
    }

    /// Computes the quotient and remainder of a polynomial division
    ///
    /// # Examples
    ///
    /// Divide two real polynomials
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::new(array![Complex::from(3.0), Complex::from(5.0), Complex::from(2.0)]);
    /// let p2 = Poly::new(array![Complex::from(2.0), Complex::from(1.0)]);
    /// let (q, r) = p1.div_rem(&p2);
    /// assert_eq!(q, Poly::new(array![Complex::from(1.5), Complex::from(1.75)]));
    /// assert_eq!(r, Poly::new(array![Complex::from(0.25)]));
    /// ```
    ///
    /// Divide two complex polynomials
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::term(Complex::new(1.0, 1.0), 2);
    /// let p2 = Poly::term(Complex::new(1.0, -1.0), 0);
    /// let (q, r) = p1.div_rem(&p2);
    /// assert_eq!(q, Poly::term(Complex::new(0.0, 1.0), 2));
    /// assert_eq!(r, Poly::new(array![]));
    /// ```
    ///
    /// # Panics
    ///
    /// Dividing by zero is not allowed! Even when using floating point numbers.
    /// Arithmetically, following the algorithm for long division here would just
    /// result in `+inf` or `-inf` quotient and `nan` remainder, but this can lead
    /// to hard to debug errors, so we prefer to outright disallow any division by zero
    /// with a more helpful error message.
    ///
    /// ```should_panic
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::new(array![Complex::from(1.0), Complex::from(2.0), Complex::from(3.0)]);
    /// let p2 = Poly::new(array![]);
    /// let (q, r) = p1.div_rem(&p2);
    /// ```
    ///
    #[must_use]
    pub fn div_rem(self, rhs: &Self) -> (Self, Self) {
        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }

        assert!(!rhs.is_zero(), "Attempted to divide polynomial by zero");

        let num: AC<T> = self.0;
        let den: AC<T> = rhs.0.clone();

        // cannot underflow because we have ensured that len is at least 1
        let num_deg = num.len() - 1;
        let den_deg = den.len() - 1;

        let scale = C::<T>::one() / den[0].clone();
        let mut quot: AC<T> = AC::<T>::zeros((num_deg - den_deg + 1).max(1));
        let mut rem: AC<T> = num;
        for k in 0..=(num_deg - den_deg) {
            let d = scale.clone() * rem[k].clone();
            quot[k] = d.clone();
            rem.slice_mut(s![k..=(k + den_deg)])
                .iter_mut()
                .zip((array![d] * den.clone()).iter())
                .for_each(|p| *p.0 = p.0.clone() - p.1.clone());
        }
        (Self(quot), Self(rem).trim_zeros())
    }
}
