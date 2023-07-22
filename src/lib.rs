#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]
use duplicate::duplicate_item;
pub use ndarray;
use ndarray::{array, s, Array, Array1, Array2, ArrayView1, ArrayViewMut1, Axis, Dimension};
pub use num_complex;
use num_complex::Complex;
use num_rational::Ratio;
pub use num_traits;
use num_traits::{Float, One, Zero};

mod scalar;
pub use scalar::Scalar;

// mod roots;
// pub use roots::Roots;

mod complex_util;
use complex_util::c_neg;
mod array_util;
use array_util::np_diag;
mod impl_num;

/// A more convenient way to write `Complex::new(...)`.
///
/// # Examples
///
/// ```
/// # use rust_poly::c;
/// use num_complex::Complex;
///
/// let c1: Complex<f32> = c!();
/// let c2 = Complex::new(0.0, 0.0);
/// let c3 = c!(1.0f32, 2.0);
/// let c4 = Complex::new(1.0, 2.0);
///
/// assert_eq!(c1, c2);
/// assert_eq!(c3, c4);
/// ```
#[macro_export]
macro_rules! c {
    () => {{
        <$crate::num_complex::Complex<_> as $crate::num_traits::Zero>::zero()
    }};
    ($re:expr, $im: expr) => {{
        $crate::num_complex::Complex::new($re, $im)
    }};
}

/// A more convenient way of writing `Poly::from(array![...])`
///
/// It takes ownership of its arguments.
///
/// It can take a list of `Scalar` or `Complex<Scalar>`. If left empty, it is
/// equivalent to `Poly::zero()`.
///
/// # Examples
///
/// Basic syntax
/// ```
/// # use rust_poly::{poly, Poly};
/// use num_traits::Zero;
/// use num_complex::Complex;
///
/// let p1: Poly<f32> = poly![];
/// let p2 = poly![1.0f32, 2.0, 3.0];
/// let p3 = poly![Complex::from(1.0), Complex::from(2.0), Complex::from(3.0)];
///
/// assert_eq!(p1, Poly::zero());
/// assert_eq!(p2, p3);
/// ```
///
/// Similarly to `vec!`, you can initialize a large polynomial where all coefficients
/// are equal like so:
/// ```
/// # use rust_poly::{poly, Poly};
/// use num_complex::Complex;
///
/// let p1 = poly![2.0; 16];
/// let p2 = poly![Complex::from(2.0); 16];
///
/// assert_eq!(p1, p2);
/// ```
///
/// You can also express complex numbers as a tuple of two scalars, mixing and matching
/// this syntax with the other syntax rules:
/// ```
/// # use rust_poly::{poly, Poly};
/// use num_complex::Complex;
///
/// let p1 = poly![(1.0, 2.0), (1.0, 2.0)];
/// let p2 = poly![(1.0, 2.0); 2];
/// let p3 = poly![Complex::new(1.0, 2.0); 2];
/// let p4 = poly![Complex::new(1.0, 2.0), Complex::new(1.0, 2.0)];
///
/// assert_eq!(p1, p2);
/// assert_eq!(p1, p3);
/// assert_eq!(p1, p4);
/// ```
#[macro_export]
macro_rules! poly {
    () => {{
        $crate::Poly::zero()
    }};
    (($re:expr, $im:expr); $n:expr) => {{
        $crate::Poly::from($crate::ndarray::Array1::from_vec(vec![$crate::c!($re, $im); $n]))
    }};
    ($elem:expr; $n:expr) => {{
        $crate::Poly::from($crate::ndarray::Array1::from_vec(vec![$elem; $n]))
    }};
    ($(($re:expr, $im:expr)),+ $(,)?) => {{
        $crate::Poly::from($crate::ndarray::array![$($crate::c!($re, $im)),*])
    }};
    ($($elems:expr),+ $(,)?) => {{
        $crate::Poly::from($crate::ndarray::array![$($elems),*])
    }};
}

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

    /// Create a polynomial from its roots
    #[must_use]
    pub fn from_roots(roots: AC<T>) -> Self {
        todo!();
    }

    /// Create a real polynomial from its real roots
    #[must_use]
    pub fn from_real_roots(roots: A<T>) -> Self {
        todo!();
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
    ///
    /// The resulting polynomial is equivalent, but is "normalized", this
    /// makes finding the degree of the polynomial as easy as just counting
    /// the coefficients.
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

    /// Removes trailing zero coefficients
    ///
    /// *For normalization, use `trim_zeros` instead!*
    ///
    /// **WARNING**: the result of this operation is not equivalent. All
    /// terms of the polynomial decrease in degree. Mostly used for specific
    /// algorithms.
    fn trim_trailing_zeros(&self) -> Self {
        let mut last: usize = self.raw_len() - 1;
        for e in self.0.iter().rev() {
            if !e.is_zero() {
                break;
            }
            last -= 1;
        }
        Self(self.0.slice(s![..=last]).to_owned())
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
        self.len() as i32 - 1
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

    /// Solve the equation `P = 0` numerically, where `P` is this polynomial.
    ///
    /// This will always succeed, but the solutions might be complex, even for
    /// real polynomials.
    ///
    /// Note that roots that are far from the origin of have multiplicity higher
    /// than 1 may have low accuracy, for better accuracy, use `roots_precise`.
    pub fn roots(&self) -> AC<T> {
        let coeffs: AC<T> = self.trim_trailing_zeros().0;
        let n = coeffs.len();
        let roots = if n > 1 {
            // build companion matrix and find its eigenvalues (the roots)
            let a = Array2::<Complex<T>>::from_diag(A::<C<T>>::ones([n-2]))
        }
    }

    /// A more numerically precise version of `roots`.
    ///
    /// It uses `roots` as the initial guess for several Newton's method iterations.
    pub fn roots_precise(&self) -> AC<T> {
        todo!()
    }

    // TODO: real polynomial decomposition into degree 1 and 2 polynomials
}

impl<T: Scalar> From<A<T>> for Poly<T> {
    fn from(value: A<T>) -> Self {
        Self::from_reals(value)
    }
}

impl<T: Scalar> From<AC<T>> for Poly<T> {
    fn from(value: AC<T>) -> Self {
        Self::new(value)
    }
}
