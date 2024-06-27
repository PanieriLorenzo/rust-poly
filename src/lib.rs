// TODO(version: v1.0.0): license/author header project-wide, see MIT guidelines
#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]

extern crate nalgebra as na;
use std::ops::Index;

use bessel::coeff;
use na::Normed;
pub use num;

use num::{Complex, One, Zero};

/// A more convenient way to write `Complex::new(...)`.
///
/// # Examples
///
/// ```
/// use rust_poly::complex;
/// use num::Complex;
///
/// let c1: Complex<f32> = complex!();
/// let c2 = Complex::new(0.0, 0.0);
/// let c3 = complex!(1.0f32, 2.0);
/// let c4 = Complex::new(1.0, 2.0);
///
/// assert_eq!(c1, c2);
/// assert_eq!(c3, c4);
/// assert_eq!(complex!(4.20), complex!(4.20, 0.0));
/// ```
#[macro_export]
macro_rules! complex {
    () => {{
        <$crate::num::Complex<_> as $crate::num::Zero>::zero()
    }};
    ($re:expr) => {{
        $crate::num::Complex::new(
            $re,
            <$crate::num::Complex<_> as $crate::num::Zero>::zero().im,
        )
    }};
    ($re:expr, $im: expr) => {{
        $crate::num::Complex::new($re, $im)
    }};
}

/// A more convenient way of writing `Poly::new(&[Complex::new(...)...])`
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
/// use rust_poly::{poly, Poly};
/// use num::Zero;
/// use num::Complex;
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
/// use num::Complex;
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
/// use rust_poly::{poly, Poly};
/// use num::Complex;
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
        $crate::Poly::from(vec![$crate::complex!($re, $im); $n])
    }};
    ($elem:expr; $n:expr) => {{
        $crate::Poly::from(vec![$elem; $n])
    }};
    ($(($re:expr, $im:expr)),+ $(,)?) => {{
        $crate::Poly::from(vec![$($crate::complex!($re, $im)),*])
    }};
    ($($elems:expr),+ $(,)?) => {{
        $crate::Poly::from(vec![$($elems),*])
    }};
}

mod scalar;
pub use scalar::Scalar;

mod complex_util;
use complex_util::{c_neg, complex_sort_mut};
mod impl_num;
mod indexing;
mod linalg_util;
pub use indexing::Get;
mod bessel;
mod casting_util;
use casting_util::usize_to_u32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<T: Scalar>(na::DVector<Complex<T>>);

pub type Poly32 = Poly<f32>;
pub type Poly64 = Poly<f64>;

impl<T: Scalar> Poly<T> {
    pub fn new(coeffs: &[Complex<T>]) -> Self {
        Self(na::DVector::from_row_slice(coeffs)).normalize()
    }

    #[must_use]
    pub const fn from_dvector(value: na::DVector<Complex<T>>) -> Self {
        Self(value)
    }

    /// The same as `Poly::new()`
    pub fn from_complex_slice(value: &[Complex<T>]) -> Self {
        Self::new(value)
    }

    #[allow(clippy::needless_pass_by_value)]
    #[must_use]
    pub fn from_complex_vec(value: Vec<Complex<T>>) -> Self {
        Self::new(value.as_slice())
    }

    pub fn from_real_slice(value: &[T]) -> Self {
        let temp_vec: Vec<_> = value.iter().map(Complex::from).collect();
        Self::new(&temp_vec)
    }

    #[allow(clippy::needless_pass_by_value)]
    #[must_use]
    pub fn from_real_vec(value: Vec<T>) -> Self {
        Self::from(value.as_slice())
    }

    /// Monic polynomial from its complex roots.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{Zero, One};
    ///
    /// let p = Poly::from_roots(&[Complex::new(-1.0, 0.0), Complex::zero(), Complex::one()]);
    /// assert_eq!(p, Poly::new(&[Complex::zero(), Complex::new(-1.0, 0.0), Complex::zero(), Complex::one()]))
    /// ```
    #[must_use]
    pub fn from_roots(roots: &[Complex<T>]) -> Self {
        if roots.is_empty() {
            return Self::one();
        }

        let mut roots: na::DVector<Complex<T>> = na::DVector::from_column_slice(roots);
        complex_sort_mut(&mut roots);

        roots
            .map(|e| Self::line(c_neg(e), Complex::<T>::one()))
            .fold(Self::one(), |acc, x| acc * x)
            .normalize()
    }

    /// Linear function as a polynomial.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{One, Zero};
    ///
    /// assert_eq!(Poly::line(Complex::one(), Complex::new(-1.0, 0.0)).eval_point(Complex::one()), Complex::zero());
    /// ```
    pub fn line(offset: Complex<T>, slope: Complex<T>) -> Self {
        if slope.is_zero() {
            return Self::new(&[offset]);
        }
        Self::new(&[offset, slope])
    }

    /// Line between two points with complex coordinates.
    ///
    /// Note that the points are determined by two complex numbers, so they are
    /// in a four dimensional space. Leave the imaginary component as zero for lines
    /// in a 2D plane.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{One, Zero};
    ///
    /// let p1 = (Complex::new(-1.0, 0.0), Complex::new(2.0, 0.0));
    /// let p2 = (Complex::new(2.0, 0.0), Complex::new(-1.0, 0.0));
    ///
    /// assert_eq!(Poly::line_from_points(p1, p2).eval_point(Complex::one()), Complex::zero());
    /// ```
    pub fn line_from_points(p1: (Complex<T>, Complex<T>), p2: (Complex<T>, Complex<T>)) -> Self {
        let slope = (p2.1 - &p1.1) / (p2.0 - &p1.0);
        let offset = p1.1 - &slope * p1.0;
        Self::line(offset, slope)
    }

    /// Create a polynomial from a single term (coefficient + degree)
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    /// use num::One;
    ///
    /// assert_eq!(Poly::term(Complex::one(), 3), poly![0.0, 0.0, 0.0, 1.0]);
    /// ```
    pub fn term(coeff: Complex<T>, degree: u32) -> Self {
        Self::line(Complex::zero(), complex!(T::one())).pow(degree) * coeff
    }

    /// Get the nth term of the polynomial as a new polynomial
    ///
    /// Will return None if out of bounds.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    /// use num::One;
    ///
    /// let p  = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p.get_term(1).unwrap(), poly![0.0, 2.0]);
    /// ```
    #[must_use]
    pub fn get_term(&self, degree: u32) -> Option<Self> {
        if degree as usize >= self.len_raw() {
            return None;
        }
        Some(Self::term(self[degree as usize].clone(), degree))
    }

    /// Iterate over coefficients, from the least significant
    pub fn iter(&self) -> std::slice::Iter<'_, na::Complex<T>> {
        self.0.as_slice().iter()
    }

    /// Iterate over coefficients, from the least significant
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, na::Complex<T>> {
        self.0.as_mut_slice().iter_mut()
    }

    /// Get the nth [Chebyshev polynomial](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
    ///
    /// ```
    /// use rust_poly::{poly, Poly};
    ///
    /// assert_eq!(Poly::cheby(2), poly![-1.0, 0.0, 2.0]);
    /// assert_eq!(Poly::cheby(3), poly![0.0, -3.0, 0.0, 4.0]);
    /// assert_eq!(Poly::cheby(4), poly![1.0, 0.0, -8.0, 0.0, 8.0])
    /// ```
    #[must_use]
    pub fn cheby(n: usize) -> Self {
        // TODO: make the first 32-ish explicit for performance
        match n {
            0 => poly![T::one()],
            1 => poly![T::zero(), T::one()],
            2 => poly![-T::one(), T::zero(), T::two()],
            3 => poly![T::zero(), -T::three(), T::zero(), T::four()],
            4 => poly![T::one(), T::zero(), -T::eight(), T::zero(), T::eight()],
            _ => poly![T::zero(), T::two()] * Self::cheby(n - 1) - Self::cheby(n - 2),
        }
    }

    /// Get the nth [Bessel polynomial](https://en.wikipedia.org/wiki/Bessel_polynomials)
    #[must_use]
    pub fn bessel(n: usize) -> Option<Self> {
        let mut poly = poly![];
        for k in 0..=n {
            let c = T::from_f64(coeff(n, k))?;
            let term = Self::term(complex!(c), usize_to_u32(k));
            poly = poly + term;
        }
        Some(poly)
    }

    #[must_use]
    pub fn reverse_bessel(n: usize) -> Option<Self> {
        let p = Self::bessel(n)?;
        let v: Vec<_> = p.iter().cloned().rev().collect();
        Some(Self::from_complex_vec(v))
    }

    fn len_raw(&self) -> usize {
        self.0.len()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert!(self.is_normalized());
        self.len_raw()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn is_normalized(&self) -> bool {
        let n = self.len_raw();
        if n == 0 {
            return true;
        }
        !self.0.index(n - 1).is_zero()
    }

    fn normalize(self) -> Self {
        if self.is_normalized() {
            return self;
        }
        let mut end = self.len_raw();
        loop {
            if end == 0 {
                return Self::zero();
            }
            if !self.0.as_slice()[end - 1].is_zero() {
                break;
            }
            end -= 1;
        }
        let ret = Self(na::DVector::from_column_slice(&self.0.as_slice()[0..end]));

        // post-condition: polynomial is now normalized
        debug_assert!(ret.is_normalized());
        ret
    }

    /// Evaluate the polynomial at a single value of `x`.
    ///
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    ///
    /// let p = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0)]);
    /// let x = Complex::new(1.0, 0.0);
    /// assert_eq!(p.eval_point(x), Complex::new(6.0, 0.0));
    /// ```
    pub fn eval_point(&self, x: Complex<T>) -> Complex<T> {
        self.eval(&na::DMatrix::<_>::from_row_slice(1, 1, &[x]))[0].clone()
    }

    /// Evaluate the polynomial for each entry of a matrix.
    #[must_use]
    pub fn eval(&self, x: &na::DMatrix<Complex<T>>) -> na::DMatrix<Complex<T>> {
        let mut c0: na::DMatrix<_> = na::DMatrix::<_>::from_element(
            x.nrows(),
            x.ncols(),
            self.0[self.len_raw() - 1].clone(),
        );
        for i in 2..=self.len_raw() {
            c0 *= x.clone();
            c0.apply(|c| *c = (*c).clone() + &self.0[self.len_raw() - i]);
        }
        c0
    }

    /// Raises a polynomial to an integer power.
    ///
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// assert_eq!(poly![1.0, 2.0, 3.0].pow(2), poly![1.0, 4.0, 10.0, 12.0, 9.0]);
    /// ```
    #[must_use]
    pub fn pow(self, pow: u32) -> Self {
        self.pow_usize(pow as usize)
    }

    #[must_use]
    pub fn pow_usize(self, pow: usize) -> Self {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        if pow == 0 {
            return Self::one();
        }

        if pow == 1 {
            return self;
        }

        // TODO: divide and conquer with powers of 2
        let mut res = self.clone();
        for _ in 2..=pow {
            res = res * self.clone();
        }
        res.normalize()
    }

    fn companion(&self) -> na::DMatrix<Complex<T>> {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        // pre-condition: poly has degree 1 or more
        assert!(
            self.len_raw() >= 2,
            "polynomials of degree 0 or less do not have a companion matrix"
        );

        if self.len_raw() == 2 {
            return na::DMatrix::from_row_slice(
                1,
                1,
                &[c_neg(self.0[0].clone()) / self.0[1].clone()],
            );
        }

        let n = self.len_raw() - 1;
        let mut mat: na::DMatrix<Complex<T>> = na::DMatrix::<Complex<T>>::zeros(n, n);

        // fill sub-diagonal with 1
        mat.view_mut((1, 0), (n - 1, n - 1))
            .fill_diagonal(Complex::<T>::one());

        // fill the rightmost column with the coefficients of the associated
        // monic polynomial
        let monic = self
            .0
            .view((0, 0), (n, 1))
            .map(|x| c_neg(x) / self.0[n].clone());
        for i in 0..n {
            mat.column_mut(n - 1)[i] = monic[i].clone();
        }
        mat
    }

    /// Find the roots of a polynomial numerically.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// let p: Poly<f64> = poly![-6.0, 11.0, -6.0, 1.0];
    /// let expected_roots = &[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0)];
    /// let temp = p.roots();    // Rust is really annoying sometimes...
    /// let calculated_roots = temp.as_slice();
    ///
    /// // assert almost equal
    /// assert!((expected_roots[0] - calculated_roots[0]).re.abs() < 0.000001);
    /// assert!((expected_roots[1] - calculated_roots[1]).re.abs() < 0.000001);
    /// assert!((expected_roots[2] - calculated_roots[2]).re.abs() < 0.000001);
    /// ```
    #[allow(clippy::missing_panics_doc)]
    #[must_use]
    pub fn roots(&self) -> Vec<Complex<T>> {
        // invariant: polynomial is normalized
        debug_assert!(self.is_normalized());

        if self.len_raw() < 2 {
            return vec![];
        }

        if self.len_raw() == 2 {
            return vec![c_neg(self.0[0].clone()) / self.0[1].clone()];
        }

        // rotated companion matrix reduces error
        let mut comp = self.companion();
        let n = comp.shape().0;
        for i in 0..n / 2 {
            comp.swap_rows(i, n - i - 1);
            comp.swap_columns(i, n - i - 1);
        }

        let mut r: na::DVector<Complex<T>> = comp.eigenvalues().expect("infallible");
        complex_sort_mut(&mut r);
        r.as_slice().to_vec()
    }

    /// Compose two polynomials, returning a new polynomial.
    ///
    /// Substitute the given polynomial `x` into `self` and expand the
    /// result into a new polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_poly::{Poly, poly};
    /// use num::{One, Complex};
    ///
    /// let f = poly![1.0, 2.0];
    /// let g = Poly::one();
    ///
    /// assert_eq!(f.clone().compose(g), f);
    #[must_use]
    pub fn compose(self, x: Self) -> Self {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(x.is_normalized());

        // TODO begin: are these checks actually making things faster?
        if self.is_zero() || x.is_zero() {
            return Self::zero();
        }

        if self.is_one() {
            return x;
        }

        if x.is_one() {
            return self;
        }
        // end

        (0..self.len_raw())
            .map(|i| Self::new(&[self.0[i].clone()]) * x.clone().pow_usize(i))
            .sum()
    }

    /// Returns true if every coefficient in the polynomial is smaller than the
    /// tolerance (using complex norm).
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{Poly, poly};
    ///
    /// assert!(poly![0.01, -0.01].almost_zero(&0.1));
    /// ```
    #[must_use]
    pub fn almost_zero(&self, tolerance: &T) -> bool {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());

        self.as_slice().iter().all(|c| c.norm() <= *tolerance)
    }

    /// Calculate the quotient and remainder uwing long division. More efficient than
    /// calculating them separately.
    ///
    /// # Panics
    /// Panics if a division by zero is attempted
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{Poly, poly};
    /// use num::{Complex, One};
    ///
    /// let c1 = poly![1.0, 2.0, 3.0];
    /// let c2 = poly![3.0, 2.0, 1.0];
    /// let expected1 = (poly![3.0], poly![-8.0, -4.0]);
    /// assert_eq!(c1.clone().div_rem(&c2).unwrap(), expected1);
    /// ```
    #[allow(clippy::cast_sign_loss)]
    #[allow(clippy::cast_possible_wrap)]
    #[must_use]
    pub fn div_rem(self, rhs: &Self) -> Option<(Self, Self)> {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(rhs.is_normalized());

        // pre-condition: don't divide by zero
        if rhs.is_zero() {
            // bail!("Attempted to divide a polynomial by zero");
            return None;
        }

        let lhs_len = self.len_raw();
        let rhs_len = self.len_raw();
        if lhs_len < rhs_len {
            return Some((Self::zero(), self));
        }
        if rhs_len == 1 {
            return Some((
                // TODO: should use checked operations
                Self(self.0 / rhs.0[rhs.len_raw() - 1].clone()),
                Self::zero(),
            ));
        }
        let len_delta = lhs_len - rhs_len;
        let scale = rhs.0[rhs.len_raw() - 1].clone();
        let rhs: na::DVector<_> = rhs
            .0
            .view_range(0..rhs.len_raw() - 1, 0..1)
            // HACK: this shouldn't be necessary, but nalgebra turns DVector into
            //       DMatrix when making a view, and needs to be politely reminded
            //       that this is a column vector.
            .column(0)
            .into();
        // TODO: useless clone of scale, it should be borrowed, but dvector does
        //       not implement Div<&_>
        // TODO: should use checked operations
        let rhs: na::DVector<_> = rhs / scale.clone();
        let mut lhs: na::DVector<_> = self.0.clone();
        let mut i = len_delta as isize;
        let mut j = (lhs_len - 1) as isize;
        while i >= 0 {
            lhs.view_range_mut(i as usize..j as usize, 0..1)
                .iter_mut()
                .zip((rhs.clone() * self.0[j as usize].clone()).iter())
                .for_each(|p| *p.0 -= p.1);
            i -= 1;
            j -= 1;
        }
        Some((
            Self(
                (lhs.view_range(j as usize + 1..lhs.len(), 0..1) / scale)
                    .column(0)
                    .into(),
            )
            .normalize(),
            Self(lhs.view_range(..(j + 1) as usize, 0..1).column(0).into()).normalize(),
        ))
    }

    fn checked_div_impl(self, rhs: &Self) -> Option<Self> {
        Some(self.div_rem(rhs)?.0)
    }

    fn checked_rem_impl(self, rhs: &Self) -> Option<Self> {
        Some(self.div_rem(rhs)?.1)
    }

    #[must_use]
    pub fn as_slice(&self) -> &[Complex<T>] {
        self.0.as_slice()
    }

    pub fn as_mut_slice(&mut self) -> &mut [Complex<T>] {
        self.0.as_mut_slice()
    }

    #[must_use]
    pub fn as_ptr(&self) -> *const Complex<T> {
        self.0.as_ptr()
    }

    pub fn as_mut_ptr(&mut self) -> *mut Complex<T> {
        self.0.as_mut_ptr()
    }

    #[must_use]
    pub fn as_view(&self) -> na::DMatrixView<Complex<T>> {
        self.0.as_view()
    }

    pub fn as_view_mut(&mut self) -> na::DMatrixViewMut<Complex<T>> {
        self.0.as_view_mut()
    }

    #[must_use]
    pub fn to_vec(&self) -> Vec<Complex<T>> {
        Vec::from(self.as_slice())
    }

    #[must_use]
    pub fn to_dvector(self) -> na::DVector<Complex<T>> {
        self.0
    }
}

impl<T: Scalar> Index<usize> for Poly<T> {
    type Output = Complex<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T: Scalar> From<na::DVector<Complex<T>>> for Poly<T> {
    fn from(value: na::DVector<Complex<T>>) -> Self {
        Self::from_dvector(value)
    }
}

impl<T: Scalar> From<&[Complex<T>]> for Poly<T> {
    fn from(value: &[Complex<T>]) -> Self {
        Self::from_complex_slice(value)
    }
}

impl<T: Scalar> From<Vec<Complex<T>>> for Poly<T> {
    fn from(value: Vec<Complex<T>>) -> Self {
        Self::from_complex_vec(value)
    }
}

impl<T: Scalar> From<&[T]> for Poly<T> {
    fn from(value: &[T]) -> Self {
        Self::from_real_slice(value)
    }
}

impl<T: Scalar> From<Vec<T>> for Poly<T> {
    fn from(value: Vec<T>) -> Self {
        Self::from_real_vec(value)
    }
}

impl<T: Scalar> From<Poly<T>> for *const Complex<T> {
    fn from(val: Poly<T>) -> Self {
        val.as_ptr()
    }
}

impl<T: Scalar> From<Poly<T>> for *mut Complex<T> {
    fn from(mut val: Poly<T>) -> Self {
        val.as_mut_ptr()
    }
}

impl<T: Scalar> From<Poly<T>> for Vec<Complex<T>> {
    fn from(val: Poly<T>) -> Self {
        val.to_vec()
    }
}

impl<T: Scalar> From<Poly<T>> for na::DVector<Complex<T>> {
    fn from(val: Poly<T>) -> Self {
        val.to_dvector()
    }
}

impl<'a, T: Scalar> IntoIterator for &'a Poly<T> {
    type IntoIter = std::slice::Iter<'a, na::Complex<T>>;
    type Item = &'a na::Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T: Scalar> IntoIterator for &'a mut Poly<T> {
    type IntoIter = std::slice::IterMut<'a, na::Complex<T>>;
    type Item = &'a mut na::Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn macro_complex() {
        assert_eq!(complex!(), Complex::<f64>::zero());
        assert_eq!(complex!(1.0, 2.0), Complex::<f64>::new(1.0, 2.0));
    }

    #[test]
    fn macro_poly() {
        assert_eq!(poly!(), Poly::<f64>::zero());
        assert_eq!(poly!(1.0), Poly::<f64>::one());
        assert_eq!(
            poly!(1.0, 2.0, 3.0),
            Poly::<f64>::new(&[
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
            ])
        );
        assert_eq!(poly!((1.0, 0.0)), Poly::<f64>::one());
        assert_eq!(
            poly!((1.0, 1.0), (2.0, 2.0), (3.0, 3.0)),
            Poly::<f64>::new(&[
                Complex::new(1.0, 1.0),
                Complex::new(2.0, 2.0),
                Complex::new(3.0, 3.0)
            ])
        );
        assert_eq!(
            poly!(2.0; 3),
            Poly::<f64>::new(&[
                Complex::new(2.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(2.0, 0.0)
            ])
        );
        assert_eq!(
            poly!((1.0, -1.0); 3),
            Poly::<f64>::new(&[
                Complex::new(1.0, -1.0),
                Complex::new(1.0, -1.0),
                Complex::new(1.0, -1.0)
            ])
        );
    }

    #[test]
    fn poly_new() {
        // trivial, but here for completeness
        Poly::new(&[Complex::new(2.0, -2.0)]);
    }

    #[test]
    fn poly_from_complex_slice() {
        let p = Poly::from_complex_slice(&[Complex::new(1.0, 2.0), Complex::new(3.0, 4.0)]);
        let e = poly!((1.0, 2.0), (3.0, 4.0));
        assert_eq!(p, e);
    }

    // TODO: test the rest of the "boring" functions

    #[test]
    fn poly_line() {
        let p = Poly::<f64>::line(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(2.0, 0.0));
        let e = poly!(1.0, 2.0);
        assert_eq!(p, e);
    }

    #[test]
    fn poly_term() {
        let p = Poly64::term(complex!(2.0), 2);
        let e = poly!(0.0, 0.0, 2.0);
        assert_eq!(p, e);
    }

    #[test]
    fn poly_bessel() {
        assert_eq!(Poly64::bessel(0).unwrap(), poly![1.0]);
        assert_eq!(Poly64::bessel(1).unwrap(), poly![1.0, 1.0]);
        assert_eq!(Poly64::bessel(2).unwrap(), poly![1.0, 3.0, 3.0]);
        assert_eq!(Poly64::bessel(3).unwrap(), poly![1.0, 6.0, 15.0, 15.0]);
        assert_eq!(
            Poly64::bessel(4).unwrap(),
            poly![1.0, 10.0, 45.0, 105.0, 105.0]
        );
        assert_eq!(
            Poly64::bessel(5).unwrap(),
            poly![1.0, 15.0, 105.0, 420.0, 945.0, 945.0]
        );
        assert_eq!(
            Poly64::bessel(6).unwrap(),
            poly![1.0, 21.0, 210.0, 1260.0, 4725.0, 10395.0, 10395.0]
        );
    }

    #[test]
    fn poly_reverse_bessel() {
        assert_eq!(Poly64::reverse_bessel(0).unwrap(), poly![1.0]);
        assert_eq!(Poly64::reverse_bessel(1).unwrap(), poly![1.0, 1.0]);
        assert_eq!(Poly64::reverse_bessel(2).unwrap(), poly![3.0, 3.0, 1.0]);
        assert_eq!(
            Poly64::reverse_bessel(3).unwrap(),
            poly![15.0, 15.0, 6.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(4).unwrap(),
            poly![105.0, 105.0, 45.0, 10.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(5).unwrap(),
            poly![945.0, 945.0, 420.0, 105.0, 15.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(6).unwrap(),
            poly![10395.0, 10395.0, 4725.0, 1260.0, 210.0, 21.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(7).unwrap(),
            poly![135135.0, 135135.0, 62370.0, 17325.0, 3150.0, 378.0, 28.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(8).unwrap(),
            poly![2027025.0, 2027025.0, 945945.0, 270270.0, 51975.0, 6930.0, 630.0, 36.0, 1.0]
        );
        assert_eq!(
            Poly64::reverse_bessel(9).unwrap(),
            poly![
                34459425.0, 34459425.0, 16216200.0, 4729725.0, 945945.0, 135135.0, 13860.0, 990.0,
                45.0, 1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(10).unwrap(),
            poly![
                654729075.0,
                654729075.0,
                310134825.0,
                91891800.0,
                18918900.0,
                2837835.0,
                315315.0,
                25740.0,
                1485.0,
                55.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(11).unwrap(),
            poly![
                13749310575.0,
                13749310575.0,
                6547290750.0,
                1964187225.0,
                413513100.0,
                64324260.0,
                7567560.0,
                675675.0,
                45045.0,
                2145.0,
                66.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(12).unwrap(),
            poly![
                316234143225.0,
                316234143225.0,
                151242416325.0,
                45831035250.0,
                9820936125.0,
                1571349780.0,
                192972780.0,
                18378360.0,
                1351350.0,
                75075.0,
                3003.0,
                78.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(13).unwrap(),
            poly![
                7905853580625.0,
                7905853580625.0,
                3794809718700.0,
                1159525191825.0,
                252070693875.0,
                41247931725.0,
                5237832600.0,
                523783260.0,
                41351310.0,
                2552550.0,
                120120.0,
                4095.0,
                91.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(14).unwrap(),
            poly![
                213458046676875.0,
                213458046676875.0,
                102776096548125.0,
                31623414322500.0,
                6957151150950.0,
                1159525191825.0,
                151242416325.0,
                15713497800.0,
                1309458150.0,
                87297210.0,
                4594590.0,
                185640.0,
                5460.0,
                105.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(15).unwrap(),
            poly![
                6190283353629375.0,
                6190283353629375.0,
                2988412653476250.0,
                924984868933125.0,
                205552193096250.0,
                34785755754750.0,
                4638100767300.0,
                496939367925.0,
                43212118950.0,
                3055402350.0,
                174594420.0,
                7936110.0,
                278460.0,
                7140.0,
                120.0,
                1.0
            ]
        );
        assert_eq!(
            Poly64::reverse_bessel(16).unwrap(),
            poly![
                191898783962510625.0,
                191898783962510625.0,
                92854250304440625.0,
                28887988983603750.0,
                6474894082531875.0,
                1109981842719750.0,
                150738274937250.0,
                16564645597500.0,
                1490818103775.0,
                110430970650.0,
                6721885170.0,
                333316620.0,
                13226850.0,
                406980.0,
                9180.0,
                136.0,
                1.0
            ]
        );
    }
}
