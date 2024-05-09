use std::fmt::Display;

use nalgebra::RealField;
use num::{Complex, Float, One, Zero};

use crate::{
    util::complex::{c_neg, complex_fmt, complex_sort_mut},
    Scalar, ScalarOps,
};

mod calculus;
mod conversions;
mod impl_num;
mod indexing;
mod internals;
mod special_funcs;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<T>(pub(crate) na::DVector<Complex<T>>);

impl<T: Scalar> Poly<T> {
    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p.shift_up(2), poly![0.0, 0.0, 1.0, 2.0, 3.0]);
    /// ```
    #[must_use]
    pub fn shift_up(&self, n: usize) -> Self {
        let mut v = vec![Complex::<T>::zero(); n];
        v.extend_from_slice(self.as_slice());
        Self::from_complex_vec(v)
    }

    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0, 4.0];
    /// assert_eq!(p.shift_down(2), poly![3.0, 4.0]);
    /// ```
    #[must_use]
    pub fn shift_down(&self, n: usize) -> Self {
        Self::from_complex_slice(&self.as_slice()[n..])
    }
}

impl<T: Scalar> Poly<T> {
    pub fn new(coeffs: &[Complex<T>]) -> Self {
        Self(na::DVector::from_row_slice(coeffs)).normalize()
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

    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert!(self.is_normalized());
        self.len_raw()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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

    /// Iterate over the terms of a polynomial
    ///
    /// ```
    /// # use rust_poly::{poly, Poly};
    ///
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p, p.terms().sum::<Poly<_>>());
    /// ```
    pub fn terms(&self) -> std::iter::Map<std::ops::Range<usize>, impl FnMut(usize) -> Self + '_> {
        debug_assert!(self.is_normalized());
        (0..self.len_raw()).map(|i| {
            self.get_term(
                i.try_into()
                    .expect("degrees above u32::MAX are not supported"),
            )
            .expect("terms are within range len_raw, this should never fail")
        })
    }
}

impl<T: Scalar + PartialOrd> Poly<T> {
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
}

impl<T: ScalarOps> Poly<T> {
    /// Evaluate the polynomial for each entry of a matrix.
    #[must_use]
    pub fn eval(&self, x: &na::DMatrix<Complex<T>>) -> na::DMatrix<Complex<T>> {
        let mut c0: na::DMatrix<_> = na::DMatrix::<_>::from_element(
            x.nrows(),
            x.ncols(),
            self.0[self.len_raw() - 1].clone(),
        );
        for i in 2..=self.len_raw() {
            c0 = c0.clone() * x.clone();
            c0.apply(|c| *c = (*c).clone() + &self.0[self.len_raw() - i]);
        }
        c0
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
}

impl<T: ScalarOps + PartialOrd> Poly<T> {
    /// Translate along x-axis (or x-plane) and y-axis (or y-plane).
    ///
    /// Using complex coordinates means you'll effectively be translating in
    /// 4D space.
    pub fn translate(mut self, x: Complex<T>, y: Complex<T>) -> Self {
        self = self.compose(Self::from_complex_slice(&[c_neg(x), Complex::<T>::one()]));
        self.0[0] += y;
        self
    }
}

impl<T: Scalar + RealField> Poly<T> {
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
}

impl<T: Scalar + Float> Poly<T> {
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
}

impl<T: Scalar + Display + PartialOrd> Display for Poly<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.iter().enumerate();
        if let Some((_, c)) = iter.next() {
            write!(f, "{}", complex_fmt(c))?;
        } else {
            return Ok(());
        }
        for (i, c) in iter {
            write!(f, " + {}*x^{}", complex_fmt(c), i)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn translate() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(
            p.translate(complex!(1.0), complex!(2.0)),
            poly![4.0, -4.0, 3.0]
        );
    }

    #[test]
    fn display() {
        let p = poly![(2.0, 0.0), (4.5, 0.0), (5.0, 1.0), (6.0, 1.5), (7.0, 2.0)];
        assert_eq!(
            p.to_string(),
            "2 + 4.5*x^1 + (5+i)*x^2 + (6+i1.5)*x^3 + (7+i2)*x^4".to_string()
        );
    }
}
