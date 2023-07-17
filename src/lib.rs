use ndarray::{array, s, Array, Array1, Axis, Dimension, ScalarOperand};
use num_traits::{Num, Zero};
use std::ops::Add;

mod scalar;
pub use scalar::Scalar;

/// polynomial as a list of coefficients of terms of descending degree
#[derive(Clone, PartialEq, Debug)]
pub struct Poly<T: Scalar>(Array1<T>);

impl<T: Scalar> Poly<T> {
    /// Create a new polynomial from a 1D array of coefficients
    pub fn new(coeffs: Array1<T>) -> Self {
        Self(coeffs).trim_zeros()
    }

    /// Creates a polynomial with a single term of degree `n`.
    ///
    /// ## Examples
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let t1 = Poly::term(1i32, 0);
    /// let t2 = Poly::term(2i32, 1);
    /// let t3 = Poly::term(3i32, 2);
    /// assert_eq!(t1, Poly::new(array![1]));
    /// assert_eq!(t2, Poly::new(array![2, 0]));
    /// assert_eq!(t3, Poly::new(array![3, 0, 0]));
    /// ```
    pub fn term(coeff: T, degree: usize) -> Self {
        let zeros: Array1<T> = Array1::<T>::zeros([degree]);
        let mut term: Array1<T> = array![coeff];
        term.append(Axis(0), zeros.view()).unwrap(); // TODO
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
    pub fn len(&self) -> usize {
        self.trim_zeros().raw_len()
    }

    /// Length of the polynomial, without trimming zeros
    fn raw_len(&self) -> usize {
        self.0.len()
    }

    /// Evaluate a polynomial at a specific input value `x`. This may be an
    /// ndarray of any dimension
    ///
    /// ## Examples
    ///
    /// Evaluate a real polynomial at real points
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// // x^2 + 2x + 1
    /// let p = Poly::new(array![1, 2, 1]);
    /// let x = array![-1, 0, 1];
    /// let y = p.eval(x);
    /// assert_eq!(y, array![0, 1, 4]);
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
    /// let y = p.eval(x);
    /// assert_eq!(y, array![Complex64::new(2.0, 3.0), Complex64::new(-2.0, 1.0)]);
    /// ```
    pub fn eval<D: Dimension>(&self, x: Array<T, D>) -> Array<T, D> {
        let mut y: Array<T, D> = Array::<T, D>::zeros(x.raw_dim());
        for pv in &self.0 {
            y = y * x.clone() + pv.clone();
        }
        y
    }
}

impl<T: Scalar> Add for Poly<T> {
    type Output = Poly<T>;

    /// Add toghether two polynomials
    ///
    /// ## Examples
    /// Add polynomials of various lengths
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let p1 = Poly::new(array![1.0, 0.0]);
    /// let p2 = Poly::new(array![1.0]);
    /// assert_eq!(p1.clone() + p1.clone(), Poly::new(array![2.0, 0.0]));
    /// assert_eq!(p2.clone() + p1.clone(), Poly::new(array![1.0, 1.0]));
    /// assert_eq!(p1 + p2, Poly::new(array![1.0, 1.0]));
    /// ```
    ///
    /// Add three terms to form a polynomial
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let t1 = Poly::term(1, 0);
    /// let t2 = Poly::term(2, 1);
    /// let t3 = Poly::term(3, 2);
    /// let sum = t1 + t2 + t3;
    /// assert_eq!(sum, Poly::new(array![3, 2, 1]));
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        if len_delta == 0 {
            Self(self.0 + rhs.0)
        } else if len_delta < 0 {
            let mut lhs: Array1<T> = Array1::<T>::zeros([len_delta.abs() as usize]);
            lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
            Self(lhs + rhs.0)
        } else {
            let mut rhs_new: Array1<T> = Array1::<T>::zeros([len_delta as usize]);
            rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
            Self(self.0 + rhs_new)
        }
    }
}
