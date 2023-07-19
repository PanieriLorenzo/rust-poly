// Implementation of traits related to numeric operations, operators and number theory

use ndarray::{array, s, Array, Array1, ArrayView1, ArrayViewMut1, Axis, Dimension};
use num_traits::{One, Zero};
use std::ops::{Add, Div, Mul, Rem, Sub};

use crate::{Poly, Scalar, A, AC, C};

impl<T: Scalar> Zero for Poly<T> {
    fn zero() -> Self {
        Self::new(array![])
    }

    fn is_zero(&self) -> bool {
        self.is_empty()
    }
}

impl<T: Scalar> One for Poly<T> {
    fn one() -> Self {
        Poly::term(C::<T>::one(), 0)
    }
}

impl<T: Scalar> Add for Poly<T> {
    type Output = Self;

    /// Add toghether two polynomials
    ///
    /// # Examples
    /// Add polynomials of various lengths
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::new(array![Complex::from(1.0), Complex::from(0.0)]);
    /// let p2 = Poly::new(array![Complex::from(1.0)]);
    /// assert_eq!(p1.clone() + p1.clone(), Poly::new(array![Complex::from(2.0), Complex::from(0.0)]));
    /// assert_eq!(p2.clone() + p1.clone(), Poly::new(array![Complex::from(1.0), Complex::from(1.0)]));
    /// assert_eq!(p1 + p2, Poly::new(array![Complex::from(1.0), Complex::from(1.0)]));
    /// ```
    ///
    /// Add three terms to form a polynomial
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let t1 = Poly::term(Complex::from(1), 0);
    /// let t2 = Poly::term(Complex::from(2), 1);
    /// let t3 = Poly::term(Complex::from(3), 2);
    /// let sum = t1 + t2 + t3;
    /// assert_eq!(sum, Poly::new(array![Complex::from(3), Complex::from(2), Complex::from(1)]));
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        // will only wrap with polynomials so large they don't fit in memory
        #[allow(clippy::cast_possible_wrap)]
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        let abs_delta = len_delta.unsigned_abs();
        match len_delta {
            0 => Self(self.0 + rhs.0),
            1.. => {
                let mut rhs_new: AC<T> = AC::<T>::zeros([abs_delta]);
                rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
                Self(self.0 + rhs_new)
            }
            _ => {
                let mut lhs: AC<T> = AC::<T>::zeros([abs_delta]);
                lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
                Self(lhs + rhs.0)
            }
        }
        .trim_zeros()
    }
}

impl<T: Scalar> Sub for Poly<T> {
    type Output = Self;

    /// Subtract one polynomial from another
    ///
    /// # Examples
    /// Subtract polynomials of various lengths
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::new(array![Complex::from(1.0), Complex::from(0.0)]);
    /// let p2 = Poly::new(array![Complex::from(1.0)]);
    /// assert_eq!(p1.clone() - p1.clone(), Poly::new(array![]));
    /// assert_eq!(p2.clone() - p1.clone(), Poly::new(array![Complex::from(-1.0), Complex::from(1.0)]));
    /// assert_eq!(p1 - p2, Poly::new(array![Complex::from(1.0), Complex::from(-1.0)]));
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        // will only wrap with polynomials so large they don't fit in memory
        #[allow(clippy::cast_possible_wrap)]
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        let abs_delta = len_delta.unsigned_abs();
        match len_delta {
            0 => Self(self.0 - rhs.0),
            1.. => {
                let mut rhs_new: AC<T> = AC::<T>::zeros([abs_delta]);
                rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
                Self(self.0 - rhs_new)
            }
            _ => {
                let mut lhs: AC<T> = AC::<T>::zeros([abs_delta]);
                lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
                Self(lhs - rhs.0)
            }
        }
        .trim_zeros()
    }
}

impl<T: Scalar> Mul for Poly<T> {
    type Output = Self;

    /// Multiplies two polynomials together
    ///
    /// # Examples
    ///
    /// Convolve two polynomials
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::new(array![Complex::from(1.0), Complex::from(2.0), Complex::from(3.0)]);
    /// let p2 = Poly::new(array![Complex::from(9.0), Complex::from(5.0), Complex::from(1.0)]);
    /// let prod = p1 * p2;
    /// assert_eq!(prod, Poly::new(array![Complex::from(9.0), Complex::from(23.0), Complex::from(38.0), Complex::from(17.0), Complex::from(3.0)]));
    /// ```
    ///
    /// Scalar multiplication
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex;
    ///
    /// let p1 = Poly::term(Complex::from(3), 0);
    /// let p2 = Poly::new(array![Complex::from(1), Complex::from(1)]);
    /// let prod1 = p1.clone() * p2.clone();
    /// let prod2 = p2.clone() * p1.clone();
    /// assert_eq!(prod1, Poly::new(array![Complex::from(3), Complex::from(3)]));
    /// assert_eq!(prod2, Poly::new(array![Complex::from(3), Complex::from(3)]));
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        Self(convolve_1d(self.0.view(), rhs.0.view())).trim_zeros()
    }
}

impl<T: Scalar> Div<&Self> for Poly<T> {
    type Output = Self;

    /// Computes the quotient of two polynomials, truncating the remainder.
    ///
    /// See also `Poly::div_rem()`.
    fn div(self, rhs: &Self) -> Self::Output {
        self.div_rem(rhs).0
    }
}

impl<T: Scalar> Rem<&Self> for Poly<T> {
    type Output = Self;

    /// Computes the remainder of the division of two polynomials.
    ///
    /// See also `Poly::div_rem()`.
    fn rem(self, rhs: &Self) -> Self::Output {
        self.div_rem(rhs).1
    }
}

// TODO: this is a slightly modified version of a ChatGPT 3.5 answer,
//     integrate it better by placing it inside `Poly::mul` and changing
//     the name of variables.
fn convolve_1d<T: Scalar>(input: ArrayView1<C<T>>, kernel: ArrayView1<C<T>>) -> AC<T> {
    let input_len = input.len();
    let kernel_len = kernel.len();
    let output_len = input_len + kernel_len - 1;

    let mut output: AC<T> = AC::<T>::zeros([output_len]);

    for i in 0..output_len {
        let mut sum = C::<T>::zero();
        for j in 0..kernel_len {
            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            let k = i as isize - j as isize;

            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            // k is guaranteed to be positive by the conditional
            #[allow(clippy::cast_sign_loss)]
            if k >= 0 && k < input_len as isize {
                sum = sum.clone() + input[k as usize].clone() * kernel[j].clone();
            }
        }
        output[i] = sum;
    }
    output
}
