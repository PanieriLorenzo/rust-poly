#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]
use ndarray::{array, s, Array, Array1, ArrayView1, ArrayViewMut1, Axis, Dimension};
use num_traits::Zero;
use std::ops::{Add, Div, Mul, Rem, Sub};

mod scalar;
pub use scalar::Scalar;

/// polynomial as a list of coefficients of terms of descending degree
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Poly<T: Scalar>(Array1<T>);

impl<T: Scalar> Poly<T> {
    /// Create a new polynomial from a 1D array of coefficients
    #[must_use]
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
    #[allow(clippy::missing_panics_doc)]
    pub fn term(coeff: T, degree: usize) -> Self {
        let zeros: Array1<T> = Array1::<T>::zeros([degree]);
        let mut term: Array1<T> = array![coeff];
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
    /// let y = p.eval(&x);
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
    /// let y = p.eval(&x);
    /// assert_eq!(y, array![Complex64::new(2.0, 3.0), Complex64::new(-2.0, 1.0)]);
    /// ```
    pub fn eval<D: Dimension>(&self, x: &Array<T, D>) -> Array<T, D> {
        let mut y: Array<T, D> = Array::<T, D>::zeros(x.raw_dim());
        for pv in &self.0 {
            y = y * x.clone() + pv.clone();
        }
        y
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
    ///
    /// let p1 = Poly::new(array![3.0, 5.0, 2.0]);
    /// let p2 = Poly::new(array![2.0, 1.0]);
    /// let (q, r) = p1.div_rem(&p2);
    /// assert_eq!(q, Poly::new(array![1.5, 1.75]));
    /// assert_eq!(r, Poly::new(array![0.25]));
    /// ```
    ///
    /// Divide two complex polynomials
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    /// use num_complex::Complex64;
    ///
    /// let p1 = Poly::term(Complex64::new(1.0, 1.0), 2);
    /// let p2 = Poly::term(Complex64::new(1.0, -1.0), 0);
    /// let (q, r) = p1.div_rem(&p2);
    /// assert_eq!(q, Poly::term(Complex64::new(0.0, 1.0), 2));
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
    ///
    /// let p1 = Poly::new(array![1.0, 2.0, 3.0]);
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

        let num: Array1<T> = self.0;
        let den: Array1<T> = rhs.0.clone();

        // cannot underflow because we have ensured that len is at least 1
        let num_deg = num.len() - 1;
        let den_deg = den.len() - 1;

        let scale = T::one() / den[0].clone();
        let mut quot: Array1<T> = Array1::zeros((num_deg - den_deg + 1).max(1));
        let mut rem: Array1<T> = num;
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

    #[must_use]
    pub fn as_ndarray(&self) -> ArrayView1<T> {
        self.0.view()
    }

    pub fn as_ndarray_mut(&mut self) -> ArrayViewMut1<T> {
        self.0.view_mut()
    }

    #[must_use]
    pub fn to_vec(&self) -> Vec<T> {
        self.0.to_vec()
    }
}

impl<T: Scalar> Zero for Poly<T> {
    fn zero() -> Self {
        Self::new(array![])
    }

    fn is_zero(&self) -> bool {
        self.is_empty()
    }
}

impl<T: Scalar> Add for Poly<T> {
    type Output = Self;

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
        // will only wrap with polynomials so large they don't fit in memory
        #[allow(clippy::cast_possible_wrap)]
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        let abs_delta = len_delta.unsigned_abs();
        match len_delta {
            0 => Self(self.0 + rhs.0),
            1.. => {
                let mut rhs_new: Array1<T> = Array1::<T>::zeros([abs_delta]);
                rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
                Self(self.0 + rhs_new)
            }
            _ => {
                let mut lhs: Array1<T> = Array1::<T>::zeros([abs_delta]);
                lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
                Self(lhs + rhs.0)
            }
        }
        .trim_zeros()
        // (if len_delta == 0 {
        //     Self(self.0 + rhs.0)
        // } else if len_delta < 0 {
        //     let mut lhs: Array1<T> = Array1::<T>::zeros([len_delta.unsigned_abs()]);
        //     lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
        //     Self(lhs + rhs.0)
        // } else {
        //     // guaranteed to be positive because we checked it in the conditional
        //     #[allow(clippy::cast_sign_loss)]
        //     let mut rhs_new: Array1<T> = Array1::<T>::zeros([len_delta as usize]);
        //     rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
        //     Self(self.0 + rhs_new)
        // })
        // .trim_zeros()
    }
}

impl<T: Scalar> Sub for Poly<T> {
    type Output = Self;

    /// Subtract one polynomial from another
    ///
    /// ## Examples
    /// Subtract polynomials of various lengths
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let p1 = Poly::new(array![1.0, 0.0]);
    /// let p2 = Poly::new(array![1.0]);
    /// assert_eq!(p1.clone() - p1.clone(), Poly::new(array![]));
    /// assert_eq!(p2.clone() - p1.clone(), Poly::new(array![-1.0, 1.0]));
    /// assert_eq!(p1 - p2, Poly::new(array![1.0, -1.0]));
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        // will only wrap with polynomials so large they don't fit in memory
        #[allow(clippy::cast_possible_wrap)]
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        let abs_delta = len_delta.unsigned_abs();
        match len_delta {
            0 => Self(self.0 - rhs.0),
            1.. => {
                let mut rhs_new: Array1<T> = Array1::<T>::zeros([abs_delta]);
                rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
                Self(self.0 - rhs_new)
            }
            _ => {
                let mut lhs: Array1<T> = Array1::<T>::zeros([abs_delta]);
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
    /// ## Examples
    ///
    /// Convolve two polynomials
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let p1 = Poly::new(array![1.0, 2.0, 3.0]);
    /// let p2 = Poly::new(array![9.0, 5.0, 1.0]);
    /// let prod = p1 * p2;
    /// assert_eq!(prod, Poly::new(array![9.0, 23.0, 38.0, 17.0, 3.0]));
    /// ```
    ///
    /// Scalar multiplication
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let p1 = Poly::term(3, 0);
    /// let p2 = Poly::new(array![1, 1]);
    /// let prod1 = p1.clone() * p2.clone();
    /// let prod2 = p2.clone() * p1.clone();
    /// assert_eq!(prod1, Poly::new(array![3, 3]));
    /// assert_eq!(prod2, Poly::new(array![3, 3]));
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
fn convolve_1d<T: Scalar>(input: ArrayView1<T>, kernel: ArrayView1<T>) -> Array1<T> {
    let input_len = input.len();
    let kernel_len = kernel.len();
    let output_len = input_len + kernel_len - 1;

    let mut output: Array1<T> = Array1::<T>::zeros([output_len]);

    for i in 0..output_len {
        let mut sum = T::zero();
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
