use ndarray::{array, s, Array, Array1, ArrayView1, ArrayViewMut1, Axis, Dimension, ScalarOperand};
use num_traits::{Num, Zero};
use std::ops::{Add, Div, Mul, Rem, Sub};

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

    /// Computes the quotient and remainder of a polynomial division
    ///
    /// ## Examples
    ///
    /// Divide two real polynomials
    ///
    /// ```
    /// # use rust_poly::Poly;
    /// use ndarray::prelude::*;
    ///
    /// let p1 = Poly::new(array![3.0, 5.0, 2.0]);
    /// let p2 = Poly::new(array![2.0, 1.0]);
    /// let (q, r) = p1.div_rem(p2);
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
    /// let (q, r) = p1.div_rem(p2);
    /// assert_eq!(q, Poly::term(Complex64::new(0.0, 1.0), 2));
    /// assert_eq!(r, Poly::new(array![]));
    /// ```
    ///
    pub fn div_rem(&self, rhs: Self) -> (Self, Self) {
        let u: Array1<T> = self.0.clone() + array![T::zero()];
        let v: Array1<T> = rhs.0 + array![T::zero()];
        let m = u.len() as isize - 1;
        let n = v.len() as isize - 1;
        let scale = T::one() / v[0].clone();
        let mut q: Array1<T> = Array1::zeros((m - n + 1).max(1) as usize);
        let mut r: Array1<T> = u.clone(); // TODO: useless assignment
        for k in 0..((m - n + 1) as usize) {
            let d = scale.clone() * r[k].clone();
            q[k] = d.clone();
            r.slice_mut(s![k..(k + n as usize + 1)])
                .iter_mut()
                .zip((array![d] * v.clone()).iter())
                .for_each(|p| *p.0 = p.0.clone() - p.1.clone());
        }
        dbg!(q.clone(), r.clone());
        (Self(q), Self(r).trim_zeros())
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
        (if len_delta == 0 {
            Self(self.0 + rhs.0)
        } else if len_delta < 0 {
            let mut lhs: Array1<T> = Array1::<T>::zeros([len_delta.abs() as usize]);
            lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
            Self(lhs + rhs.0)
        } else {
            let mut rhs_new: Array1<T> = Array1::<T>::zeros([len_delta as usize]);
            rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
            Self(self.0 + rhs_new)
        })
        .trim_zeros()
    }
}

impl<T: Scalar> Sub for Poly<T> {
    type Output = Poly<T>;

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
        let len_delta = self.raw_len() as isize - rhs.raw_len() as isize;
        (if len_delta == 0 {
            Self(self.0 - rhs.0)
        } else if len_delta < 0 {
            let mut lhs: Array1<T> = Array1::<T>::zeros([len_delta.abs() as usize]);
            lhs.append(Axis(0), self.0.view()).unwrap(); // TODO
            Self(lhs - rhs.0)
        } else {
            let mut rhs_new: Array1<T> = Array1::<T>::zeros([len_delta as usize]);
            rhs_new.append(Axis(0), rhs.0.view()).unwrap(); // TODO
            Self(self.0 - rhs_new)
        })
        .trim_zeros()
    }
}

impl<T: Scalar> Mul for Poly<T> {
    type Output = Poly<T>;

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

impl<T: Scalar> Div for Poly<T> {
    type Output = Poly<T>;

    /// Computes the quotient of two polynomials, truncating the remainder.
    ///
    /// See also `Poly::div_rem()`.
    fn div(self, rhs: Self) -> Self::Output {
        self.div_rem(rhs).0
    }
}

impl<T: Scalar> Rem for Poly<T> {
    type Output = Poly<T>;

    /// Computes the remainder of the division of two polynomials.
    ///
    /// See also `Poly::div_rem()`.
    fn rem(self, rhs: Self) -> Self::Output {
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
            let k = i as isize - j as isize;
            if k >= 0 && k < input_len as isize {
                sum = sum.clone() + input[k as usize].clone() * kernel[j].clone();
            }
        }
        output[i] = sum;
    }

    output
}
