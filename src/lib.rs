use ndarray::{Array, Array1, Dimension, ScalarOperand};
use num_traits::Num;

/// polynomial as a list of coefficients of terms of descending degree
#[derive(Clone, PartialEq, Debug)]
pub struct Poly<T: Num + Clone + ScalarOperand>(Array1<T>);

impl<T: Num + Clone + ScalarOperand> Poly<T> {
    pub fn new(coeffs: Array1<T>) -> Self {
        Self(coeffs)
    }

    /// evaluate a polynomial at a specific input value `x`. This may be an
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
