use itertools::Itertools;
use num::Complex;

use crate::{Poly, RealScalar};

impl<T: RealScalar> Poly<T> {
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

    /// Iterate over coefficients, from the least significant
    pub fn iter(&self) -> std::slice::Iter<'_, Complex<T>> {
        self.0.as_slice().iter()
    }

    /// Iterate over coefficients, from the least significant
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Complex<T>> {
        self.0.as_mut_slice().iter_mut()
    }

    #[must_use]
    pub fn to_vec(&self) -> Vec<Complex<T>> {
        Vec::from(self.as_slice())
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

    #[must_use]
    pub fn from_real_slice(value: &[T]) -> Self {
        Self::from_real_iterator(value.iter().cloned())
    }

    #[allow(clippy::needless_pass_by_value)]
    #[must_use]
    pub fn from_real_vec(value: Vec<T>) -> Self {
        Self::from_real_slice(value.as_slice())
    }

    #[must_use]
    pub fn from_real_iterator(coeffs: impl Iterator<Item = T>) -> Self {
        let res = Self::from_complex_iterator(coeffs.map(|x| Complex::from(x)));
        debug_assert!(res.is_normalized());
        res
    }

    #[must_use]
    pub fn from_complex_iterator(coeffs: impl Iterator<Item = Complex<T>>) -> Self {
        Self(coeffs.collect_vec()).normalize()
    }
}

impl<T: RealScalar> From<&[Complex<T>]> for Poly<T> {
    fn from(value: &[Complex<T>]) -> Self {
        Self::from_complex_slice(value)
    }
}

impl<T: RealScalar> From<Vec<Complex<T>>> for Poly<T> {
    fn from(value: Vec<Complex<T>>) -> Self {
        Self::from_complex_vec(value)
    }
}

// TODO: these are ambiguous and lead to making a `Poly<Complex<Complex<T>>>`
//       but making a `Real` trait would preventing making blanket impls for
//       numeric types.
// impl<T: Scalar> From<&[T]> for Poly<T> {
//     fn from(value: &[T]) -> Self {
//         Self::from_real_slice(value)
//     }
// }

// impl<T: Scalar> From<Vec<T>> for Poly<T> {
//     fn from(value: Vec<T>) -> Self {
//         Self::from_real_vec(value)
//     }
// }

impl<T: RealScalar> From<Poly<T>> for *const Complex<T> {
    fn from(val: Poly<T>) -> Self {
        val.as_ptr()
    }
}

impl<T: RealScalar> From<Poly<T>> for *mut Complex<T> {
    fn from(mut val: Poly<T>) -> Self {
        val.as_mut_ptr()
    }
}

impl<T: RealScalar> From<Poly<T>> for Vec<Complex<T>> {
    fn from(val: Poly<T>) -> Self {
        val.to_vec()
    }
}

impl<'a, T: RealScalar> IntoIterator for &'a Poly<T> {
    type IntoIter = std::slice::Iter<'a, Complex<T>>;
    type Item = &'a Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T: RealScalar> IntoIterator for &'a mut Poly<T> {
    type IntoIter = std::slice::IterMut<'a, Complex<T>>;
    type Item = &'a mut Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}
