use num::Complex;

use crate::{Poly, Scalar};

impl<T> Poly<T> {
    #[must_use]
    pub const fn from_dvector(value: na::DVector<Complex<T>>) -> Self {
        Self(value)
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
    pub fn to_dvector(self) -> na::DVector<Complex<T>> {
        self.0
    }

    /// Iterate over coefficients, from the least significant
    pub fn iter(&self) -> std::slice::Iter<'_, na::Complex<T>> {
        self.0.as_slice().iter()
    }

    /// Iterate over coefficients, from the least significant
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, na::Complex<T>> {
        self.0.as_mut_slice().iter_mut()
    }
}

impl<T: Clone> Poly<T> {
    #[must_use]
    pub fn to_vec(&self) -> Vec<Complex<T>> {
        Vec::from(self.as_slice())
    }
}

impl<T: Scalar> Poly<T> {
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
        // TODO: avoid unnecessary initialization of Vec
        let temp_vec: Vec<_> = value.iter().map(Complex::from).collect();
        Self::new(&temp_vec)
    }

    #[allow(clippy::needless_pass_by_value)]
    #[must_use]
    pub fn from_real_vec(value: Vec<T>) -> Self {
        Self::from_real_slice(value.as_slice())
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

impl<'a, T> IntoIterator for &'a Poly<T> {
    type IntoIter = std::slice::Iter<'a, na::Complex<T>>;
    type Item = &'a na::Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut Poly<T> {
    type IntoIter = std::slice::IterMut<'a, na::Complex<T>>;
    type Item = &'a mut na::Complex<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}
