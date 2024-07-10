use super::{Complex, Poly, RealScalar};

impl<T: RealScalar> Poly<T> {
    /// Index from the end, useful for porting algorithm that use the descending convention
    pub(crate) fn coeffs_descending(&self, idx: usize) -> &Complex<T> {
        &self.0[self.len_raw() - idx - 1]
    }

    pub(crate) fn coeffs_descending_mut(&mut self, idx: usize) -> &mut Complex<T> {
        let n = self.len_raw();
        &mut self.0[n - idx - 1]
    }
}

impl<T: RealScalar> Poly<T> {
    /// Return a slice containing the coefficients in ascending order of degree
    ///
    /// This is an alias for [`Poly::as_slice`] for API consistency.
    #[inline]
    #[must_use]
    pub fn coeffs(&self) -> &[Complex<T>] {
        self.as_slice()
    }

    /// Return a mutable slice containing the coefficient in ascending order of degree
    ///
    /// This is an alias for [`Poly::as_mut_slice()`] for API consistency.
    #[inline]
    pub fn coeffs_mut(&mut self) -> &mut [Complex<T>] {
        self.as_mut_slice()
    }
}
