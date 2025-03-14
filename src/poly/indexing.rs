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
