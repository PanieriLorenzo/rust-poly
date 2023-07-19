use crate::{ComplexScalar, Poly, RealScalar};
use ndarray::Array1;
use num_complex::Complex;

pub trait Roots {
    type Output;
    fn roots(&self) -> Array1<Self::Output>;
}

impl<R: RealScalar> Roots for Poly<R> {
    type Output = Complex<R>;

    fn roots(&self) -> Array1<Self::Output> {
        todo!()
    }
}
