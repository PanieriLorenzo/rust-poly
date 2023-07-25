// internal utilities for dealing with Complex annoyiances

use num_complex::Complex;
use num_traits::Zero;

use crate::Scalar;

// neg operator for Complex, as it does not implement std::ops::Neg
#[inline(always)]
pub(crate) fn c_neg<T: Scalar>(x: Complex<T>) -> Complex<T> {
    Complex::<T>::zero() - x
}
