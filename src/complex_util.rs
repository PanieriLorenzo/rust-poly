// internal utilities for dealing with Complex annoyiances

use num_traits::Zero;

use crate::{Scalar, C};

// neg operator for Complex, as it does not implement std::ops::Neg
#[inline(always)]
pub(crate) fn c_neg<T: Scalar>(x: C<T>) -> C<T> {
    C::<T>::zero() - x
}
