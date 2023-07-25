use crate::Scalar;

#[inline(always)]
pub(crate) fn neg<T: Scalar>(x: T) -> T {
    T::zero() - x
}
