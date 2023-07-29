use crate::Scalar;

#[inline(always)]
pub fn neg<T: Scalar>(x: T) -> T {
    T::zero() - x
}
