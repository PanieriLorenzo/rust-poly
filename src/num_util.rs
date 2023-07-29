use crate::Scalar;

pub fn neg<T: Scalar>(x: T) -> T {
    T::zero() - x
}
