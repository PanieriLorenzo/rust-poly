use num::cast;

use crate::Scalar;

pub(crate) fn usize_to_i32(x: usize) -> i32 {
    x.try_into().unwrap_or(i32::MAX)
}

pub(crate) fn usize_to_u32(x: usize) -> u32 {
    x.try_into().unwrap_or(u32::MAX)
}

pub(crate) fn usize_to_f64(x: usize) -> f64 {
    cast(x).unwrap_or(f64::INFINITY)
}

/// Convert usize to Scalar without any primitive casting
#[deprecated(note = "use num::FromPrimitive::from_usize instead")]
#[must_use]
pub fn usize_to_scalar<T: Scalar>(x: usize) -> T {
    T::from_usize(x).unwrap()
}
