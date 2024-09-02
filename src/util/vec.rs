use std::ops::{Add, Div};

use num::{FromPrimitive, Zero};

pub fn slice_mean<T: Zero + Add<Output = T> + Div<Output = T> + FromPrimitive + Clone>(
    v: &[T],
) -> T {
    let num = v.len();
    v.iter().fold(T::zero(), |a, b| a + b.clone()) / T::from_usize(num).expect("overflow")
}
