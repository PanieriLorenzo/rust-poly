use ndarray::ScalarOperand;
use num_complex::{Complex32, Complex64};
use num_traits::Num;

pub trait Scalar: Num + Clone + ScalarOperand + core::fmt::Debug {}

macro_rules! impl_scalar {
    ($($t:ty),* $(,)?) => {
        $(
            impl Scalar for $t {}
        )*
    };
}

impl_scalar!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f32, f64, Complex32, Complex64
);
