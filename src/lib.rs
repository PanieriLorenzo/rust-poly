#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]

mod scalar;
pub use scalar::Scalar;

// mod roots;
// pub use roots::Roots;

mod complex_util;
use complex_util::c_neg;
mod array_util;
use array_util::np_diag;
mod impl_num;

extern crate nalgebra as na;
