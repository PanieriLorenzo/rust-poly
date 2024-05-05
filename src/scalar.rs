use num::Num;
use std::ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign};

/// The trait bounds necessary to provide the basic functionality of this crate.
pub trait Scalar: Clone + PartialEq + std::fmt::Debug + Num + 'static {}
impl<T: Clone + PartialEq + std::fmt::Debug + Num + 'static> Scalar for T {}

// TODO: these are required by nalgebra for things that shouldn't require them.
//       perhaps in the future they can be dropped?
/// Trait bounds necessary to provide more advanced mathematical features.
#[allow(clippy::module_name_repetitions)]
pub trait ScalarOps: Scalar + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign {}
impl<T: Scalar + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign> ScalarOps for T {}
