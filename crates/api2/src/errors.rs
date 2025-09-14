// TODO: move these into the correct module, this top-level is for top-level public errors.

pub const DEGREE_TOO_LARGE: &str = "degrees larger than i32::MAX are not supported";
pub const CAST_OVERFLOW: &str = "overflow while casting";
pub const INFALLIBLE_CONVERSION: &str = "conversion should never fail";
