//! Macros for reducing doc comment boilerplate.

/// Documents `from_f64` panics.
macro_rules! panic_t_from_f64 {
    () => {
        r"If `T` is any primitive type, this function does not panic. However, if `T` does not implement [`num::FromPrimitive::from_f64`] correctly, this function might panic under certain circumstances (casting of extreme values, usually)\n\n"
    };
}
pub(crate) use panic_t_from_f64;

/// Documents generic [`num::FromPrimitive`] panics.
///
/// You must provide a string containing the integer type with the largest range
/// of values that should be supported.
macro_rules! panic_t_from_int {
    ($ty_str:expr) => {
        concat!(r"If `T` is `f32` or `f64`, this function does not panic. However, if `T` cannot represent all primitive integers smaller than [`", $ty_str, r"::MAX`], this might panic under certain circumstances (casting of extreme values, usually)\n\n")
    };
}
pub(crate) use panic_t_from_int;

/// Documents panics due to absurdly large polynomials
macro_rules! panic_absurd_size {
    () => {
        r"May theoretically panic for absurdly large polynomials, however such polynomials will likely not fit in memory anyway.\n\n"
    }
}
pub(crate) use panic_absurd_size;

/// Default explanation for [`roots::Error::NoConverge`] errors.
macro_rules! errors_no_converge {
    () => {
        r"- `NoConverge`: the root finder did not converge within the given constraints. The best guess so far is returned with the error.\n"
    };
}
pub(crate) use errors_no_converge;
