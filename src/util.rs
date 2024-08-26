//! Internal utilities, not part of the API

pub(crate) mod big_float;
pub(crate) mod casting;
pub(crate) mod complex;
pub(crate) mod doc_macros;
pub(crate) mod float;
pub(crate) mod iterator;
pub(crate) mod linalg;
pub(crate) mod luts;

// re-exported by crate root
#[doc(hidden)]
pub mod __testing;
