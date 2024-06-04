//! Internal utilities, not part of the API

// TODO: I don't like how much stuff we are forced to make "public" but hidden
//       we should figure out how to make these private

pub mod casting;
pub mod complex;
pub mod linalg;
pub mod luts;

pub mod testing;

pub(crate) mod float;
