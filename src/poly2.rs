//! The new API, will replace the [`poly`] module.

use crate::{num::Zero, util::doc_macros::panic_absurd_size};

pub mod aliases;
pub mod poly_base;

/// A generic polynomial
///
/// This trait makes very few assumptions. May be univariate or multivariate,
/// integer, real or complex-valued, owned or borrowed, etc...
///
/// Must have a valid way to represent a zero polynomial, this is expressed by
/// implementing [`crate::num::Zero`]. A zero polynomial behaves like the neutral
/// element of addition over polynomials. A zero polynomial has degree -1 by
/// convention.
pub trait Poly<T>: Zero {
    /// The owned version of this polynomial. May be `Self` for representations
    /// that are already owned.
    type OwnedRepr;

    /// Return the degree as `usize`.
    ///
    /// Note that unlike [`Poly::degree`], this will saturate at 0 for zero
    /// polynomials. As the degree of zero polynomials is undefined.
    fn degree_usize(&self) -> usize;

    /// The degree of a polynomial (the maximum exponent)
    ///
    /// Note that this will return `-1` for zero polynomials. The degree of
    /// zero polynomials is undefined, but we use the `-1` convention adopted
    /// by some authors \[which authors?\].
    ///
    /// # Panics
    #[doc = panic_absurd_size!()]
    fn degree(&self) -> i64 {
        if self.is_zero() {
            return -1;
        }
        self.degree_usize()
            .try_into()
            .expect("usize did not fit into i64")
    }
}
