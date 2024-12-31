//! The new API, will replace the [`poly`] module.

use std::ops::{Add, Div, Mul, Rem, Sub};

use num::CheckedDiv;

use crate::{
    num::Zero,
    poly2,
    storage::{BaseStore, OwnedStore, OwnedUniStore, UniStore},
    util::doc_macros::panic_absurd_size,
};

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
pub trait Poly<T>:
    ToOwned
    + Add<Self, Output = Self::Owned>
    + Mul<Self, Output = Self::Owned>
    + Sub<Self, Output = Self::Owned>
    + Div<Self, Output = Self::Owned>
    + Rem<Self, Output = Self::Owned>
where
    Self::Owned: OwnedPoly<T>,
    Self: Sized,
    <<Self as poly2::Poly<T>>::BackingStorage as ToOwned>::Owned: OwnedStore<T>,
{
    type BackingStorage: BaseStore<T>;

    /// Implementation detail
    #[doc(hidden)]
    fn _as_store(&self) -> &Self::BackingStorage;

    /// Implementation detail
    #[doc(hidden)]
    fn _from_store(store: Self::BackingStorage) -> Self;

    fn zero() -> Self::Owned;

    fn one() -> Self::Owned;

    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

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

    /// Raises a polynomial to an integer power.
    ///
    /// # Caveats
    /// We adopt the convention that $0^0=1$, even though some authors leave this
    /// case undefined. We believe this to be more useful as it naturally arises
    /// in integer exponentiation when defined as repeated multiplication.
    #[must_use]
    fn pow(&self, pow: u32) -> Self::Owned {
        self.pow_usize(pow as usize)
    }

    /// Same as [`Poly::pow`], but takes a `usize` exponent.
    fn pow_usize(&self, pow: usize) -> Self::Owned;

    fn terms(&self) -> impl Iterator<Item = Self::Owned>;

    /// Return an iterator over the coefficients in ascending order of degree
    #[inline]
    #[must_use]
    fn coeffs<'a>(&'a self) -> impl Iterator<Item = &'a T>
    where
        T: 'a,
    {
        self._as_store().iter()
    }

    /// Polynomial composition.
    fn compose(&self, other: &Self) -> Self::Owned;

    fn eval(&self, x: T) -> T;

    /// Evaluate the polynomial for each entry of a slice. May be faster than
    /// repeatedly calling [`Poly::eval`].
    fn eval_multiple(&self, points: &[T], out: &mut [T])
    where
        T: Clone,
    {
        // TODO: parallelize this loop
        for (y, x) in out.iter_mut().zip(points) {
            *y = self.eval(x.clone());
        }
    }
}

/// Univariate polynomial
pub trait UniPoly<T>: Poly<T>
where
    Self::Owned: OwnedPoly<T>,
    Self::BackingStorage: UniStore<T>,
    // HACK: can't this diamond be flattened?
    <Self::BackingStorage as ToOwned>::Owned: OwnedUniStore<T>,
    <Self::Owned as Poly<T>>::BackingStorage: OwnedUniStore<T>,
{
    fn shift_up(&self, n: usize) -> Self::Owned
    where
        T: Zero + Clone,
    {
        let mut v = <Self::Owned as Poly<T>>::BackingStorage::zeros(&[n]);
        v.extend(self.coeffs().cloned());
        Self::Owned::_from_store(v)
    }
}

pub trait OwnedPoly<T>: Poly<T, Owned = Self> {}
