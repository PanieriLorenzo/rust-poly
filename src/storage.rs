//! Backing store traits

use num::Zero;

pub mod impl_storage;

/// The base storage trait, makes no assumptions about ownership or lifetimes.
///
/// This is the trait of borrowed and immutable storage, which is also the base
/// for mutable and owned representations. It is also the backing store for
/// multi-variate polynomials, if these are ever implemented.
///
/// It is based on the design of the `ndarray` crate, but only contains what is
/// needed for polynomials.
pub trait BaseStore<T>: ToOwned
where
    Self::Owned: OwnedStore<T>,
{
    /// Non-owning iterator
    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a T>
    where
        T: 'a;

    fn shape(&self) -> &[usize];

    fn ndim(&self) -> usize;
}

/// Mutable storage, either borrowed or owned.
pub trait MutStore<T>: BaseStore<T>
where
    Self::Owned: OwnedStore<T>,
{
}

/// Univariate storage.
pub trait UniStore<T>: BaseStore<T>
where
    Self::Owned: OwnedUniStore<T>,
{
}

/// Owned storage.
pub trait OwnedStore<T>: MutStore<T, Owned = Self> {
    fn empty() -> Self;

    fn zeros(shape: &[usize]) -> Self
    where
        T: Zero;

    fn from_iter(shape: &[usize], values: impl IntoIterator<Item = T>) -> Self;
}

/// Growable uni-dimensional storage, e.g. [`Vec`].
pub trait OwnedUniStore<T>: UniStore<T> + OwnedStore<T> {
    fn push(&mut self, val: T);

    fn extend(&mut self, values: impl IntoIterator<Item = T>) {
        for val in values {
            self.push(val);
        }
    }

    // fn from_iter(values: impl IntoIterator<Item = T>) -> Self
    // where
    //     Self: Sized,
    // {
    //     let mut v = Self::empty();
    //     v.extend(values);
    //     v
    // }
}
