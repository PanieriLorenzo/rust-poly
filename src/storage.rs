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
    /// Non-owning iterator over all elements in the container
    ///
    /// For multi-dimensional
    /// containers, the order of the coefficients is undefined for now, but will
    /// be stabilized once/if multi-variate polynomials are implemented.
    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a T>
    where
        T: 'a;

    /// Slice into the backing store's data, the data should be stored as one contiguous region.
    ///
    /// For multi-dimensional
    /// containers, the order of the coefficients is undefined for now, but will
    /// be stabilized once/if multi-variate polynomials are implemented.
    fn as_slice(&self) -> &[T];

    fn shape(&self) -> Box<[usize]>;

    fn ndim(&self) -> usize;
}

/// Mutable storage, either borrowed or owned.
pub trait MutStore<T>: BaseStore<T>
where
    Self::Owned: OwnedStore<T>,
{
    /// Slice into the backing store's data, the data should be stored as one contiguous region
    fn as_mut_slice(&mut self) -> &mut [T];
}

/// Univariate storage.
pub trait UniStore<T>: BaseStore<T>
where
    Self::Owned: OwnedUniStore<T>,
{
}

/// Owned storage.
pub trait OwnedStore<T>: MutStore<T, Owned = Self> + Sized {
    fn empty() -> Self;

    fn zeros(shape: &[usize]) -> Self
    where
        T: Zero;

    /// Construct a store from a sequence of values and the shape, i.e. dimensions,
    /// that the values should be arranged in. Returns `None` if the values do
    /// not fit in the shape.
    fn from_iter(shape: &[usize], values: impl IntoIterator<Item = T>) -> Option<Self>;
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
