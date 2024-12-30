//! Backing store traits

pub mod impl_storage;

/// The base storage trait, makes no assumptions about ownership or lifetimes.
///
/// This is the trait of borrowed and immutable storage, which is also the base
/// for mutable and owned representations. It is also the backing store for
/// multi-variate polynomials, if these are ever implemented.
///
/// It is based on the design of the `ndarray` crate, but only contains what is
/// needed for polynomials.
pub trait BaseStore {
    type Elem;
}

/// Univariate storage.
pub trait UniStore {}

/// Growable storage.
pub trait DynStore {}
