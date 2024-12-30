//! Backing representation for all crate-provided impls of [`Poly2`] trait.

use crate::storage::BaseStore;

pub struct PolyBase<S>
where
    S: BaseStore,
{
    storage: S,
}
