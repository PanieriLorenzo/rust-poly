use std::marker::PhantomData;

use crate::storage::{BaseStore, OwnedStore};

// HACK: T should have been an associated type on BaseStore, but Rust has an ICE
//       when dealing with recursive associated type trait bounds, so it gets very
//       dicey to specify the bounds correctly, hence why there's a _phantom.
pub struct RootsBase<T, S>
where
    S: BaseStore<T>,
    <S as ToOwned>::Owned: OwnedStore<T>,
{
    _phantom: PhantomData<T>,
    storage: S,
}
