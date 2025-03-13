use crate::storage::{BaseStore, OwnedStore};

use super::roots_base::RootsBase;

pub struct RootsOutput<T, S>
where
    S: BaseStore<T>,
    <S as ToOwned>::Owned: OwnedStore<T>,
{
    roots: RootsBase<T, S>,
    stop_reason: StopReason,
}

#[non_exhaustive]
pub enum StopReason {
    NoConverge,
    Other(&'static str),
}
