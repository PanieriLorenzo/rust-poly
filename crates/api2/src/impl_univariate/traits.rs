//! Trait implementations for univariate polynomials

use std::borrow::Borrow;

use crate::{
    base::{BasePoly, UnivariateMarker},
    storage_traits::{BasicStorage, OwnedStorage},
};

impl<S: OwnedStorage> FromIterator<S::T> for BasePoly<S, UnivariateMarker> {
    fn from_iter<T: IntoIterator<Item = S::T>>(iter: T) -> Self {
        // SAFETY: univariate polynomials trivially map to a contiguous storage
        // with no additional checks required.
        let mut ret = unsafe { Self::from_parts(S::from_iter(iter), 1) };
        ret.normalize_inner();
        ret
    }
}

impl<S: BasicStorage, M> PartialEq for BasePoly<S, M> {
    fn eq(&self, other: &Self) -> bool {
        self.storage.borrow() == other.storage.borrow() && self.vars == other.vars
    }
}
