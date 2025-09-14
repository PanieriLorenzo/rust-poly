use num_traits::Zero;

use crate::{
    scalar_traits::BasicScalar,
    storage_traits::{BasicStorage, MutStorage, OwnedStorage},
};

impl<T: BasicScalar, const N: usize> BasicStorage for [T; N] {
    type T = T;
}

impl<T: BasicScalar> BasicStorage for &[T] {
    type T = T;
}

impl<T: BasicScalar> BasicStorage for Vec<T> {
    type T = T;
}

impl<T: BasicScalar> OwnedStorage for Vec<T> {
    fn truncate(&mut self, new_length: usize) {
        self.truncate(new_length);
    }

    fn truncate_start(&mut self, new_length: usize) {
        // TODO: consider `remove(0)` strategy for small vectors / removing only
        // a few elements.
        let len_delta = self.len().saturating_sub(new_length);
        let range = len_delta..self.len();
        *self = self.drain(range).collect();
    }
}

impl<T: BasicScalar> MutStorage for Vec<T> {}
