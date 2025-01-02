use num::Zero;

use super::{BaseStore, MutStore, OwnedStore, OwnedUniStore, UniStore};

impl<T: Clone> BaseStore<T> for Vec<T> {
    #[inline]
    #[must_use]
    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a T>
    where
        T: 'a,
    {
        self.as_slice().iter()
    }
}

impl<T: Clone> UniStore<T> for Vec<T> {}

impl<T: Clone> OwnedStore<T> for Vec<T> {
    fn empty() -> Self {
        vec![]
    }

    fn zeros(shape: &[usize]) -> Self
    where
        T: Zero,
    {
        debug_assert_eq!(shape.len(), 1, "Vec is 1 dimensional");
        vec![T::zero(); shape[0]]
    }
}

impl<T: Clone> OwnedUniStore<T> for Vec<T> {
    #[inline]
    fn push(&mut self, val: T) {
        self.push(val);
    }
}

impl<T: Clone> MutStore<T> for Vec<T> {}
