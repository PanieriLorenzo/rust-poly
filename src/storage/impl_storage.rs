use num::Zero;

use super::{BaseStore, DynStore, DynUniStore, OwnedStore, UniStore};

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

impl<T> UniStore for Vec<T> {}

impl<T: Clone> OwnedStore<T> for Vec<T> {
    fn zeros(shape: &[usize]) -> Self
    where
        T: Zero,
    {
        debug_assert_eq!(shape.len(), 1, "Vec is 1 dimensional");
        vec![T::zero(); shape[0]]
    }
}

impl<T: Clone> DynStore<T> for Vec<T> {}

impl<T: Clone> DynUniStore<T> for Vec<T> {
    #[inline]
    fn push(&mut self, val: T) {
        self.push(val);
    }
}
