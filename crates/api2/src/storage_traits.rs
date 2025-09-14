use std::{
    borrow::{Borrow, BorrowMut},
    fmt::Debug,
    ops::{Add, Index, Range},
};

use num_traits::Zero;

use crate::scalar_traits::BasicScalar;

/// A contiguous container, i.e. not a sparse container.
pub trait BasicStorage: Borrow<[Self::T]> + Debug {
    type T: BasicScalar;
}

/// A [`BasicStorage`] that can also be borrowed mutably, though not dynamically sized.
pub trait MutStorage: BasicStorage + BorrowMut<[Self::T]> {}

/// A [`BasicStorage`] that owns its data
pub trait OwnedStorage: MutStorage + FromIterator<Self::T> {
    /// Truncate to a new length by removing items from the end.
    fn truncate(&mut self, new_length: usize);

    /// Truncate to a new length by removing items from the beginning.
    fn truncate_start(&mut self, new_length: usize);
}
