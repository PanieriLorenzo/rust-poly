use std::marker::PhantomData;

use crate::storage_traits::BasicStorage;

/// Mark polynomial as univariate, which is more efficient than a multivariate
/// polynomial with 1 variable.
#[derive(Debug)]
pub enum UnivariateMarker {}

#[derive(Debug)]
pub struct BasePoly<S, M> {
    pub(crate) storage: S,
    pub(crate) vars: usize,
    _phantom_marker: PhantomData<M>,
}

impl<S, M> BasePoly<S, M> {
    /// # Safety
    /// This method is actually safe, but we reserve the right to make it unsafe
    /// for performance reasons at any moment, if we decide that unchecked indexing
    /// of the internal storage is necessary for performance reasons.
    ///
    /// Creating an univariate polynomial this way is always safe, but a multivariate
    /// polynomial should have a storage with size equal to `(degree + 1) * vars`.
    pub unsafe fn from_parts(storage: S, vars: usize) -> Self {
        Self {
            storage,
            vars,
            _phantom_marker: PhantomData,
        }
    }
}
