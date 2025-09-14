mod inner;
mod traits;

use crate::{
    base::{BasePoly, UnivariateMarker},
    root_finders::{
        ComplexRootsBuilder, ComplexRootsInitialGuessStrategy, ComplexRootsMultiplicityPolicy,
    },
    scalar_traits::NonIntegerScalar,
    storage_traits::{BasicStorage, OwnedStorage},
};

impl<S: BasicStorage> BasePoly<S, UnivariateMarker> {
    /// Find all complex roots of this polynomial.
    pub fn complex_roots(&self) -> ComplexRootsBuilder<'_, S>
    where
        S::T: NonIntegerScalar,
    {
        ComplexRootsBuilder {
            poly: self,
            active_guesses: vec![],
            initial_guess_strategy: ComplexRootsInitialGuessStrategy::None,
            finder_strategies: vec![],
            multiplicity_policy: ComplexRootsMultiplicityPolicy::None,
        }
    }
}

impl<S: OwnedStorage> BasePoly<S, UnivariateMarker> {}
