mod aberth_ehrlich;
mod errors;
mod initial_guess;
mod parallel_newton;

use crate::{
    aliases::{C, R},
    base::{BasePoly, UnivariateMarker},
    roots::Roots,
    scalar_traits::{BasicScalar, NonIntegerScalar, RealScalar},
    storage_traits::BasicStorage,
};

use self::{
    aberth_ehrlich::aberth_ehrlich,
    errors::{BreakReason, ContinueReason, RootsError},
    initial_guess::{initial_guess_annulus, initial_guess_uniform},
};

pub enum ComplexRootsInitialGuessStrategy<T> {
    /// Don't do aything, use whichever roots were provided with [`ComplexRootsBuilder::initial_guess_pool`].
    None,

    /// Guesses are uniformly distributed in a square within the given ranges.
    RandomUniform {
        seed: u64,
        min_re: T,
        max_re: T,
        min_im: T,
        max_im: T,
    },

    RandomAnnulus {
        seed: u64,
        bias: T,
        perturbation: T,
    },
    // TODO: Hull {},
    // TODO: GridSearch {},
}

pub enum ComplexRootsFinderStrategyKind {
    /// Use aberth-ehrlich method
    AberthEhrlich,
    ParallelNewton,
}

pub enum ComplexRootsMultiplicityPolicy<T> {
    None,
    BroadcastBest { detection_epsilon: T },
    BroadcastAverage { detection_epsilon: T },
    KeepBest { detection_epsilon: T },
    KeepAverage { detection_epsilon: T },
}

pub struct ComplexRootsFinderStrategy<T> {
    kind: ComplexRootsFinderStrategyKind,
    min_iter: usize,
    max_iter: Option<usize>,
    stop_epsilon: T,
}

pub struct ComplexRootsBuilder<'a, S: BasicStorage>
where
    S::T: NonIntegerScalar,
{
    pub(crate) poly: &'a BasePoly<S, UnivariateMarker>,
    pub(crate) active_guesses: Vec<C<S>>,
    pub(crate) initial_guess_strategy: ComplexRootsInitialGuessStrategy<R<S>>,
    pub(crate) finder_strategies: Vec<ComplexRootsFinderStrategy<R<S>>>,
    pub(crate) multiplicity_policy: ComplexRootsMultiplicityPolicy<R<S>>,
}

impl<S: BasicStorage> ComplexRootsBuilder<'_, S>
where
    S::T: NonIntegerScalar,
{
    /// Provide some fixed guesses, if you know approximately where the roots are.
    pub fn initial_guess_pool(mut self, guesses: &[C<S>]) -> Self {
        self.active_guesses.extend_from_slice(guesses);
        self
    }

    /// How to generate the initial guesses for the root finders to improve.
    pub fn initial_guess_strategy(
        mut self,
        strategy: ComplexRootsInitialGuessStrategy<R<S>>,
    ) -> Self {
        self.initial_guess_strategy = strategy;
        self
    }

    /// Add a root finder strategy. Call this method multiple times to create a
    /// chain of root finder strategies, which are applied in sequence.
    pub fn add_finder_strategy(mut self, strategy: ComplexRootsFinderStrategy<R<S>>) -> Self {
        self.finder_strategies.push(strategy);
        self
    }

    /// What to do with roots with higher multiplicity.
    pub fn multiplicity_policy(mut self, policy: ComplexRootsMultiplicityPolicy<R<S>>) -> Self {
        self.multiplicity_policy = policy;
        self
    }

    pub fn finish(self) -> Result<Roots<C<S>>, RootsError<C<S>>> {
        let Self {
            poly,
            mut active_guesses,
            initial_guess_strategy,
            finder_strategies,
            multiplicity_policy,
        } = self;
        let mut poly = poly.into_owned_poly_inner();
        let mut roots = vec![];

        // find trivial roots
        poly.trim_inner();
        poly.zero_roots_inner(&mut roots);
        poly.linear_roots_inner(&mut roots);
        poly.quadratic_roots_inner(&mut roots);
        if poly.degree_raw_inner() <= 0 {
            // we already found all the roots
            return Ok(Roots(roots.into()));
        }

        // generate initial guesses
        let needed_guesses = poly.degree_raw_inner().saturating_sub(active_guesses.len());
        match initial_guess_strategy {
            ComplexRootsInitialGuessStrategy::None => (),
            ComplexRootsInitialGuessStrategy::RandomUniform {
                seed,
                min_re,
                max_re,
                min_im,
                max_im,
            } => initial_guess_uniform(
                seed,
                min_re,
                max_re,
                min_im,
                max_im,
                &mut active_guesses,
                needed_guesses,
            ),
            ComplexRootsInitialGuessStrategy::RandomAnnulus {
                seed,
                bias,
                perturbation,
            } => initial_guess_annulus(
                seed,
                &poly,
                bias,
                perturbation,
                &mut active_guesses,
                needed_guesses,
            ),
        }
        while active_guesses.len() > poly.degree_raw_inner() {
            active_guesses.pop();
        }
        if active_guesses.len() < poly.degree_raw_inner() {
            return Err(RootsError::InsufficientGuesses);
        }

        // run root finding
        let mut last_status = ContinueReason::Ok;
        for strategy in finder_strategies {
            let flow = match strategy.kind {
                ComplexRootsFinderStrategyKind::AberthEhrlich => aberth_ehrlich(
                    &mut poly,
                    &mut active_guesses,
                    strategy.min_iter,
                    strategy.max_iter,
                    strategy.stop_epsilon,
                ),
                ComplexRootsFinderStrategyKind::ParallelNewton => todo!(),
            };
            match flow {
                std::ops::ControlFlow::Continue(c) => last_status = c,
                std::ops::ControlFlow::Break(b) => return Err(b.into()),
            }
        }
        if last_status != ContinueReason::Ok {
            // last stage did not converge
            return Err(RootsError::from_continue_reason(
                last_status,
                Roots(active_guesses.into()),
            ));
        }

        // handle roots with multiplicity
        match multiplicity_policy {
            ComplexRootsMultiplicityPolicy::None => (),
            ComplexRootsMultiplicityPolicy::BroadcastBest { detection_epsilon } => todo!(),
            ComplexRootsMultiplicityPolicy::BroadcastAverage { detection_epsilon } => todo!(),
            ComplexRootsMultiplicityPolicy::KeepBest { detection_epsilon } => todo!(),
            ComplexRootsMultiplicityPolicy::KeepAverage { detection_epsilon } => todo!(),
        }

        Ok(Roots(active_guesses.into()))
    }
}
