use thiserror::Error;

use crate::{roots::Roots, storage_traits::BasicStorage};

mod sealed {
    pub trait Sealed {}
}

#[derive(Debug, Error)]
pub enum RootsError<T> {
    #[error("reached max_iter without converging")]
    MaxIter(Roots<T>),

    #[error("could not improve roots any further")]
    NoProgress(Roots<T>),

    #[error("not enough initial guesses were provided")]
    InsufficientGuesses,

    #[error("cannot run root finder with duplicate guesses")]
    DuplicateGuesses,
}

impl<T> RootsError<T> {
    pub(crate) fn from_continue_reason(reason: ContinueReason, roots: Roots<T>) -> Self {
        match reason {
            ContinueReason::Ok => unreachable!("reason is not an error"),
            ContinueReason::MaxIter => Self::MaxIter(roots),
        }
    }
}

impl<T> From<BreakReason> for RootsError<T> {
    fn from(value: BreakReason) -> Self {
        match value {
            BreakReason::DuplicateGuesses => Self::DuplicateGuesses,
        }
    }
}

/// Extension trait for `Result<Roots<T>, RootsError<Roots<T>>>`
pub trait RootsResultExt<T>: sealed::Sealed {
    /// Get roots even if an error was raised. Returns [`None`] if an error
    /// occurred without producing any roots.
    fn ignore_errors(self) -> Option<Roots<T>>;
}

impl<T> sealed::Sealed for Result<Roots<T>, RootsError<T>> {}

impl<T> RootsResultExt<T> for Result<Roots<T>, RootsError<T>> {
    fn ignore_errors(self) -> Option<Roots<T>> {
        match self {
            Ok(roots) | Err(RootsError::MaxIter(roots) | RootsError::NoProgress(roots)) => {
                Some(roots)
            }
            _ => None,
        }
    }
}

#[derive(PartialEq, Eq)]
pub(crate) enum ContinueReason {
    Ok,
    MaxIter,
}

#[derive(PartialEq, Eq)]
pub(crate) enum BreakReason {
    DuplicateGuesses,
}
