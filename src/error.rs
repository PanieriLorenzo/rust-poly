use std::fmt;

use thiserror::Error;

#[derive(Debug, Error)]
pub(crate) enum ErrorKind {
    /// Use this for when the user-provided maxiter is reached
    #[error("did not converge within the given number of iterations")]
    MaxIterUser,

    /// Use this for when an internally defined maxiter is reached
    #[error("did not converge")]
    SlowConvergence,

    /// Use this when the outermost defined maxiter is reached, but this wasn't
    /// specified by the user
    #[error("did not converge")]
    MaxIterOuter,

    /// Use this for errors that only appear on pathological inputs
    #[error("did not converge")]
    Pathological,

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

/// The top-level error type for this crate.
#[derive(Debug, Error)]
pub struct Error {
    #[source]
    pub(crate) source: ErrorKind,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.source.fmt(f)
    }
}

impl Error {
    pub(crate) fn max_iter_user() -> Self {
        Self {
            source: ErrorKind::MaxIterUser,
        }
    }

    pub(crate) fn max_iter_inner() -> Self {
        Self {
            source: ErrorKind::SlowConvergence,
        }
    }

    pub(crate) fn max_iter_outer() -> Self {
        Self {
            source: ErrorKind::MaxIterOuter,
        }
    }

    pub(crate) fn pathological() -> Self {
        Self {
            source: ErrorKind::Pathological,
        }
    }

    /// Maps [`ErrorKind::MaxIterOuter`] to [`ErrorKind::MaxIterInner`]
    pub(crate) fn map_inner(self) -> Self {
        match self.source {
            ErrorKind::MaxIterOuter => Self::max_iter_inner(),
            _ => self,
        }
    }
}
