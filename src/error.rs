use std::fmt;

use thiserror::Error;

#[derive(Debug, Error)]
pub(crate) enum ErrorKind {
    /// Use this for when the user-provided maxiter is reached
    #[error("did not converge within the given number of iterations")]
    MaxIterUser,

    /// Use this for when an internally defined maxiter is reached
    #[error("did not converge")]
    MaxIterInner,

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

/// The top-level error type for this crate.
#[derive(Debug, Error)]
pub struct Error {
    do_recover: bool,

    #[source]
    source: ErrorKind,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.source.fmt(f)
    }
}

impl Error {
    pub(crate) fn max_iter_user(do_recover: bool) -> Self {
        Self {
            do_recover,
            source: ErrorKind::MaxIterUser,
        }
    }

    pub(crate) fn max_iter_inner(do_recover: bool) -> Self {
        Self {
            do_recover,
            source: ErrorKind::MaxIterInner,
        }
    }
}
