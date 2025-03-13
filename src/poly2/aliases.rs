use super::{poly_base::PolyBase, roots_base::RootsBase};

/// A growable univariate polynomial.
pub type DynUniPoly<T> = PolyBase<T, Vec<T>>;

/// A growable set of univariate polynomial roots.
pub type DynUniRoots<T> = RootsBase<T, Vec<T>>;
