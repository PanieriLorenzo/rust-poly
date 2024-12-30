use super::poly_base::PolyBase;

/// A growable univariate polynomial.
pub type DynUniPoly<T> = PolyBase<Vec<T>>;