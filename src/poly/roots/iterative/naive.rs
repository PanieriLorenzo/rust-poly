/// Naive Newton's method.
///
/// It's strongly recommended to use a different root finder, as naive Newton
/// iteration is not robust at all. Use this only if dealing with small well-behaved
/// polynomials were the simplicity gives a performance advantage, or as a baseline
/// for benchmarking purposes.
///
/// Unlike [`roots::Madsen`] this has none of the tweaks that improve convergence
/// on pathological roots. This will often get stuck on any slightly challenging
/// root.
pub struct Naive {}
