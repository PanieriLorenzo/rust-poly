/// An iterator that counts up until it reaches max, at which point it saturates
///
/// This is an endless iterator.

#[inline]
pub fn saturating_counter() -> impl Iterator<Item = usize> {
    (0..usize::MAX).chain(std::iter::repeat(usize::MAX))
}
