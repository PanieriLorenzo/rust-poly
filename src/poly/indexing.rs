use std::ops::{
    Bound, Index, Range, RangeBounds, RangeFrom, RangeFull, RangeInclusive, RangeTo,
    RangeToInclusive,
};

use super::{Complex, Poly, RealScalar};

mod sealed {
    pub trait Sealed {}
}

pub trait Get<I, T: RealScalar>: sealed::Sealed {
    fn get(&self, idx: I) -> Option<Poly<T>>;
}

impl<T: RealScalar> sealed::Sealed for Poly<T> {}

impl<T: RealScalar> Poly<T> {
    /// Index from the end, useful for porting algorithm that use the descending convention
    pub(crate) fn coeff_descending(&self, idx: usize) -> &Complex<T> {
        &self.0[self.len_raw() - idx - 1]
    }

    pub(crate) fn coeff_descending_mut(&mut self, idx: usize) -> &mut Complex<T> {
        let n = self.len_raw();
        &mut self.0[n - idx - 1]
    }

    /// Implementation for all range-based indexing (because Rust is super annoying
    /// around the different range iterators, see [#3550](https://github.com/rust-lang/rfcs/pull/3550))
    fn get_range_inner(&self, idx_start: Bound<&usize>, idx_end: Bound<&usize>) -> Option<Self> {
        debug_assert!(self.is_normalized());

        let start = match idx_start {
            Bound::Included(x) => *x,
            Bound::Excluded(_) => panic!("range start can't be exclusive"),
            Bound::Unbounded => 0,
        };
        let end = match idx_end {
            Bound::Included(x) => *x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.len_raw(),
        };
        if start >= self.len_raw() || end > self.len_raw() {
            return None;
        }
        let coeffs = &self.as_slice()[start..end];
        Some(Self::from_complex_slice(coeffs).shift_up(start))
    }
}

impl<T: RealScalar> Get<usize, T> for Poly<T> {
    fn get(&self, idx: usize) -> Option<Self> {
        debug_assert!(self.is_normalized());
        self.terms().nth(idx)
    }
}

impl<T: RealScalar> Get<isize, T> for Poly<T> {
    #[allow(clippy::cast_possible_wrap)]
    #[allow(clippy::cast_sign_loss)]
    fn get(&self, idx: isize) -> Option<Self> {
        debug_assert!(self.is_normalized());

        if idx >= 0 {
            return self.get(idx as usize);
        }

        // find the index from the end
        let idx = self.len_raw() as isize + idx;
        // if negative, it means index was out of range
        if idx < 0 {
            return None;
        }

        self.get(idx as usize)
    }
}

macro_rules! impl_get_for_bounds {
    ($r:ty) => {
        impl<T: RealScalar> Get<$r, T> for Poly<T> {
            fn get(&self, idx: $r) -> Option<Poly<T>> {
                self.get_range_inner(idx.start_bound(), idx.end_bound())
            }
        }
    };
}

impl_get_for_bounds!(Range<usize>);
impl_get_for_bounds!(RangeInclusive<usize>);
impl_get_for_bounds!(RangeFrom<usize>);
impl_get_for_bounds!(RangeTo<usize>);
impl_get_for_bounds!(RangeToInclusive<usize>);
impl_get_for_bounds!(RangeFull);

// TODO: should be defined in terms of Get or vice-versa
impl<T: RealScalar> Index<usize> for Poly<T> {
    type Output = Complex<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_get_isize() {
        let p = poly![0.0, 1.0, 2.0, 3.0];
        assert_eq!(p.get(0isize).unwrap(), poly![0.0]);
        assert_eq!(p.get(1isize).unwrap(), poly![0.0, 1.0]);
        assert_eq!(p.get(-1isize).unwrap(), poly![0.0, 0.0, 0.0, 3.0]);
        assert_eq!(p.get(-2isize).unwrap(), poly![0.0, 0.0, 2.0]);
        assert_eq!(p.get(-4isize).unwrap(), poly![0.0]);

        // out of bounds
        assert!(p.get(4isize).is_none());
        assert!(p.get(-5isize).is_none());
    }

    #[test]
    fn test_get_usize() {
        let p = poly![0.0, 1.0, 2.0, 3.0];
        assert_eq!(p.get(0usize).unwrap(), poly![0.0]);
        assert_eq!(p.get(1usize).unwrap(), poly![0.0, 1.0]);
        assert_eq!(p.get(3usize).unwrap(), poly![0.0, 0.0, 0.0, 3.0]);

        // out of bounds
        assert!(p.get(4usize).is_none());
    }

    #[test]
    fn get_range() {
        let p = poly![1.0, 2.0, 3.0, 4.0];
        assert_eq!(p.get(1..3), Some(poly![0.0, 2.0, 3.0]));
    }

    #[test]
    fn get_range_inclusive() {
        let p = poly![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(p.get(2..=4), p.get(2..5));
    }

    #[test]
    fn get_range_full() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(p.get(..), Some(p));
    }
}
