use std::ops::{
    Bound, Index, Range, RangeBounds, RangeFrom, RangeFull, RangeInclusive, RangeTo,
    RangeToInclusive,
};

use itertools::Itertools;

use super::{Complex, Poly, Scalar};

mod sealed {
    pub trait Sealed {}
}

pub trait Get<I, T>: sealed::Sealed {
    type Output;

    fn get(&self, idx: I) -> Option<Self::Output>;
}

impl<T> sealed::Sealed for Poly<T> {}

use super::*;

impl<T: Scalar> Poly<T> {
    /// Implementation for all range-based indexing (because Rust is super annoying
    /// around the different range iterators, see [#3550](https://github.com/rust-lang/rfcs/pull/3550))
    fn get_range_inner(&self, idx_start: Bound<&usize>, idx_end: Bound<&usize>) -> Poly<T> {
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
        let coeffs = &self.as_slice()[start..end];
        Poly::from_complex_slice(coeffs).shift_up(start)
    }
}

impl<T: Scalar> Get<usize, T> for Poly<T> {
    type Output = Complex<T>;

    fn get(&self, idx: usize) -> Option<Complex<T>> {
        debug_assert!(self.is_normalized());

        if idx >= self.len_raw() {
            return None;
        }

        Some(self.0[idx].clone())
    }
}

impl<T: Scalar> Get<isize, T> for Poly<T> {
    type Output = Complex<T>;

    #[allow(clippy::cast_possible_wrap)]
    #[allow(clippy::cast_sign_loss)]
    fn get(&self, idx: isize) -> Option<Complex<T>> {
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

// impl<T: Scalar> Get<Range<usize>, T> for Poly<T> {
//     type Output = Poly<T>;

//     fn get(&self, idx: Range<usize>) -> Option<Poly<T>> {
//
//     }
// }

macro_rules! impl_get_for_bounds {
    ($r:ty) => {
        impl<T: Scalar> Get<$r, T> for Poly<T> {
            type Output = Poly<T>;

            fn get(&self, idx: $r) -> Option<Poly<T>> {
                Some(self.get_range_inner(idx.start_bound(), idx.end_bound()))
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
impl<T: Scalar> Index<usize> for Poly<T> {
    type Output = Complex<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[cfg(test)]
mod test {
    use num::{complex::Complex64, One, Zero};
    use numeric_constant_traits::{Three, Two};

    use super::*;

    #[test]
    fn test_get_isize() {
        let p = poly![0.0, 1.0, 2.0, 3.0];
        assert_eq!(p.get(0isize).unwrap(), Complex64::zero());
        assert_eq!(p.get(1isize).unwrap(), Complex64::one());
        assert_eq!(p.get(-1isize).unwrap(), Complex64::three());
        assert_eq!(p.get(-2isize).unwrap(), Complex64::two());
        assert_eq!(p.get(-4isize).unwrap(), Complex64::zero());

        // out of bounds
        assert!(p.get(4isize).is_none());
        assert!(p.get(-5isize).is_none());
    }

    #[test]
    fn test_get_usize() {
        let p = poly![0.0, 1.0, 2.0, 3.0];
        assert_eq!(p.get(0usize).unwrap(), Complex64::zero());
        assert_eq!(p.get(1usize).unwrap(), Complex64::one());
        assert_eq!(p.get(3usize).unwrap(), Complex64::three());

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
