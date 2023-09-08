use super::{Complex, Poly, Scalar};

pub trait Get<I, T: Scalar> {
    fn get(&self, idx: I) -> Option<Complex<T>>;
}

impl<T: Scalar> Get<usize, T> for Poly<T> {
    fn get(&self, idx: usize) -> Option<Complex<T>> {
        debug_assert!(self.is_normalized());

        if idx >= self.len_raw() {
            return None;
        }

        Some(self.0[idx].clone())
    }
}

impl<T: Scalar> Get<isize, T> for Poly<T> {
    fn get(&self, idx: isize) -> Option<Complex<T>> {
        if idx >= 0 {
            return self.get(idx as usize);
        }

        // find the index from the end
        let idx = self.len_raw() as isize + idx;
        // if negative, it means index was out of range
        if idx < 0 {
            return None;
        }

        return self.get(idx as usize);
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
        assert!(p.get(-5isize).is_none())
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
}
