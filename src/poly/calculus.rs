use itertools::chain;
use num::{Complex, Zero};

use crate::{Poly, Scalar, __util::casting::usize_to_scalar};

impl<T: Scalar> Poly<T> {
    /// Derivative
    #[must_use]
    pub fn diff(self) -> Self {
        debug_assert!(self.is_normalized());

        // derivative of constant is zero
        if self.degree_raw() == 0 {
            return Self::from_real_slice(&[T::zero()]);
        }

        let coeffs: Vec<_> = (0..self.len())
            .map(|x| T::from_usize(x).expect("degree to high to convert to T"))
            .map(|x| Complex::new(x, T::zero()))
            .zip(self.0.iter())
            .map(|(n, c)| n * c)
            .skip(1) // shift degrees down
            .collect();
        Self::from_complex_vec(coeffs).normalize()
    }

    /// Antiderivative (with C=0)
    pub fn integral(self) -> Self {
        debug_assert!(self.is_normalized());

        let coeffs: Vec<_> = chain(
            [Complex::<T>::zero()],
            (1..=self.len())
                .map(usize_to_scalar::<T>)
                .map(|x| Complex::new(x, T::zero()))
                .zip(self.0.iter())
                .map(|(n, c)| c / n),
        )
        .collect();
        Self::from_complex_vec(coeffs).normalize()
    }
}

#[cfg(test)]
mod test {
    use crate::Poly64;

    #[test]
    fn diff() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(p.diff(), poly![2.0, 6.0]);
    }

    /// This was a bug
    #[test]
    fn diff1() {
        let one = poly![1.0];
        assert_eq!(one.diff().degree(), 0);
    }

    #[test]
    fn integral() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(p.integral(), poly![0.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn integral_diff() {
        let p = poly![1.0, 2.0, 3.0];
        let q = p.clone().integral().diff();
        assert_eq!(p, q);
    }
}
