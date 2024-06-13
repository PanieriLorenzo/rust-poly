use crate::{
    num::{Complex, Float, Zero},
    Poly, ScalarOps,
    __util::complex::{c_min, c_neg},
};

impl<T: ScalarOps + Float> Poly<T> {
    /// [ref](https://doi.org/10.1007/BF01933524)
    pub(crate) fn initial_guess_smallest(&self) -> Complex<T> {
        debug_assert!(self.is_normalized());
        debug_assert!(self.len_raw() >= 2);

        let small = Float::recip(T::from_u16(1_000).expect("overflow"));
        let p_diff = self.clone().diff();
        let mut pz = self.eval_point(Complex::zero());
        let mut pdz = p_diff.eval_point(Complex::zero());

        // avoid divide by zero
        if pdz.norm() < small {
            pz += small;
            pdz += small;
        }

        let theta = (c_neg(pz) / pdz).arg();
        let mut iter_coeffs = self.0.iter();
        let a0 = iter_coeffs.next().expect("infallible");

        let mut guess = iter_coeffs
            .zip(1..)
            .map(|(ak, k)| {
                Complex::i()
                    .scale(theta)
                    .exp()
                    .scale((a0 / ak).norm())
                    .powf(T::one() / T::from_usize(k).expect("overflow"))
            })
            .reduce(c_min)
            .expect("infallible")
            .scale(Float::recip(T::from_u8(2).expect("overflow")));

        if guess.im.is_zero() {
            // add a small constant because some methods can't converge to
            // complex roots if the initial guess is real
            guess += Complex::i().scale(Float::recip(T::from_u16(1_000).expect("overflow")));
        }
        guess
    }

    pub(crate) fn initial_guess_lower_bound(&self) -> T {
        todo!()
    }

    pub(crate) fn initial_guess_uppwer_bound(&self) -> T {
        todo!()
    }
}
