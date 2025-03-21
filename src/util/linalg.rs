use num::{Complex, Zero};

use crate::RealScalar;

pub(crate) fn convolve_1d<T: RealScalar>(
    input: &[Complex<T>],
    kernel: &[Complex<T>],
) -> Vec<Complex<T>> {
    let input_len = input.len();
    let kernel_len = kernel.len();

    debug_assert!(input_len + kernel_len > 0);
    let output_len = input_len + kernel_len - 1;

    // TODO: could probably use a collect_vec to make it more idiomatic
    let mut output = vec![Complex::zero(); output_len];

    for (i, o) in output.iter_mut().enumerate() {
        let mut sum = Complex::<T>::zero();
        for (j, ker) in kernel.iter().enumerate() {
            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            let k = i as isize - j as isize;

            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            // k is guaranteed to be positive by the conditional
            #[allow(clippy::cast_sign_loss)]
            if k >= 0 && k < input_len as isize {
                sum += input[k as usize].clone() * ker.clone();
            }
        }
        *o = sum;
    }
    output
}
