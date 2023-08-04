use nalgebra::DVector;
use num_complex::Complex;
use num_traits::Zero;

use crate::Scalar;

pub fn convolve_1d<T: Scalar>(
    input: &DVector<Complex<T>>,
    kernel: &DVector<Complex<T>>,
) -> DVector<Complex<T>> {
    let input_len = input.len();
    let kernel_len = kernel.len();
    let output_len = input_len + kernel_len - 1;

    let mut output: DVector<Complex<T>> = DVector::<Complex<T>>::zeros(output_len);

    for i in 0..output_len {
        let mut sum = Complex::<T>::zero();
        for j in 0..kernel_len {
            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            let k = i as isize - j as isize;

            // will only wrap with polynomials so large they don't fit in memory
            #[allow(clippy::cast_possible_wrap)]
            // k is guaranteed to be positive by the conditional
            #[allow(clippy::cast_sign_loss)]
            if k >= 0 && k < input_len as isize {
                sum = sum.clone() + input[k as usize].clone() * kernel[j].clone();
            }
        }
        output[i] = sum;
    }
    output
}
