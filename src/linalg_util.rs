use nalgebra::DMatrix;
use num_complex::Complex;

use crate::Scalar;

pub fn reverse_mut<T: Scalar>(mat: &mut DMatrix<Complex<T>>) {
    let n = mat.shape().0;
    for i in 0..n / 2 {
        mat.swap_rows(i, n - i - 1);
    }
    let n = mat.shape().1;
    for i in 0..n / 2 {
        mat.swap_columns(i, n - i - 1);
    }
}
