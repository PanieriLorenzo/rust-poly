use std::{cmp, ops::Div};

use na::{dvector, matrix, ComplexField, DMatrix, DMatrixViewMut, Matrix2, RealField};
use nalgebra::DVector;
use num::{Complex, One, Zero};

use crate::Scalar;

pub(crate) fn set_subdiagonal<T: Clone>(matrix: &mut DMatrixViewMut<T>, k: i128, values: &[T]) {
    let m = matrix.nrows() as i128;
    let n = matrix.ncols() as i128;
    for (i, s) in ((-k).max(0)..(n - k).min(m)).enumerate() {
        matrix.row_mut(s as usize)[(s + k) as usize] = values[i].clone();
    }
}

pub(crate) fn convolve_1d<T: Scalar>(
    input: &DVector<Complex<T>>,
    kernel: &DVector<Complex<T>>,
) -> DVector<Complex<T>> {
    let input_len = input.len();
    let kernel_len = kernel.len();

    debug_assert!(input_len + kernel_len > 0);
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
                sum = sum + input[k as usize] * kernel[j];
            }
        }
        output[i] = sum;
    }
    output
}

/// Create a householder matrix from a column vector
fn col_2_householder<T: Scalar + RealField>(mut col: DVector<Complex<T>>) -> DMatrix<Complex<T>> {
    type C<T> = Complex<T>;

    debug_assert!(
        !col.is_empty(),
        "can't construct householder matrix from empty vector"
    );

    let denom = col[0] + col[0].signum() * col.dot(&col).sqrt();

    debug_assert_ne!(
        denom,
        C::zero(),
        "can't construct householder matrix from column whose first entry is zero"
    );

    col.apply(|z| *z = (*z).div(denom));

    // ensure first element is 1
    col[0] = C::one();

    let two = C::<T>::one() + C::<T>::one();
    let norm_sq = col.norm_squared();
    let conj_transpose = col.conjugate().transpose();
    DMatrix::identity(col.len(), col.len()) - (col * conj_transpose) * (two / norm_sq)
}

/// Crate a householder vector from a column vector
fn col_2_householder_vec<T: Scalar + RealField>(col: &mut DVector<Complex<T>>) {
    debug_assert!(
        !col.is_empty(),
        "can't construct householder vector from empty vector"
    );

    let denom = col[0] + col[0].signum() * col.dot(col).sqrt();

    debug_assert_ne!(
        denom,
        Complex::zero(),
        "can't construct householder vector from column whose first entry is zero"
    );

    col.apply(|z| *z = (*z).div(denom));

    // ensure first element is 1
    col[0] = Complex::one();

    let norm = col.norm();
    col.apply(|z| *z = (*z).div(norm));
}

fn upper_hessenberg<T: Scalar + RealField>(mut this: DMatrixViewMut<Complex<T>>) {
    let n = this.nrows();
    debug_assert_eq!(n, this.ncols(), "matrix must be square");

    for i in 0..n - 2 {
        let mut h_holder_vec: DVector<_> = this.view((i + 1, i), (n - i - 1, 1)).column(0).into();
        col_2_householder_vec(&mut h_holder_vec);

        {
            // apply holder to the left
            let two = Complex::<T>::one() + Complex::one();
            let mut block = this.view_mut((i + 1, i), (n - i - 1, n - i));
            block -= &h_holder_vec * (h_holder_vec.transpose() * &block) * two;
        }

        {
            // apply holder to the right
            let two = Complex::<T>::one() + Complex::one();
            let mut block = this.view_mut((0, i + 1), (n, n - i - 1));
            block -= (&block * &h_holder_vec) * h_holder_vec.transpose() * two;
        }
    }

    // flush lower left triangle to zero
    for i in 0..n - 2 {
        for j in i + 2..n {
            unsafe {
                *this.get_unchecked_mut((j, i)) = Complex::zero();
            }
        }
    }
}

/// Compute the complex Givens rotation matrix for 2D vectors
///
/// Ref: D. Bindel, J. Demmel, W. Kahan, O. Marques "On Computing Givens Rotations
/// Reliably and Efficiently" [DOI](https://doi.org/10.1145/567806.567809)
#[allow(clippy::similar_names)]
fn complex_givens_rot_2d<T: Scalar + RealField>(
    x: Complex<T>,
    y: Complex<T>,
) -> Matrix2<Complex<T>> {
    let sign = |z: Complex<T>| {
        if z.is_zero() {
            Complex::one()
        } else {
            z.signum()
        }
    };

    let (cos, sin) = match (x.is_zero(), y.is_zero()) {
        (_, true) => (Complex::one(), Complex::zero()),
        (true, false) => (Complex::zero(), sign(y.conj())),
        (false, false) => {
            let r_mod = x.hypot(y);
            let c = Complex::from_real(x.abs() / r_mod);
            let s = sign(x) * y.conj() / r_mod;
            (c, s)
        }
    };

    matrix![cos, sin; -sin.conj(), cos.conj()]
}

// TODO: document what LAPACK routine this corresponds to
// TODO: one-letter bindings
/// Port of [`rulinalg::matrix::decomposition::eigen::balance_matrix`](https://github.com/AtheMathmo/rulinalg/blob/0ea49678d2dfa0e0a0df9cd26f49f6330aef80c4/src/matrix/decomposition/eigen.rs#L12)
#[allow(clippy::many_single_char_names)]
fn balance_matrix<T: Scalar + RealField>(
    mut this: DMatrixViewMut<Complex<T>>,
    _max_iter: usize,
) -> Result<(), ()> {
    // TODO: lots of single letter variables, make it more readable
    const MAX_ITER_INNER: usize = 100;

    let n = this.nrows();
    let radix = T::one() + T::one();

    // prevent infinite loop
    let max_iter_outer = 100.max(n);

    debug_assert_eq!(n, this.ncols(), "matrix must be square");

    // TODO: store in workspace
    let mut d = DMatrix::<Complex<T>>::identity(n, n);
    let mut converged = false;

    for _ in 0..max_iter_outer {
        if converged {
            return Ok(());
        }

        converged = true;

        for i in 0..n {
            let mut c: T = this.column(i).norm();
            let mut r: T = this.row(i).norm();

            let s = c * c + r * r;
            let mut f = T::one();

            'for_else: {
                for _ in 0..MAX_ITER_INNER {
                    match c.partial_cmp(&(r / radix)) {
                        Some(cmp::Ordering::Equal | cmp::Ordering::Greater) => break 'for_else,
                        None => return Err(()),
                        _ => {}
                    }
                    c *= radix;
                    r /= radix;
                    f *= radix;
                }
                // else: did not converge
                return Err(());
            }

            'for_else: {
                for _ in 0..MAX_ITER_INNER {
                    match c.partial_cmp(&(r * radix)) {
                        Some(cmp::Ordering::Less) => break 'for_else,
                        None => return Err(()),
                        _ => {}
                    }
                    c /= radix;
                    r *= radix;
                    f /= radix;
                }
                // else: did not converge
                return Err(());
            }

            if (c * c + r * r) < T::from_f64(0.95).expect("infallible") * s {
                converged = false;
                d.row_mut(i)[i] = d.row(i)[i].scale(f);

                for j in 0..n {
                    this.row_mut(j)[i] = this.row(j)[i].scale(f);
                    this.row_mut(i)[j] = this.row(i)[j].scale(T::one() / f);
                }
            }
        }
    }
    Err(())
}

/// Francis shift algorithm
// TODO: single-char variables
#[allow(clippy::many_single_char_names)]
pub(crate) fn eigen_francis_shift<T: Scalar + RealField>(
    mut this: DMatrixViewMut<Complex<T>>,
    epsilon: T,
    _max_iter: usize,
    max_iter_per_deflation: usize,
) -> Result<Vec<Complex<T>>, Vec<Complex<T>>> {
    // TODO: tune this to input
    const MAX_ITER_BALANCE: usize = 100;

    let n = this.nrows();

    // prevent infinite loop
    let max_total_iter: usize = max_iter_per_deflation * n;

    debug_assert!(n > 2, "use eigen_2x2 for 2x2 matrices");
    debug_assert_eq!(n, this.ncols(), "matrix must be square");

    // TODO: optional step, i.e. make a tuning option
    upper_hessenberg(this.as_view_mut());
    let mut h = this;
    balance_matrix(h.as_view_mut(), MAX_ITER_BALANCE).map_err(|()| vec![])?;

    // the final index of the active sub-matrix
    let mut p = n - 1;

    // keeps track of total (inner) iterations per each deflation step
    let mut piter = 0;

    for _ in 0..max_total_iter {
        if p <= 1 {
            break;
        }

        if piter >= max_iter_per_deflation {
            // TODO: return eigenvalues so far
            return Err(vec![]);
        }

        // TODO: what is q?
        let q = p - 1;
        let s = h[(q, q)] + h[(p, p)];
        let t = h[(q, q)] * h[(p, p)] - h[(q, p)] * h[(p, q)];

        // TODO: what are these?
        let mut x = h.row(0)[0] * h.row(0)[0] + h.row(0)[1] * h.row(1)[0] - h.row(0)[0] * s + t;
        let mut y = h.row(1)[0] * (h.row(0)[0] + h.row(1)[1] - s);
        let mut z = h.row(1)[0] * h.row(2)[1];

        for k in 0..p - 1 {
            let householder = col_2_householder(dvector![x, y, z]);

            // apply householder transformation to block (on the left)
            {
                let r = 1.max(k) - 1;
                let mut h_block = h.view_mut((k, r), (3, n - r));
                let transformed = &householder * &h_block;
                // TODO: this could avoid an unnecessary copy if nalgebra was
                //       able to do in-place matmul
                h_block.copy_from(&transformed);
            }

            // apply householder transformation to block (on the right)
            {
                let r = cmp::min(k + 4, p + 1);
                let mut h_block = h.view_mut((0, k), (r, 3));
                let transformed = &h_block * householder.transpose();
                h_block.copy_from(&transformed);
            }

            x = h.row(k + 1)[k];
            y = h.row(k + 2)[k];

            if k < p - 2 {
                z = h.row(k + 3)[k];
            }
        }

        let givens_mat = complex_givens_rot_2d(x, y);

        {
            // apply givens rotation to block (on the left)
            let mut h_block = h.view_mut((q, p - 2), (2, n - p + 2));
            let transformed = givens_mat * &h_block;
            h_block.copy_from(&transformed);
        }

        {
            // apply givens rotation to block (on the right)
            let mut h_block = h.view_mut((0, q), (p + 1, 2));
            let transformed = &h_block * givens_mat.transpose();
            h_block.copy_from(&transformed);
        }

        // check for convergence
        if h[(p, q)].abs() < epsilon * (h[(q, q)].abs() + h[(p, p)].abs()) {
            // deflation step:
            // if h[p,q] < epsilon we flush it to zero and shrink the problem
            h[(p, q)] = Complex::<T>::zero();
            p -= 1;
            piter = 0;
        } else if h[(p - 1, q - 1)].abs() < epsilon * (h[(q - 1, q - 1)].abs() + h[(q, q)].abs()) {
            // deflation step for 2x2 blocks
            h[(p - 1, q - 1)] = Complex::<T>::zero();
            p -= 2;
            piter = 0;
        } else {
            // if no deflation happened, increment counter for current deflation
            piter += 1;
        }
    }

    Ok(h.diagonal().as_slice().to_vec())
}

#[cfg(test)]
mod test {
    use na::{dmatrix, DMatrix};
    use num::complex::{Complex64, ComplexFloat};

    use super::{balance_matrix, eigen_francis_shift, set_subdiagonal};

    #[test]
    fn test_balance_matrix() {
        let mut m = dmatrix![0.0, 1.0, 2.0; 3.0, 4.0, 5.0; 6.0, 7.0, 8.0].cast::<Complex64>();
        balance_matrix(m.as_view_mut(), 100).unwrap();
        assert_eq!(
            m,
            dmatrix![0.0, 2.0, 4.0; 1.5, 4.0, 5.0; 3.0, 7.0, 8.0].cast::<Complex64>()
        );
    }

    #[test]
    fn test_egien_3x3() {
        let mut m =
            dmatrix![17.0, 22.0, 27.0; 22.0, 29.0, 36.0; 27.0, 36.0, 45.0].cast::<Complex64>();
        let eigs = eigen_francis_shift(m.as_view_mut(), 1E-9, 100, 100).unwrap();
        let eig_1 = 90.4026;
        let eig_2 = 0.5973;
        let eig_3 = 0.0;

        assert!(eigs.iter().any(|x| (x - eig_1).abs() < 1e-4));
        assert!(eigs.iter().any(|x| (x - eig_2).abs() < 1e-4));
        assert!(eigs.iter().any(|x| (x - eig_3).abs() < 1e-4));
    }

    #[test]
    fn test_set_subdiagonal() {
        let mut m = DMatrix::zeros(3, 4);
        set_subdiagonal(&mut m.as_view_mut(), 1, &[1.0f64, 2.0, 3.0]);
        println!("{m}");
    }
}
