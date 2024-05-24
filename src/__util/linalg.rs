use std::{
    cmp,
    hint::black_box,
    ops::{Div, MulAssign},
};

use na::{
    dvector,
    givens::{self, GivensRotation},
    matrix, ComplexField, DMatrix, DMatrixView, DMatrixViewMut, DVectorView, Dim, Matrix2, OMatrix,
    RealField,
};
use nalgebra::DVector;
use num::{complex::Complex64, Complex, One, Zero};

use crate::Scalar;

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
                sum = sum.clone() + input[k as usize].clone() * kernel[j].clone();
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
        col.len() > 0,
        "can't construct householder matrix from empty vector"
    );

    let denom = col[0].clone() + col[0].clone().signum() * col.dot(&col).sqrt();

    debug_assert_ne!(
        denom,
        C::zero(),
        "can't construct householder matrix from column whose first entry is zero"
    );

    col.apply(|z| *z = z.clone().div(denom.clone()));

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
        col.len() > 0,
        "can't construct householder matrix from empty vector"
    );

    let denom = col[0].clone() + col[0].clone().signum() * col.dot(&col).sqrt();

    debug_assert_ne!(
        denom,
        Complex::zero(),
        "can't construct householder matrix from column whose first entry is zero"
    );

    col.apply(|z| *z = z.clone().div(denom.clone()));

    // ensure first element is 1
    col[0] = Complex::one();

    let norm = col.norm();
    col.apply(|z| *z = z.clone().div(norm.clone()))
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
            let r_mod = x.clone().hypot(y.clone());
            let c = Complex::from_real(x.clone().abs() / r_mod.clone());
            let s = sign(x) * y.conj() / r_mod;
            (c, s)
        }
    };

    matrix![cos.clone(), sin; -sin.conj(), cos.conj()]
}

fn get_diag<T: Clone>(this: DMatrixView<T>, idx: usize) -> T {
    this.row(idx)[idx].clone()
}

fn set_diag<T>(mut this: DMatrixViewMut<T>, idx: usize, val: T) {
    this.row_mut(idx)[idx] = val;
}

// TODO: document what LAPACK routine this corresponds to
// TODO: this should take a view not a reference
/// Port of [`rulinalg::matrix::decomposition::eigen::balance_matrix`](https://github.com/AtheMathmo/rulinalg/blob/0ea49678d2dfa0e0a0df9cd26f49f6330aef80c4/src/matrix/decomposition/eigen.rs#L12)
fn balance_matrix<T: Scalar + RealField>(mut this: DMatrixViewMut<Complex<T>>) {
    // TODO: lots of single letter variables, make it more readable

    let n = this.nrows();
    let radix = T::one() + T::one();

    debug_assert_eq!(n, this.ncols(), "matrix must be square");

    // TODO: store in workspace
    let mut d = DMatrix::<Complex<T>>::identity(n, n);
    let mut converged = false;

    // TODO: max iter
    while !converged {
        converged = true;

        for i in 0..n {
            let mut c: T = this.column(i).norm();
            let mut r: T = this.row(i).norm();

            let s = c.clone() * c.clone() + r.clone() * r.clone();
            let mut f = T::one();

            // TODO: max iter
            while c < r.clone() / radix.clone() {
                c = c * radix.clone();
                r = r / radix.clone();
                f = f * radix.clone();
            }

            // TODO: max iter
            while c >= r.clone() * radix.clone() {
                c = c / radix.clone();
                r = r * radix.clone();
                f = f / radix.clone();
            }

            if (c.clone() * c.clone() + r.clone() * r.clone())
                < T::from_f64(0.95).expect("infallible") * s
            {
                converged = false;
                d.row_mut(i)[i] = d.row(i)[i].clone().scale(f.clone());

                // TODO: separate this to a separate procedure with a descriptive name
                for j in 0..n {
                    this.row_mut(j)[i] = this.row(j)[i].clone().scale(f.clone());
                    this.row_mut(i)[j] = this.row(i)[j].clone().scale(T::one() / f.clone());
                }
            }
        }
    }
}

/// Francis shift algorithm
pub(crate) fn eigen_francis_shift<T: Scalar + RealField>(
    mut this: DMatrixViewMut<Complex<T>>,
    epsilon: T,
) -> Result<Vec<Complex<T>>, ()> {
    let n = this.nrows();

    debug_assert!(n > 2, "use eigen_2x2 for 2x2 matrices");
    debug_assert_eq!(n, this.ncols(), "matrix must be square");

    // TODO: optional step, i.e. make a tuning option
    upper_hessenberg(this.as_view_mut());
    let mut h = this;
    let _h_dbg: Vec<_> = h.iter().collect();
    black_box(_h_dbg);
    balance_matrix(h.as_view_mut());
    let _h_dbg: Vec<_> = h.iter().collect();
    black_box(_h_dbg);

    // the final index of the active matrix (???)
    let mut p = n - 1;

    // TODO: maxiter
    let mut niter = 0;
    while p > 1 {
        if niter >= 100 {
            panic!();
        }
        niter += 1;
        // TODO: what are these?
        let q = p - 1;
        let s = h[(q, q)].clone() + h[(p, p)].clone();
        let t = h[(q, q)].clone() * h[(p, p)].clone() - h[(q, p)].clone() * h[(p, q)].clone();

        // TODO: what are these?
        let mut x = h.row(0)[0].clone() * h.row(0)[0].clone()
            + h.row(0)[1].clone() * h.row(1)[0].clone()
            - h.row(0)[0].clone() * s.clone()
            + t;
        let mut y = h.row(1)[0].clone() * (h.row(0)[0].clone() + h.row(1)[1].clone() - s);
        let mut z = h.row(1)[0].clone() * h.row(2)[1].clone();

        for k in 0..p - 1 {
            let householder = col_2_householder(dvector![x.clone(), y.clone(), z.clone()]);

            // apply householder transformation to block (on the left)
            {
                let r = 1.max(k) - 1;
                let mut h_block = h.view_mut((k, r), (3, n - r));
                let transformed = &householder * &h_block;
                // TODO: this could avoid an unnecessary copy if nalgebra was
                //       able to do in-place matmul
                h_block.copy_from(&transformed);
            }
            let _h_dbg: Vec<_> = h.iter().collect();
            black_box(_h_dbg);

            // apply householder transformation to block (on the right)
            {
                let r = cmp::min(k + 4, p + 1);
                let mut h_block = h.view_mut((0, k), (r, 3));
                let transformed = &h_block * householder.transpose();
                h_block.copy_from(&transformed);
            }
            let _h_dbg: Vec<_> = h.iter().collect();
            black_box(_h_dbg);

            x = h.row(k + 1)[k].clone();
            y = h.row(k + 2)[k].clone();

            if k < p - 2 {
                z = h.row(k + 3)[k].clone();
            }
        }

        let givens_mat = complex_givens_rot_2d(x, y);

        {
            // apply givens rotation to block (on the left)
            let mut h_block = h.view_mut((q, p - 2), (2, n - p + 2));
            let transformed = &givens_mat * &h_block;
            h_block.copy_from(&transformed);
        }
        let _h_dbg: Vec<_> = h.iter().collect();
        black_box(_h_dbg);

        {
            // apply givens rotation to block (on the right)
            let mut h_block = h.view_mut((0, q), (p + 1, 2));
            let transformed = &h_block * givens_mat.transpose();
            h_block.copy_from(&transformed);
        }
        let _h_dbg: Vec<_> = h.iter().collect();
        black_box(_h_dbg);

        // check for convergence
        if h[(p, q)].clone().abs()
            < epsilon.clone() * (h[(q, q)].clone().abs() + h[(p, p)].clone().abs())
        {
            // if h[p,q] < epsilon we flush it to zero and shrink the problem
            h[(p, q)] = Complex::<T>::zero();
            p -= 1;
        } else if h[(p - 1, q - 1)].clone().abs()
            < epsilon.clone() * (h[(q - 1, q - 1)].clone().abs() + h[(q, q)].clone().abs())
        {
            h[(p - 1, q - 1)] = Complex::<T>::zero();
            p -= 2;
        }
        let _h_dbg: Vec<_> = h.iter().collect();
        black_box(_h_dbg);
    }
    let _h_dbg: Vec<_> = h.iter().collect();
    black_box(_h_dbg);

    Ok(h.diagonal().as_slice().to_vec())
}

#[cfg(test)]
mod test {
    use na::{dmatrix, matrix};
    use num::complex::{Complex64, ComplexFloat};

    use super::{balance_matrix, eigen_francis_shift};

    #[test]
    fn test_balance_matrix() {
        let mut m = dmatrix![0.0, 1.0, 2.0; 3.0, 4.0, 5.0; 6.0, 7.0, 8.0].cast::<Complex64>();
        balance_matrix(m.as_view_mut());
        assert_eq!(
            m,
            dmatrix![0.0, 2.0, 4.0; 1.5, 4.0, 5.0; 3.0, 7.0, 8.0].cast::<Complex64>()
        );
    }

    #[test]
    fn test_egien_3x3() {
        let mut m =
            dmatrix![17.0, 22.0, 27.0; 22.0, 29.0, 36.0; 27.0, 36.0, 45.0].cast::<Complex64>();
        let eigs = eigen_francis_shift(m.as_view_mut(), 1E-9).unwrap();
        let eig_1 = 90.4026;
        let eig_2 = 0.5973;
        let eig_3 = 0.0;

        assert!(eigs.iter().any(|x| (x - eig_1).abs() < 1e-4));
        assert!(eigs.iter().any(|x| (x - eig_2).abs() < 1e-4));
        assert!(eigs.iter().any(|x| (x - eig_3).abs() < 1e-4));
    }
}
