use anyhow::ensure;
use na::{ComplexField, DMatrix, DVector, RealField};

use crate::{
    num::{Complex, Float, One, Zero},
    Poly, Scalar, ScalarOps,
    __util::linalg::set_subdiagonal,
};

use super::RootFinder;

mod francis_qr;
pub use francis_qr::FrancisQR;

/// What shape of companion matrix to use
///
/// Frobenius's companion matrix is the classic companion matrix $C$, of the form:
/// ```python
/// [
///     [0, 0, ..., 0, -c_0],
///     [1, 0, ..., 0, -c_1],
///     [0, 1, ..., 0, -c_2],
///     [:, :, ..., :,   : ],
///     [0, 0, ..., 1, -c_n]
/// ]
/// ```
///
/// The trasposed form is simply $C^T$.
///
/// The rotated form is a 180 degree rotation (not mirroring) of $C$.
///
/// Schmeisser is a companion matrix proposed by [Gerhard Schmeisser 1993](https://doi.org/10.1016/0024-3795(93)90268-S).
pub enum CompoanionMatrixType {
    Frobenius,
    FrobeniusTransposed,
    FrobeniusRotated,
    Schmeisser,
}

pub struct EigenState<T: Scalar> {
    matrix: Option<DMatrix<Complex<T>>>,
}

impl<T: Scalar> EigenState<T> {
    fn new() -> Self {
        Self { matrix: None }
    }
}

pub struct EigenConfig {
    companion_matrix_type: CompoanionMatrixType,
}

impl EigenConfig {
    fn new() -> Self {
        Self {
            companion_matrix_type: CompoanionMatrixType::Frobenius,
        }
    }
}

pub trait EigenvalueRootFinder<T: ScalarOps + RealField>: RootFinder<T> {
    fn with_companion_matrix_type(mut self, matrix_type: CompoanionMatrixType) -> Self {
        self.eigen_config_mut().companion_matrix_type = matrix_type;
        self
    }

    fn eigen_state_mut(&mut self) -> &mut EigenState<T>;

    fn eigen_state(&self) -> &EigenState<T>;

    fn eigen_config_mut(&mut self) -> &mut EigenConfig;

    fn eigen_config(&self) -> &EigenConfig;

    fn init_matrix(&mut self) -> anyhow::Result<()> {
        match self.eigen_config().companion_matrix_type {
            CompoanionMatrixType::Frobenius => {
                self.state_mut().poly.make_monic();
                self.eigen_state_mut().matrix = Some(self.state().poly.companion());
            }
            CompoanionMatrixType::FrobeniusTransposed => {
                self.state_mut().poly.make_monic();
                self.eigen_state_mut().matrix = Some(self.state().poly.companion());
                self.eigen_state_mut()
                    .matrix
                    .as_mut()
                    .map(|m: &mut DMatrix<_>| m.transpose_mut());
            }
            CompoanionMatrixType::FrobeniusRotated => {
                self.state_mut().poly.make_monic();
                self.eigen_state_mut().matrix = Some(self.state().poly.companion());
                let n = self
                    .eigen_state()
                    .matrix
                    .as_ref()
                    .expect("infallible")
                    .nrows();
                for i in 0..n / 2 {
                    self.eigen_state_mut()
                        .matrix
                        .as_mut()
                        .expect("infallible")
                        .swap_rows(i, n - i - 1);
                    self.eigen_state_mut()
                        .matrix
                        .as_mut()
                        .expect("infallible")
                        .swap_columns(i, n - i - 1);
                }
            }
            CompoanionMatrixType::Schmeisser => {
                self.eigen_state_mut().matrix = Some(schmeisser(&self.state().poly)?);
            }
        }
        Ok(())
    }
}

/// The modified Euclidean algorithm proposed by [Schmeisser 1993](https://doi.org/10.1016/0024-3795(93)90268-S).
fn schmeisser<T: ScalarOps + Float + RealField>(
    poly: &Poly<T>,
) -> anyhow::Result<DMatrix<Complex<T>>> {
    debug_assert!(poly.is_normalized());
    let mut u = poly.clone();
    u.make_monic();
    let n = u.degree_raw();

    let f1 = u.clone();
    let f2 = u
        .diff()
        .scaled(Complex::new(T::from_usize(n).expect("overflow"), T::zero()).recip());

    let mut f_i = f1;
    let mut f_ip1 = f2;
    let mut cs = vec![];
    let mut qs = vec![];

    for _ in 1..n {
        let (q, r) = f_i
            .clone()
            .div_rem(&f_ip1)
            .expect("zero divisor should have been taken care of by now");
        //let r = -r;

        if r.is_zero() {
            cs.push(Complex::zero());
            let mut f_ip1_diff = f_ip1.clone().diff();
            f_ip1_diff.make_monic();
            f_i = f_ip1;
            f_ip1 = f_ip1_diff;
        } else {
            cs.push((-r.last()).sqrt());
            let mut f_ip2 = r;
            f_ip2.make_monic();
            f_i = f_ip1;
            f_ip1 = f_ip2;
        }

        qs.push(-q.eval_point(Complex::zero()))
    }

    ensure!(
        f_ip1.is_one(),
        "cannot compute Schmeisser companion matrix due to pathological input"
    );

    qs.push(-f_i.eval_point(Complex::zero()));

    let mut t = DMatrix::from_diagonal(&DVector::from(qs));
    set_subdiagonal(&mut t.as_view_mut(), 1, &cs);
    set_subdiagonal(&mut t.as_view_mut(), -1, &cs);

    Ok(t)
}

#[cfg(test)]
mod test {
    use na::DMatrix;

    use crate::{Poly, __util::linalg::set_subdiagonal, num::Complex};

    use super::schmeisser;

    #[test]
    fn test_schmeisser() {
        // based on [Schmeisser 1993](https://doi.org/10.1016/0024-3795(93)90268-S).
        let p = Poly::from_roots(&[
            complex!(0.0, 0.0),
            complex!(1.0, 0.0),
            complex!(-1.0, 0.0),
            complex!(2.0, 0.0),
            complex!(-2.0, 0.0),
        ]);
        let mut expected = DMatrix::zeros(5, 5);
        let cs = &[
            std::f64::consts::SQRT_2,
            (7.0f64 / 5.0).sqrt(),
            (36.0f64 / 35.0).sqrt(),
            (4.0f64 / 7.0).sqrt(),
        ];
        set_subdiagonal(&mut expected.as_view_mut(), -1, cs);
        set_subdiagonal(&mut expected.as_view_mut(), 1, cs);
        let expected = expected.cast::<Complex<_>>();
        let got = schmeisser(&p).unwrap();
        let diff = &expected - &got;
        assert!(diff.sum().norm() < 1E-14);
    }
}
