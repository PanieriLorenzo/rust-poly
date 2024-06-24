use crate::{
    num::{Complex, Float, One, Zero},
    roots::{
        self, eigenvalue::schmeisser, Error::NoConverge, FinderConfig, FinderHistory, FinderState,
        RootFinder,
    },
    Poly, Scalar, ScalarOps,
    __util::linalg::eigen_francis_shift,
};
use na::RealField;

use super::{CompanionMatrixType, EigenConfig, EigenState, EigenvalueRootFinder};

/// Find all roots of a polynomial using a naive Newton-Raphson approach.
///
/// You are likely looking for [`roots::newton`] instead, which is a more
/// robust implementation of the same algorithm.
///
/// This algorithm is not very stable, and will often fail to find fairly ordinary
/// roots. So only use this approach for specific purposes, namely:
/// - You know its not going to get stuck and you want a slight performance improvement
/// - You are benchmarking your custom root finder against a classical Newton-Raphson approach.
pub fn francis_qr<T: ScalarOps + RealField>(
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    companion_matrix_type: CompanionMatrixType,
) -> roots::Result<T> {
    let mut eval_counter = 0;
    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let mut roots = vec![];

    let trivial_roots = poly.trivial_roots(epsilon);
    eval_counter += trivial_roots.1;
    roots.extend(trivial_roots.0.iter());

    debug_assert!(poly.is_normalized());
    if poly.degree_raw() == 0 {
        log::debug!("{{evaluations: {eval_counter}}}");
        return Ok(roots);
    }

    let small_num = T::small_safe();
    let mut needs_unshifting = false;

    // leading zero coefficients can be problematic, so shifting to
    // avoid (the householder reflectors cannot be constructed if the
    // first element is zero)
    if poly.0[0].norm() < small_num {
        needs_unshifting = true;
        *poly = poly.clone().translate(Complex::one(), Complex::zero());
    }

    poly.make_monic();

    let mut companion = match companion_matrix_type {
        // TODO: companion shouldn't be a method
        CompanionMatrixType::Frobenius => poly.companion(),
        CompanionMatrixType::FrobeniusTransposed => poly.companion().transpose(),
        CompanionMatrixType::FrobeniusRotated => {
            let mut m = poly.companion();
            let n = m.nrows();
            for i in 0..n / 2 {
                m.swap_rows(i, n - i - 1);
                m.swap_columns(i, n - i - 1);
            }
            m
        }
        CompanionMatrixType::Schmeisser => schmeisser(poly)?,
    };

    match eigen_francis_shift(
        companion.as_view_mut(),
        epsilon,
        usize::MAX,
        max_iter.unwrap_or(usize::MAX),
    ) {
        Ok(mut v) => {
            if needs_unshifting {
                for r in &mut v {
                    *r -= Complex::one();
                }
            }
            *poly = Poly::one();
            roots.extend(v);
            log::debug!("{{evaluations: {eval_counter}}}");
            Ok(roots)
        }
        Err(mut v) => {
            if needs_unshifting {
                for r in &mut v {
                    *r -= Complex::one();
                }
                *poly = poly.clone().translate(-Complex::one(), Complex::zero());
            }
            // note that the roots in `roots` are clean, so we don't return
            // them in the error result
            Err(NoConverge(v))
        }
    }
}

pub struct FrancisQR<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    history: Option<FinderHistory<T>>,
    eigen_config: EigenConfig,
    eigen_state: EigenState<T>,
}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for FrancisQR<T> {
    fn from_poly(poly: crate::Poly<T>) -> Self {
        Self {
            state: FinderState::new(poly),
            config: FinderConfig::new(),
            history: None,
            eigen_config: EigenConfig::new(),
            eigen_state: EigenState::new(),
        }
    }

    fn run(&mut self) -> std::result::Result<(), crate::roots::Error<()>> {
        debug_assert!(self.state().poly.is_normalized());

        // handle trivial cases
        let epsilon = self.config().epsilon;
        let (mut roots, _) = self.state_mut().poly.trivial_roots(epsilon);
        if self.state().poly.degree_raw() == 0 {
            self.state_mut().clean_roots.extend(roots.iter());
            return Ok(());
        }

        let small_num = T::small_safe();
        let mut needs_unshifting = false;
        // leading zero coefficients can be problematic, so shifting to
        // avoid (the householder reflectors cannot be constructed if the
        // first element is zero)
        if self.state().poly.0[0].norm() < small_num {
            needs_unshifting = true;
            self.state_mut().poly = self
                .state()
                .poly
                .clone()
                .translate(Complex::one(), Complex::zero());
        }

        self.init_matrix()?;
        let epsilon = self.config().epsilon;
        let max_iter = self.config().max_iter;
        match eigen_francis_shift(
            self.eigen_state_mut()
                .matrix
                .as_mut()
                .expect("infallible")
                .as_view_mut(),
            epsilon,
            usize::MAX,
            max_iter,
        ) {
            Ok(mut v) => {
                if needs_unshifting {
                    for r in &mut v {
                        *r -= Complex::one();
                    }
                }
                self.state_mut().clean_roots.extend(v.iter());
                self.state_mut().poly = Poly::one();
                roots.extend(v);
                Ok(())
            }
            Err(mut v) => {
                if needs_unshifting {
                    for r in &mut v {
                        *r -= Complex::one();
                    }
                    self.state_mut().poly = self
                        .state()
                        .poly
                        .clone()
                        .translate(-Complex::one(), Complex::zero());
                }
                self.state_mut().dirty_roots.extend(v.iter());
                // note that the roots in `roots` are clean, so we don't return
                // them in the error result
                Err(NoConverge(()))
            }
        }
    }

    fn state_mut(&mut self) -> &mut crate::roots::FinderState<T> {
        &mut self.state
    }

    fn state(&self) -> &crate::roots::FinderState<T> {
        &self.state
    }

    fn config(&mut self) -> &mut crate::roots::FinderConfig<T> {
        &mut self.config
    }

    fn history(&mut self) -> &mut Option<crate::roots::FinderHistory<T>> {
        &mut self.history
    }
}

impl<T: ScalarOps + RealField + Float> EigenvalueRootFinder<T> for FrancisQR<T> {
    fn eigen_state_mut(&mut self) -> &mut EigenState<T> {
        &mut self.eigen_state
    }

    fn eigen_state(&self) -> &EigenState<T> {
        &self.eigen_state
    }

    fn eigen_config_mut(&mut self) -> &mut EigenConfig {
        &mut self.eigen_config
    }

    fn eigen_config(&self) -> &EigenConfig {
        &self.eigen_config
    }
}

#[cfg(test)]
mod test {
    use num::One;

    use crate::{
        roots::{
            eigenvalue::{CompanionMatrixType, EigenvalueRootFinder},
            RootFinder,
        },
        Poly, Poly64,
        __util::testing::check_roots,
    };

    use super::FrancisQR;

    #[test]
    fn francis_qr() {
        let expected_roots = vec![complex!(1.0), complex!(1.0), complex!(2.0), complex!(3.0)];
        let p = Poly::from_roots(&expected_roots);
        let roots = FrancisQR::from_poly(p)
            .with_epsilon(1E-14)
            .with_max_iter(100)
            .with_companion_matrix_type(CompanionMatrixType::Schmeisser)
            .roots()
            .unwrap();
        assert!(check_roots(roots, expected_roots, 1E-14));
    }

    #[test]
    pub fn degree_0() {
        let mut p = Poly64::one();
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(roots.is_empty());
        assert!(p.is_one());
    }

    #[test]
    fn degree_1() {
        let roots_expected = vec![complex!(1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_2() {
        let roots_expected = vec![complex!(1.0), complex!(2.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-12));
    }

    #[test]
    fn degree_3() {
        let roots_expected = vec![complex!(1.0), complex!(2.0), complex!(3.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-8));
    }

    #[test]
    fn degree_3_complex() {
        let roots_expected = vec![complex!(1.0), complex!(0.0, 1.0), complex!(0.0, -1.0)];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusTransposed,
        )
        .unwrap();
        assert!(check_roots(roots, roots_expected, 2.0));
    }

    #[test]
    fn degree_5_multiplicity_3() {
        let roots_expected = vec![
            complex!(1.0),
            complex!(2.0),
            complex!(2.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(
            check_roots(roots.clone(), roots_expected, 1E-2),
            "{roots:?}"
        );
    }

    #[test]
    fn degree_5_2_zeros() {
        let roots_expected = vec![
            complex!(0.0),
            complex!(0.0),
            complex!(1.0),
            complex!(2.0),
            complex!(3.0),
        ];
        let mut p = crate::Poly::from_roots(&roots_expected);
        let roots = super::francis_qr(
            &mut p,
            Some(1E-14),
            Some(100),
            CompanionMatrixType::FrobeniusRotated,
        )
        .unwrap();
        assert!(check_roots(roots, roots_expected, 1E-8));
    }
}
