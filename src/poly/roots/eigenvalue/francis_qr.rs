use crate::{
    num::{Float, One},
    roots::{Error::NoConverge, FinderConfig, FinderHistory, FinderState, RootFinder},
    Poly, Scalar, ScalarOps,
    __util::linalg::eigen_francis_shift,
};
use na::RealField;

use super::{EigenConfig, EigenState, EigenvalueRootFinder};

pub struct FrancisQR<T: Scalar> {
    state: FinderState<T>,
    config: FinderConfig<T>,
    history: Option<FinderHistory<T>>,
    eigen_config: EigenConfig,
    eigen_state: EigenState<T>,
}

impl<T: ScalarOps + Float + RealField> RootFinder<T> for FrancisQR<T> {
    fn from_poly(poly: crate::Poly<T>) -> Self {
        FrancisQR {
            state: FinderState::new(poly),
            config: FinderConfig::new(),
            history: None,
            eigen_config: EigenConfig::new(),
            eigen_state: EigenState::new(),
        }
    }

    fn roots(&mut self) -> crate::roots::Result<T> {
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
            Ok(v) => {
                self.state_mut().clean_roots.extend(v.iter());
                self.state_mut().poly = Poly::one();
                Ok(v)
            }
            Err(v) => {
                self.state_mut().dirty_roots.extend(v.iter());
                Err(NoConverge(v))
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
    use crate::{
        roots::{
            eigenvalue::{CompoanionMatrixType, EigenvalueRootFinder},
            RootFinder,
        },
        Poly,
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
            .with_companion_matrix_type(CompoanionMatrixType::Schmeisser)
            .roots()
            .unwrap();
        assert!(check_roots(roots, expected_roots, 1E-14));
    }
}
