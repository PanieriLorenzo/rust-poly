use crate::num::{Complex, Float};
use na::RealField;

use crate::ScalarOps;

use super::RootFinder;

mod naive;
pub use naive::Naive;
mod newton;
pub use newton::Newton;

pub enum DeflationStrategy {
    LongDivision,
    DeflateForward,
    DeflateBackward,
    DeflateComposite,
}

#[allow(clippy::module_name_repetitions)]
pub trait IterativeRootFinder<T: ScalarOps + PartialOrd + Float + RealField>:
    RootFinder<T>
{
    /// Find one root, without shrinkage
    ///
    /// # Errors
    /// Finder did not converge
    fn next_root(&mut self) -> super::Result<T>;

    fn with_deflation_strategy(&mut self, _strat: DeflationStrategy) {
        todo!()
    }

    /// Find multiple roots by shrinking
    ///
    /// # Errors
    /// Finder did not converge
    ///
    /// # Panics
    /// If `n` is larger than the degree of the polynomial
    fn next_n_roots(&mut self, n: usize) -> super::Result<T> {
        debug_assert!(self.state().poly.is_normalized());
        assert!(
            n <= self.state().poly.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

        for i in 0..n {
            if self.state().poly.degree_raw() == 0 {
                // we found all the roots early
                break;
            }
            if let Some(history_handle) = self.history() {
                history_handle.roots_history.push(vec![]);
            }
            let rs = self.next_root()?;
            self.state_mut().clean_roots.extend(rs.iter().copied());
            if i != (n - 1) {
                // if polynomial was second degree, next_root returns two roots
                for r in rs {
                    self.state_mut().poly = self.state_mut().poly.clone().deflate_composite(r);
                }
            }
        }

        Ok(self.state_mut().clean_roots.clone())
    }

    /// Garwick & Ward stopping criterion (see [Nikolajsen 2014](https://doi.org/10.1098/rsos.140206))
    // TODO: with specialization use Nikolajsen 2014 if T is f64 or f32, but right
    //       now this is fine for all real-like, including fractions and infinite
    //       precision
    fn stop_iteration(&self, z: Complex<T>, z_old: Complex<T>, z_old_old: Complex<T>) -> bool {
        let delta_z = (z - z_old).norm();
        let delta_z_old = (z_old - z_old_old).norm();
        let z_norm = z_old.norm();
        let em3 = T::from_f64(1E-3).expect("overflow");
        let em4 = T::from_f64(1E-4).expect("overflow");
        let em7 = T::from_f64(1E-7).expect("overflow");
        (z_norm < em4 && delta_z <= em7 || z_norm >= em4 && delta_z / z_norm <= em3)
            && delta_z >= delta_z_old
    }
}
