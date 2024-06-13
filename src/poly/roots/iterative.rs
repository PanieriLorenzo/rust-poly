use crate::num::{Complex, Float, Zero};
use na::RealField;

use crate::ScalarOps;

use super::RootFinder;

pub mod naive;
pub mod newton;

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

        // eliminate zero roots
        for _ in 0..n {
            if self.state().poly.eval_point(Complex::zero()).norm() < self.config().epsilon {
                self.state_mut().poly = self.state().poly.clone().deflate_composite(Complex::zero())
            } else {
                break;
            }
        }

        // trivial cases
        match self.state_mut().poly.degree_raw() {
            0 => return Ok(vec![]),
            1 => return Ok(self.state_mut().poly.clone().linear()),
            2 => return Ok(self.state_mut().poly.clone().quadratic()[..n].into()),
            _ => {}
        }

        for i in 0..n {
            if let Some(history_handle) = self.history() {
                history_handle.roots_history.push(vec![]);
            }
            let r = self.next_root()?;
            self.state_mut().clean_roots.extend(r.iter().copied());
            if i != (n - 1) {
                //self.state().poly = self.state().poly.clone() / Poly::from_roots(&r);
                self.state_mut().poly = self.state_mut().poly.clone().deflate_composite(r[0]);
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
