use na::RealField;
use num::Float;

use crate::ScalarOps;

use super::RootFinder;

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
    fn next_n_roots(&mut self, n: usize) -> super::Result<T> {
        debug_assert!(self.state().poly.is_normalized());
        assert!(
            i32::try_from(n).expect("overflow") <= self.state().poly.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

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
}
