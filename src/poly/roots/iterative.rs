use na::RealField;
use num::Float;

use crate::{Poly, ScalarOps};

use super::RootFinder;

pub mod newton;

pub enum DeflationStrategy {
    LongDivision,
    DeflateForward,
    DeflateBackward,
    DeflateComposite,
}

pub trait IterativeRootFinder<T: ScalarOps + PartialOrd + Float + RealField>:
    RootFinder<T>
{
    /// Find one root, without shrinkage
    fn next_root(&mut self) -> super::Result<T>;

    fn with_deflation_strategy(&mut self, strat: DeflationStrategy) {
        todo!()
    }

    /// Find multiple roots by shrinking
    fn next_n_roots(&mut self, n: usize) -> super::Result<T> {
        debug_assert!(self.state().poly.is_normalized());
        assert!(
            n as i32 <= self.state().poly.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

        // trivial cases
        match self.state().poly.degree_raw() {
            0 => return Ok(vec![]),
            1 => return Ok(self.state().poly.clone().linear()),
            2 => return Ok(self.state().poly.clone().quadratic()[..n].into()),
            _ => {}
        }

        for i in 0..n {
            if let Some(stats_handle) = self.statistics() {
                stats_handle.roots_history.push(vec![]);
                stats_handle.roots_err_sq_history.push(vec![]);
            }
            let r = self.next_root()?;
            self.state().clean_roots.extend(r.iter().cloned());
            if i != (n - 1) {
                //self.state().poly = self.state().poly.clone() / Poly::from_roots(&r);
                self.state().poly = self.state().poly.clone().deflate_composite(r[0]);
            }
            if let Some(stats_handle) = self.statistics() {
                stats_handle.roots_err_sq.push(
                    *stats_handle
                        .roots_err_sq_history
                        .last()
                        .expect("at least one iteration")
                        .last()
                        .expect("at least one root"),
                );
            }
        }

        if let Some(stats_handle) = self.statistics() {
            let mut accumulator = T::zero();
            for e in &stats_handle.roots_err_sq {
                accumulator += *e;
            }
            accumulator /= T::from_usize(stats_handle.roots_err_sq.len())
                .expect("cannot reliably calculate stats for massive polynomial (>> 2^16 coeffs)");
            stats_handle.rmse = Some(Float::sqrt(accumulator));
        }
        Ok(self.state().clean_roots.clone())
    }
}
