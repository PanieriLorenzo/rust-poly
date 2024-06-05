use crate::{Poly, ScalarOps};

use super::RootFinder;

pub mod newton;

pub trait IterativeRootFinder<T: ScalarOps + PartialOrd>: RootFinder<T> {
    /// Find one root, without shrinkage
    fn next_root(&mut self) -> super::Result<T>;

    /// Find multiple roots by shrinking
    fn next_n_roots(&mut self, n: usize) -> super::Result<T> {
        debug_assert!(self.state().poly.is_normalized());
        assert!(
            n as i32 <= self.state().poly.degree_raw(),
            "for a polynomial of degree D, there can't be more than D roots"
        );

        for i in 0..n {
            let r = self.next_root()?;
            self.state().clean_roots.extend(r.iter().cloned());
            if i != (n - 1) {
                self.state().poly = self.state().poly.clone() / Poly::from_roots(&r);
            }
        }

        Ok(self.state().clean_roots.clone())
    }
}
