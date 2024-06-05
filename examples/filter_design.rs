//! Complex polynomials are used extensively in control theory and signal processing.
//! In this example, we use Poly to find the roots of a large reverse-Bessel polynomial,
//! this is a common step in the design of Bessel filters. The `try_roots` method
//! is very stable and robust, so we can use it to design filters of very high
//! order. In this example, it's an order 48 filter.
//!
//! Then the roots are paired up in pairs, to form the denominator of the transfer
//! functions of cascaded biquads. Note that we are skipping the stage of scaling
//! the roots to obtain the required cutoff frequency. Here the cutoff is unscaled
//! and implicit.
//!
//! A full implementation would require additional steps for scaling and discretizing,
//! those are skipped here.

use rust_poly::Poly64;

fn main() {
    let num_poles = 48;
    let roots = Poly64::reverse_bessel(num_poles)
        .unwrap()
        // epsilon is the largest allowed error in the roots
        //
        // max_iter is the max number of iterations before giving up
        //
        // max_tries is how many times the solver can fail (run out of iterations)
        // before giving up entirely. When the solver gets stuck, the problem
        // is shrunk by finding a single root with a different method, this unstucks
        // the solver for some problematic polynomials. In this case its not
        // necessary so we put 1 try.
        .try_roots(1E-14, 1000, 1, None, None, None)
        .unwrap();
    roots
        .chunks(2)
        .map(Poly64::from_roots)
        .for_each(|stage| println!("{stage}"));
}
