//! This example shows how the accuracy and performance of root finders can be tuned

use rust_poly::{
    poly,
    roots::{NewtonFinder, RootFinder},
};

fn main() {
    // some arbitrary polynomial
    let poly = poly![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

    // set up the root finder. Note the `collect_stats` argument. This enables
    // collection of stats, which normally are disabled as they slow things down
    let mut solver = NewtonFinder::from_poly(poly)
        .with_epsilon(1E-8)
        .with_max_iter(1000)
        .collect_stats();

    // we don't actually care about the roots as we're just benchmarking
    let _ = solver.roots().unwrap();

    // now we can inspect the statistics
    println!("{:?}", solver.statistics().as_ref().unwrap().roots_err_sq);
}
