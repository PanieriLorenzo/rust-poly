//! This example shows how the accuracy and performance of root finders can be tuned

use rust_poly::{
    poly,
    roots::{Newton, RootFinder},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // some arbitrary polynomial
    let poly = poly![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

    // set up the root finder. Note the `collect_history` argument. This enables
    // collection of all intermediate guesses, which normally are disabled as
    // they slow things down
    let mut solver = Newton::from_poly(poly)
        .with_epsilon(1E-8)
        .with_max_iter(1000)
        .collect_history();

    // we don't actually care about the roots as we're just benchmarking
    let _ = solver.roots()?;

    // now we can inspect the history
    println!(
        "total iterations: {:?}",
        solver.history().as_ref().unwrap().total_iter()
    );

    Ok(())
}
