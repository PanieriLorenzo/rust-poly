//! This example showcases how the different initial guesses are distributed,
//! by plotting them using the Plotly crate.

use std::iter;

use itertools::Itertools;
use num::{complex::Complex64, Complex};
use plotly::{
    color::NamedColor,
    common::{Marker, MarkerSymbol, Mode},
    Plot, Scatter,
};
use rust_poly::{
    complex,
    roots::{initial_guesses_circle, initial_guesses_random},
    Poly64,
};

fn main() {
    const DEGREE: usize = 100;
    let (p, roots_re, roots_im) = random_poly(DEGREE);
    example_circle::<DEGREE>(
        &p,
        &roots_re,
        &roots_im,
        1.0,
        0.0,
        "temp/guesses_circle_max_bias.svg",
    );
    example_circle::<DEGREE>(
        &p,
        &roots_re,
        &roots_im,
        0.0,
        0.0,
        "temp/guesses_circle_min_bias.svg",
    );
    example_circle::<DEGREE>(
        &p,
        &roots_re,
        &roots_im,
        0.5,
        0.0,
        "temp/guesses_circle_mid_bias.svg",
    );
    example_circle::<DEGREE>(
        &p,
        &roots_re,
        &roots_im,
        0.5,
        1.0,
        "temp/guesses_circle_max_perturb.svg",
    );
    example_circle::<DEGREE>(
        &p,
        &roots_re,
        &roots_im,
        0.03,
        0.01,
        "temp/guesses_circle_realistic.svg",
    );
    example_random::<DEGREE>(&p, &roots_re, &roots_im, "temp/guesses_random.svg");
}

/// Build and plot an example of the initial_guesses_circle function.
fn example_circle<const DEGREE: usize>(
    p: &Poly64,
    roots_re: &[f64],
    roots_im: &[f64],
    bias: f64,
    perturbation: f64,
    location: &str,
) {
    let mut guesses = [complex!(0.0); DEGREE];
    initial_guesses_circle(p, bias, 1, perturbation, &mut guesses);
    plot(roots_re, roots_im, &guesses, location);
}

/// Build and plot an example of the initial_guesses_random function.
fn example_random<const DEGREE: usize>(
    p: &Poly64,
    roots_re: &[f64],
    roots_im: &[f64],
    location: &str,
) {
    let mut guesses = [complex!(0.0); DEGREE];
    initial_guesses_random(p.clone(), 1, &mut guesses);
    plot(roots_re, roots_im, &guesses, location);
}

/// Crate a random polynomial of degree D with roots contained in the unit circle,
/// returning the polynomial, as well as the real and imaginary components of
/// the roots.
fn random_poly(degree: usize) -> (Poly64, Vec<f64>, Vec<f64>) {
    let mut re_rng = fastrand::Rng::with_seed(1);
    let mut im_rng = re_rng.fork();
    let (roots_re, roots_im, roots): (Vec<_>, Vec<_>, Vec<_>) =
        iter::from_fn(move || Some((re_rng.f64(), im_rng.f64())))
            .take(degree)
            .map(|(rho, theta)| {
                let root = Complex::from_polar(rho, theta * std::f64::consts::TAU);
                (root.re, root.im, root)
            })
            .multiunzip();
    (Poly64::from_roots(&roots), roots_re, roots_im)
}

/// Plot roots and guesses using Plotly
fn plot(roots_re: &[f64], roots_im: &[f64], guesses: &[Complex64], location: &str) {
    let mut plot = Plot::new();
    let trace = Scatter::new(roots_re.to_owned(), roots_im.to_owned())
        .mode(Mode::Markers)
        .marker(
            Marker::new()
                .symbol(MarkerSymbol::CircleXOpen)
                .size(15)
                .color(NamedColor::Blue),
        )
        .name("actual roots");
    plot.add_trace(trace);

    let (guesses_re, guesses_im) = guesses
        .iter()
        .map(|z| (z.re, z.im))
        .unzip::<_, _, Vec<_>, Vec<_>>();
    let trace = Scatter::new(guesses_re, guesses_im)
        .mode(Mode::Markers)
        .marker(
            Marker::new()
                .symbol(MarkerSymbol::SquareXOpen)
                .size(12)
                .color(NamedColor::Red),
        )
        .name("guesses");
    plot.add_trace(trace);

    plot.write_image(location, plotly::ImageFormat::SVG, 800, 600, 1.0);
}
