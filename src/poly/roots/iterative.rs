use crate::{
    num::{Complex, Float},
    Poly, Scalar,
};
use na::RealField;

use crate::ScalarOps;

mod naive;
pub use naive::naive;
mod newton;
pub use newton::newton;
mod halley;
pub use halley::halley;

pub type NextRootFun<T> =
    fn(
        poly: &Poly<T>,
        epsilon: T,
        max_iter: Option<usize>,
        initial_guess: Option<Complex<T>>,
    ) -> std::result::Result<(Vec<Complex<T>>, u128), super::Error<Vec<Complex<T>>>>;

pub(crate) fn deflate<T: ScalarOps + RealField>(
    next_root_fun: NextRootFun<T>,
    poly: &mut Poly<T>,
    epsilon: Option<T>,
    max_iter: Option<usize>,
    initial_guesses: &[Complex<T>],
) -> super::Result<T> {
    let mut eval_counter = 0;
    let epsilon = epsilon.unwrap_or(T::tiny_safe());
    let mut roots = vec![];
    let mut initial_guesses = initial_guesses.iter().copied();

    // until we've found all roots
    loop {
        let trivial_roots = poly.trivial_roots(epsilon);
        eval_counter += trivial_roots.1;
        roots.extend(trivial_roots.0.iter());

        debug_assert!(poly.is_normalized());
        if poly.degree_raw() == 0 {
            log::debug!("{{evaluations: {eval_counter}}}");
            return Ok(roots);
        }

        let next_guess = initial_guesses.next();
        let (next_roots, num_evals) = next_root_fun(poly, epsilon, max_iter, next_guess)?;
        let root = next_roots[0];
        eval_counter += num_evals;
        roots.push(root);
        // TODO: deflate_composite should borrow instead
        *poly = poly.clone().deflate_composite(root);
    }
}

pub enum DeflationStrategy {
    LongDivision,
    DeflateForward,
    DeflateBackward,
    DeflateComposite,
}

/// Garwick & Ward stopping criterion (see [Nikolajsen 2014](https://doi.org/10.1098/rsos.140206))
// TODO: with specialization use Nikolajsen 2014 if T is f64 or f32, but right
//       now this is fine for all real-like, including fractions and infinite
//       precision
fn stopping_criterion_garwick<T: ScalarOps>(
    z: Complex<T>,
    z_old: Complex<T>,
    z_old_old: Complex<T>,
) -> bool {
    let delta_z = (z - z_old).norm();
    let delta_z_old = (z_old - z_old_old).norm();
    let z_norm = z_old.norm();
    let em3 = T::from_f64(1E-3).expect("overflow");
    let em4 = T::from_f64(1E-4).expect("overflow");
    let em7 = T::from_f64(1E-7).expect("overflow");
    (z_norm < em4 && delta_z <= em7 || z_norm >= em4 && delta_z / z_norm <= em3)
        && delta_z >= delta_z_old
}

/// This struct lazily computes derivatives upon request, so that they are only
/// used by methods that require them.
///
/// Many methods use derivatives, but they don't always need to compute them.
/// For ease of implementation, all methods that require derivatives can use
/// this one type.
pub(crate) struct LazyDerivatives<'a, T: Scalar> {
    zeroth: &'a Poly<T>,
    first_and_higher: Vec<Poly<T>>,
}

impl<'a, T: Scalar> LazyDerivatives<'a, T> {
    pub fn new(poly: &'a Poly<T>) -> Self {
        Self {
            zeroth: poly,
            first_and_higher: vec![],
        }
    }

    pub fn get_nth_derivative(&mut self, n: usize) -> &Poly<T> {
        if n == 0 {
            return self.zeroth;
        }

        if self.first_and_higher.is_empty() {
            self.first_and_higher.push(self.zeroth.clone().diff());
        }

        for _ in self.first_and_higher.len()..n {
            let next_diff = self
                .first_and_higher
                .last()
                .expect("infallible")
                .clone()
                .diff();
            self.first_and_higher.push(next_diff);
        }
        &self.first_and_higher[n - 1]
    }
}

/// Estimate root multiplicity using Lagouanelle 1966
fn multiplicity_lagouanelle<T: Scalar>(
    px: Complex<T>,
    pdx: Complex<T>,
    pddx: Complex<T>,
) -> Complex<T> {
    let pdx_2 = pdx * pdx;
    pdx_2 / (pdx_2 - px * pddx)
}

/// Speed up convergence using Madsen 1973
fn line_search_accelerate<T: ScalarOps>(
    poly: &Poly<T>,
    guess: Complex<T>,
    delta: Complex<T>,
) -> (Complex<T>, u128) {
    let mut eval_count = 0;

    let mut guess_best = guess - delta;
    let mut px_best_norm = poly.eval(guess_best).norm();
    eval_count += 1;
    for p in 2..=poly.degree_raw() {
        let step_size = T::from_usize(p).expect("overflow");
        let guess_new = guess - delta.scale(step_size);
        let px_new = poly.eval(guess_new);
        eval_count += 1;
        let px_new_norm = px_new.norm();

        if px_new_norm >= px_best_norm {
            // stop searching as soon as it stops improving
            break;
        }
        log::trace!(
            "accelerated {{step_size: \"{step_size:?}\", before: \"{guess_best:?}\", after: \"{guess_new:?}\"}}"
        );
        px_best_norm = px_new_norm;
        guess_best = guess_new;
    }

    (guess_best, eval_count)
}

/// Slow down convergence using Madsen 1973
fn line_search_decelerate<T: ScalarOps>(
    poly: &Poly<T>,
    guess: Complex<T>,
    delta: Complex<T>,
) -> (Complex<T>, u128) {
    // arbitrary constants
    const MAX_STEPS: u32 = 2;
    const ROTATION_RADIANS: f64 = 0.643_501_108_793_284_4 /* atan(0.75) */;
    // TODO: when const trait methods are supported, this should be
    //       made fully const.
    let rotation = Complex::from_polar(T::one(), T::from_f64(ROTATION_RADIANS).expect("overflow"));

    let mut eval_count = 0;

    let mut delta_best = delta;
    let mut px_best_norm = poly.eval(guess - delta_best).norm();
    eval_count += 1;
    for p in 1..=MAX_STEPS {
        let step_size = T::from_i32(2i32.pow(p)).expect("overflow").recip();
        let delta_new = delta.scale(step_size);
        let guess_new = guess - delta_new;
        let px_new = poly.eval(guess_new);
        eval_count += 1;
        let px_new_norm = px_new.norm();
        if px_new_norm >= px_best_norm {
            // stop searching as soon as it stops improving
            return (guess_new, eval_count);
        }
        log::trace!("decelerated {{step_size: \"{step_size:?}\", after: \"{guess_new:?}\"}}");
        px_best_norm = px_new_norm;
        delta_best = delta_new;
    }

    (guess - delta_best * rotation, eval_count)
}

#[cfg(test)]
mod test {

    use na::Complex;
    use num::complex::ComplexFloat;

    use crate::Poly64;

    use super::LazyDerivatives;

    #[test]
    fn initial_guess_smallest() {
        assert!(
            (poly![24.0, -14.0, -13.0, 2.0, 1.0].initial_guess_smallest()
                - Complex::new(0.68, 0.0))
            .norm()
                < 0.01
        );
    }

    #[test]
    fn lazy_derivative() {
        let poly = poly![1.0, 2.0, 3.0, 4.0];
        let mut lazy = LazyDerivatives::new(&poly);
        assert_eq!(*lazy.get_nth_derivative(0), poly![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(lazy.first_and_higher.len(), 0);
        assert_eq!(*lazy.get_nth_derivative(1), poly![2.0, 6.0, 12.0]);
        assert_eq!(lazy.first_and_higher.len(), 1);
        assert_eq!(*lazy.get_nth_derivative(2), poly![6.0, 24.0]);
        assert_eq!(lazy.first_and_higher.len(), 2);
        assert_eq!(*lazy.get_nth_derivative(3), poly![24.0]);
        assert_eq!(lazy.first_and_higher.len(), 3);
        assert_eq!(*lazy.get_nth_derivative(4), poly![0.0]);
        assert_eq!(lazy.first_and_higher.len(), 4);
    }

    #[test]
    fn lazy_derivative_out_of_order() {
        let poly = poly![1.0, 2.0, 3.0, 4.0];
        let mut lazy = LazyDerivatives::new(&poly);
        assert_eq!(*lazy.get_nth_derivative(2), poly![6.0, 24.0]);
        assert_eq!(lazy.first_and_higher.len(), 2);
        assert_eq!(*lazy.get_nth_derivative(0), poly![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(lazy.first_and_higher.len(), 2);
        assert_eq!(*lazy.get_nth_derivative(3), poly![24.0]);
        assert_eq!(lazy.first_and_higher.len(), 3);
        assert_eq!(*lazy.get_nth_derivative(4), poly![0.0]);
        assert_eq!(lazy.first_and_higher.len(), 4);
        assert_eq!(*lazy.get_nth_derivative(1), poly![2.0, 6.0, 12.0]);
        assert_eq!(lazy.first_and_higher.len(), 4);
    }
}
