use crate::{num::Complex, Poly, Poly2, RealScalar};

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

/// Garwick & Ward stopping criterion (see [Nikolajsen 2014](https://doi.org/10.1098/rsos.140206))
// TODO: with specialization use Nikolajsen 2014 if T is f64 or f32, but right
//       now this is fine for all real-like, including fractions and infinite
//       precision
fn stopping_criterion_garwick<T: RealScalar>(
    z: Complex<T>,
    z_old: Complex<T>,
    z_old_old: Complex<T>,
) -> bool {
    let delta_z = (z - z_old.clone()).norm_sqr();
    let delta_z_old = (z_old.clone() - z_old_old).norm_sqr();
    let z_norm = z_old.norm_sqr();
    let em3 = T::from_f64(1E-3).expect("overflow");
    let em4 = T::from_f64(1E-4).expect("overflow");
    let em7 = T::from_f64(1E-7).expect("overflow");
    (z_norm < em4 && delta_z <= em7 || z_norm >= em4 && delta_z.clone() / z_norm <= em3)
        && delta_z >= delta_z_old
}

/// This struct lazily computes derivatives upon request, so that they are only
/// used by methods that require them.
///
/// Many methods use derivatives, but they don't always need to compute them.
/// For ease of implementation, all methods that require derivatives can use
/// this one type.
pub struct LazyDerivatives<'a, T: RealScalar> {
    zeroth: &'a Poly<T>,
    first_and_higher: Vec<Poly<T>>,
}

impl<'a, T: RealScalar> LazyDerivatives<'a, T> {
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
#[allow(clippy::similar_names)]
fn multiplicity_lagouanelle<T: RealScalar>(
    px: Complex<T>,
    pdx: Complex<T>,
    pddx: Complex<T>,
) -> Complex<T> {
    let pdx_2 = pdx.clone() * pdx;
    pdx_2.clone() / (pdx_2 - px * pddx)
}

/// Speed up convergence using Madsen 1973
fn line_search_accelerate<T: RealScalar>(
    poly: &Poly<T>,
    guess: Complex<T>,
    delta: Complex<T>,
) -> (Complex<T>, u128) {
    let mut eval_count = 0;

    let mut guess_best = guess.clone() - delta.clone();
    let mut px_best_norm = poly.eval(guess_best.clone()).norm_sqr();
    eval_count += 1;
    for p in 2..=poly.degree_raw() {
        let step_size = T::from_usize(p).expect("overflow");
        let guess_new = guess.clone() - delta.clone().scale(step_size.clone());
        let px_new = poly.eval(guess_new.clone());
        eval_count += 1;
        let px_new_norm = px_new.norm_sqr();

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
fn line_search_decelerate<T: RealScalar>(
    poly: &Poly<T>,
    guess: Complex<T>,
    delta: Complex<T>,
) -> (Complex<T>, u128) {
    // arbitrary constants
    const MAX_STEPS: u32 = 2;
    // TODO: when const trait methods are supported, this should be
    //       made fully const.
    // a rotation of about 53 degrees
    let rotation = Complex::new(
        T::from_f64(0.6).expect("overflow"),
        T::from_f64(0.8).expect("overflow"),
    );

    let mut eval_count = 0;

    let mut delta_best = delta.clone();
    let mut px_best_norm = poly.eval(guess.clone() - delta_best.clone()).norm_sqr();
    eval_count += 1;
    for p in 1..=MAX_STEPS {
        let step_size = T::from_i32(2i32.pow(p)).expect("overflow").recip();
        let delta_new = delta.scale(step_size.clone());
        let guess_new = guess.clone() - delta_new.clone();
        let px_new = poly.eval(guess_new.clone());
        eval_count += 1;
        let px_new_norm = px_new.norm_sqr();
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

    use num::Complex;

    use crate::roots::initial_guess::initial_guess_smallest;

    use super::LazyDerivatives;

    #[test]
    fn test_initial_guess_smallest() {
        assert!(
            (initial_guess_smallest(&poly![24.0, -14.0, -13.0, 2.0, 1.0])
                - Complex::new(0.68, 0.0))
            .norm()
                < 0.2
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
