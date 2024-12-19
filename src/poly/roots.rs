use crate::{
    util::{
        complex::{c_from_f128, c_neg, c_sqrt, c_to_f128},
        vec::slice_mean,
    },
    Poly, RealScalar,
};
use anyhow::anyhow;
use f128::f128;
use itertools::Itertools;
use num::{Complex, FromPrimitive, One, Zero};

mod single_root;
pub use single_root::{halley, naive, newton};
mod all_roots;
pub use all_roots::{aberth_ehrlich, deflate, halley_deflate, naive_deflate, newton_deflate};
mod many_roots;
pub use many_roots::{halley_parallel, naive_parallel, newton_parallel, parallel};
mod initial_guess;
pub use initial_guess::{initial_guess_smallest, initial_guesses_circle};

#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error<T> {
    #[error("root finder did not converge within the given constraints")]
    NoConverge(T),

    #[error("unexpected error while running root finder")]
    Other(#[from] anyhow::Error),
}

// TODO: make a type that contains results with some extra info and an `.unpack_roots` method.

pub type Result<T> = std::result::Result<Vec<Complex<T>>, Error<Vec<Complex<T>>>>;

// TODO: use everywhere
pub type Roots<T> = Vec<Complex<T>>;

pub enum PolishingMode<T> {
    None,
    StandardPrecision {
        epsilon: T,
        min_iter: usize,
        max_iter: usize,
    },
    #[cfg(target_arch = "x86_64")]
    HighPrecision {
        epsilon: T,
        min_iter: usize,
        max_iter: usize,
    },
}

pub enum MultiplesHandlingMode<T> {
    None,
    BroadcastBest { detection_epsilon: T },
    BroadcastAverage { detection_epsilon: T },
    KeepBest { detection_epsilon: T },
    KeepAverage { detection_epsilon: T },
}

pub enum InitialGuessMode<T> {
    GuessPoolOnly,
    RandomAnnulus { bias: T, perturbation: T, seed: u64 },
    // TODO: Hull {},
    // TODO: GridSearch {},
}

impl<T: RealScalar> Poly<T> {
    /// A convenient way of finding roots, with a pre-configured root finder.
    /// Should work well for most real polynomials of low degree.
    ///
    /// Use [`Poly::roots_expert`] if you need more control over performance or accuracy.
    ///
    /// # Errors
    /// - Solver did not converge within `max_iter` iterations
    /// - Some other edge-case was encountered which could not be handled (please
    ///   report this, as we can make this solver more robust!)
    pub fn roots(&self, epsilon: T, max_iter: usize) -> Result<T> {
        self.roots_expert(
            epsilon.clone(),
            max_iter,
            0,
            PolishingMode::StandardPrecision {
                epsilon: epsilon.clone(),
                min_iter: 0,
                max_iter,
            },
            MultiplesHandlingMode::BroadcastBest {
                // TODO: tune ratio
                detection_epsilon: epsilon * T::from_f64(1.5).expect("overflow"),
            },
            &[],
            InitialGuessMode::RandomAnnulus {
                bias: T::from_f64(0.5).expect("overflow"),
                perturbation: T::from_f64(0.5).expect("overflow"),
                seed: 1,
            },
        )
    }

    /// Highly configurable root finder.
    ///
    /// [`Poly::roots`] will often be good enough, but you may know something
    /// about the polynomial you are factoring that allows you to tweak the
    /// settings.
    ///
    /// # Errors
    /// - Solver did not converge within `max_iter` iterations
    /// - Some other edge-case was encountered which could not be handled (please
    ///   report this, as we can make this solver more robust!)
    /// - The combination of parameters that was provided is invalid
    pub fn roots_expert(
        &self,
        epsilon: T,
        max_iter: usize,
        _min_iter: usize,
        polishing_mode: PolishingMode<T>,
        multiples_handling_mode: MultiplesHandlingMode<T>,
        initial_guess_pool: &[Complex<T>],
        initial_guess_mode: InitialGuessMode<T>,
    ) -> Result<T> {
        debug_assert!(self.is_normalized());

        let mut this = self.clone();

        let mut roots: Vec<Complex<T>> = this.zero_roots(epsilon.clone());

        match this.degree_raw() {
            1 => {
                roots.extend(this.linear_roots());
                return Ok(roots);
            }
            2 => {
                roots.extend(this.quadratic_roots());
                return Ok(roots);
            }
            _ => {}
        }

        this.make_monic();

        debug_assert!(this.is_normalized());
        let mut initial_guesses = Vec::with_capacity(this.degree_raw());
        for guess in initial_guess_pool.iter().cloned() {
            initial_guesses.push(guess);
        }

        // fill remaining guesses with zeros and prepare to replace with computed
        // initial guesses
        let delta = this.degree_raw() - initial_guesses.len();
        for _ in 0..delta {
            initial_guesses.push(Complex::<T>::zero());
        }
        let remaining_guesses_view =
            &mut initial_guesses[initial_guess_pool.len()..this.degree_raw()];

        match initial_guess_mode {
            InitialGuessMode::GuessPoolOnly => {
                if initial_guess_pool.len() < this.degree_raw() {
                    return Err(Error::Other(anyhow!("not enough initial guesses, you must provide one guess per root when using GuessPoolOnly")));
                }
            }
            InitialGuessMode::RandomAnnulus {
                bias,
                perturbation,
                seed,
            } => {
                initial_guesses_circle(&this, bias, seed, perturbation, remaining_guesses_view);
            } // TODO: InitialGuessMode::Hull {} => todo!(),
              // TODO: InitialGuessMode::GridSearch {} => todo!(),
        }

        log::trace!("{initial_guesses:?}");

        roots.extend(aberth_ehrlich(
            &mut this,
            Some(epsilon.clone()),
            Some(max_iter),
            &initial_guesses,
        )?);

        // further polishing of roots
        let roots: Roots<T> = match polishing_mode {
            PolishingMode::None => Ok(roots),
            PolishingMode::StandardPrecision {
                epsilon,
                min_iter,
                max_iter,
            } => newton_parallel(&mut this, Some(epsilon), Some(max_iter), &roots),

            #[cfg(target_arch = "x86_64")]
            PolishingMode::HighPrecision {
                epsilon,
                min_iter,
                max_iter,
            } => {
                let mut this = this.clone().cast_to_f128();
                let roots = roots.iter().cloned().map(|z| c_to_f128(z)).collect_vec();
                newton_parallel(
                    &mut this,
                    Some(f128::from(epsilon.to_f64().expect("overflow"))),
                    Some(max_iter),
                    &roots,
                )
                .map(|v| v.into_iter().map(|z| c_from_f128::<T>(z)).collect_vec())
                .map_err(|e| match e {
                    Error::NoConverge(v) => {
                        Error::NoConverge(v.into_iter().map(|z| c_from_f128::<T>(z)).collect_vec())
                    }
                    Error::Other(o) => Error::Other(o),
                })
            }
        }?;

        match multiples_handling_mode {
            MultiplesHandlingMode::None => Ok(roots),
            MultiplesHandlingMode::BroadcastBest { detection_epsilon } => Ok(best_multiples(
                &this,
                group_multiples(roots, detection_epsilon),
                true,
            )),
            MultiplesHandlingMode::BroadcastAverage { detection_epsilon } => Ok(average_multiples(
                &this,
                group_multiples(roots, detection_epsilon),
                true,
            )),
            MultiplesHandlingMode::KeepBest { detection_epsilon } => Ok(best_multiples(
                &this,
                group_multiples(roots, detection_epsilon),
                false,
            )),
            MultiplesHandlingMode::KeepAverage { detection_epsilon } => Ok(average_multiples(
                &this,
                group_multiples(roots, detection_epsilon),
                false,
            )),
        }
    }
}

// private
impl<T: RealScalar> Poly<T> {
    fn zero_roots(&mut self, epsilon: T) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());

        let mut roots = vec![];
        for _ in 0..self.degree_raw() {
            if self.eval(Complex::zero()).norm_sqr() < epsilon {
                roots.push(Complex::zero());
                // deflating zero roots can be accomplished simply by shifting
                *self = self.shift_down(1);
            } else {
                break;
            }
        }

        roots
    }

    fn linear_roots(&mut self) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());
        debug_assert_eq!(self.degree_raw(), 1);

        self.trim();
        if self.degree_raw() < 1 {
            return vec![];
        }

        let a = self.0[1].clone();
        let b = self.0[0].clone();

        // we found all the roots
        *self = Self::one();

        vec![-b / a]
    }

    /// Quadratic formula
    fn quadratic_roots(&mut self) -> Vec<Complex<T>> {
        debug_assert!(self.is_normalized());
        debug_assert_eq!(self.degree_raw(), 2);

        // trimming trailing almost zeros to avoid overflow
        self.trim();
        if self.degree_raw() == 1 {
            return self.linear_roots();
        }
        if self.degree_raw() == 0 {
            return vec![];
        }

        let a = self.0[2].clone();
        let b = self.0[1].clone();
        let c = self.0[0].clone();
        let four = Complex::<T>::from_u8(4).expect("overflow");
        let two = Complex::<T>::from_u8(2).expect("overflow");

        // TODO: switch to different formula when b^2 and 4c are very close due
        //       to loss of precision
        let plus_minus_term = c_sqrt(b.clone() * b.clone() - four * a.clone() * c);
        let x1 = (plus_minus_term.clone() - b.clone()) / (two.clone() * a.clone());
        let x2 = (c_neg(b.clone()) - plus_minus_term) / (two * a);

        // we found all the roots
        *self = Self::one();

        vec![x1, x2]
    }
}

/// Find roots that are within a given tolerance from each other and group them
fn group_multiples<T: RealScalar>(roots: Roots<T>, epsilon: T) -> Vec<Roots<T>> {
    // groups with their respective mean
    let mut groups: Vec<(Roots<T>, Complex<T>)> = vec![];

    let mut roots = roots;

    while !roots.is_empty() {
        // now for each root we find a group whose median is within tolerance,
        // if we don't find any we add a new group with the one root
        // if we do, we add the point, update the mean
        'roots_loop: for root in roots.drain(..) {
            for group in &mut groups {
                if (group.1.clone() - root.clone()).norm_sqr() <= epsilon {
                    group.0.push(root.clone());
                    group.1 = slice_mean(&group.0);
                    continue 'roots_loop;
                }
            }
            groups.push((vec![root.clone()], root));
        }

        // now we loop through all the groups and through each element in each group
        // and remove any elements that now lay outside of the tolerance
        for group in &mut groups {
            // hijacking retain to avoid having to write a loop where we delete
            // things from the collection we're iterating from.
            group.0.retain(|r| {
                if (r.clone() - group.1.clone()).norm_sqr() <= epsilon {
                    true
                } else {
                    roots.push(r.clone());
                    false
                }
            });
        }

        // finally we prune empty groups
        groups.retain(|g| !g.0.is_empty());
    }

    groups.into_iter().map(|(r, _)| r).collect_vec()
}

fn best_multiples<T: RealScalar>(
    poly: &Poly<T>,
    groups: Vec<Roots<T>>,
    do_broadcast: bool,
) -> Roots<T> {
    // find the best root in each group
    groups
        .into_iter()
        .flat_map(|group| {
            let len = group.len();
            let best = group
                .into_iter()
                .map(|root| (root.clone(), poly.eval(root).norm_sqr()))
                .reduce(|(a_root, a_eval), (b_root, b_eval)| {
                    if a_eval < b_eval {
                        (a_root, a_eval)
                    } else {
                        (b_root, b_eval)
                    }
                })
                .expect("empty groups not allowed")
                .0;
            if do_broadcast {
                vec![best; len]
            } else {
                vec![best]
            }
        })
        .collect_vec()
}

fn average_multiples<T: RealScalar>(
    poly: &Poly<T>,
    groups: Vec<Roots<T>>,
    do_broadcast: bool,
) -> Roots<T> {
    groups
        .into_iter()
        .flat_map(|group| {
            let len_usize = group.len();
            debug_assert!(len_usize > 0);
            let len = T::from_usize(len_usize).expect("infallible");
            let sum: Complex<T> = group.into_iter().sum();
            let avg = sum / len;
            if do_broadcast {
                vec![avg; len_usize]
            } else {
                vec![avg]
            }
        })
        .collect_vec()
}

#[cfg(test)]
mod test {
    use num::complex::ComplexFloat;

    use crate::Poly64;

    /// See [#3](https://github.com/PanieriLorenzo/rust-poly/issues/3)
    #[test]
    fn roots_of_reverse_bessel() {
        let poly = Poly64::reverse_bessel(2).unwrap();
        let roots = poly.roots(1E-10, 1000).unwrap();
        assert!((roots[0].re() - -1.5).abs() < 0.01);
        assert!((roots[0].im().abs() - 0.866).abs() < 0.01);
        assert!((roots[1].re() - -1.5).abs() < 0.01);
        assert!((roots[1].im().abs() - 0.866).abs() < 0.01);
    }
}
