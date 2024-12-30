use std::fmt::Display;

use itertools::Itertools;
use num::{Complex, One, Zero};

use crate::{
    util::{
        complex::{c_neg, complex_fmt, complex_sort_mut},
        doc_macros::panic_absurd_size,
    },
    RealScalar,
};

mod base;
mod calculus;
mod conversions;
mod impl_num;
mod indexing;
pub mod roots;
mod special_funcs;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<T: RealScalar>(pub(crate) Vec<Complex<T>>);

impl<T: RealScalar> crate::poly2::Poly<T> for Poly<T> {
    fn degree_usize(&self) -> usize {
        debug_assert!(self.is_normalized());
        self.degree_raw()
    }
}

impl<T: RealScalar> Poly<T> {
    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p.shift_up(2), poly![0.0, 0.0, 1.0, 2.0, 3.0]);
    /// ```
    #[must_use]
    pub fn shift_up(&self, n: usize) -> Self {
        let mut v = vec![Complex::<T>::zero(); n];
        v.extend_from_slice(self.as_slice());
        Self::from_complex_vec(v)
    }

    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0, 4.0];
    /// assert_eq!(p.shift_down(2), poly![3.0, 4.0]);
    /// ```
    #[must_use]
    pub fn shift_down(&self, n: usize) -> Self {
        Self::from_complex_slice(&self.as_slice()[n..])
    }
}

impl<T: RealScalar> Poly<T> {
    pub fn new(coeffs: &[Complex<T>]) -> Self {
        Self(coeffs.to_owned()).normalize()
    }

    /// Complex constant as a polynomial of degree zero
    #[inline(always)]
    #[must_use]
    pub fn constant(c: Complex<T>) -> Self {
        Self::new(&[c])
    }

    /// Linear function as a polynomial.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{One, Zero};
    ///
    /// assert_eq!(Poly::line(Complex::one(), Complex::new(-1.0, 0.0)).eval_point(Complex::one()), Complex::zero());
    /// ```
    pub fn line(offset: Complex<T>, slope: Complex<T>) -> Self {
        if slope.is_zero() {
            return Self::new(&[offset]);
        }
        Self::new(&[offset, slope])
    }

    /// Line between two points with complex coordinates.
    ///
    /// Note that the points are determined by two complex numbers, so they are
    /// in a four dimensional space. Leave the imaginary component as zero for lines
    /// in a 2D plane.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{One, Zero};
    ///
    /// let p1 = (Complex::new(-1.0, 0.0), Complex::new(2.0, 0.0));
    /// let p2 = (Complex::new(2.0, 0.0), Complex::new(-1.0, 0.0));
    ///
    /// assert_eq!(Poly::line_from_points(p1, p2).eval_point(Complex::one()), Complex::zero());
    /// ```
    pub fn line_from_points(p1: (Complex<T>, Complex<T>), p2: (Complex<T>, Complex<T>)) -> Self {
        let slope = (p2.1 - p1.1.clone()) / (p2.0 - p1.0.clone());
        let offset = p1.1.clone() - slope.clone() * p1.0;
        Self::line(offset, slope)
    }

    /// Create a polynomial from a single term (coefficient + degree)
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    /// use num::One;
    ///
    /// assert_eq!(Poly::term(Complex::one(), 3), poly![0.0, 0.0, 0.0, 1.0]);
    /// ```
    pub fn term(coeff: Complex<T>, degree: u32) -> Self {
        Self::line(Complex::zero(), complex!(T::one())).pow(degree) * coeff
    }

    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert!(self.is_normalized());
        self.len_raw()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Raises a polynomial to an integer power.
    ///
    /// # Caveats
    /// We adopt the convention that $0^0=1$, even though some authors leave this
    /// case undefined. We believe this to be more useful as it naturally arises
    /// in integer exponentiation when defined as repeated multiplication (as we
    /// implement it).
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    ///
    /// assert_eq!(poly![1.0, 2.0, 3.0].pow(2), poly![1.0, 4.0, 10.0, 12.0, 9.0]);
    /// ```
    #[must_use]
    pub fn pow(self, pow: u32) -> Self {
        self.pow_usize(pow as usize)
    }

    /// Same as [`Poly::pow`], but takes a `usize` exponent.
    #[must_use]
    pub fn pow_usize(self, pow: usize) -> Self {
        // invariant: poly is normalized
        debug_assert!(self.is_normalized());

        if pow == 0 {
            return Self::one();
        }

        if pow == 1 {
            return self;
        }

        if self.is_zero() {
            return self;
        }

        // TODO: divide and conquer with powers of 2
        let mut res = self.clone();
        for _ in 2..=pow {
            res = res * self.clone();
        }
        res.normalize()
    }

    /// Compute the conjugate polynomial, that is a polynomial where every
    /// coefficient is conjugated.
    ///
    /// To evaluate a conjugate polynomial, you must evaluate it at the conjugate
    /// of the input, i.e. `poly.conj().eval(z.conj())`
    #[must_use]
    pub fn conj(&self) -> Self {
        Self(self.0.iter().cloned().map(|z| z.conj()).collect_vec()).normalize()
    }

    /// Get the nth term of the polynomial as a new polynomial
    ///
    /// Will return None if out of bounds.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{poly, Poly};
    /// use num::Complex;
    /// use num::One;
    ///
    /// let p  = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p.get_term(1).unwrap(), poly![0.0, 2.0]);
    /// ```
    #[must_use]
    pub fn get_term(&self, degree: u32) -> Option<Self> {
        if degree as usize >= self.len_raw() {
            return None;
        }
        Some(Self::term(self.as_slice()[degree as usize].clone(), degree))
    }

    /// Iterate over the terms of a polynomial
    ///
    /// ```
    /// # use rust_poly::{poly, Poly};
    ///
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p, p.terms().sum::<Poly<_>>());
    /// ```
    ///
    /// # Panics
    /// On polynomials with a degree higher than `u32::MAX`
    pub fn terms(&self) -> std::iter::Map<std::ops::Range<usize>, impl FnMut(usize) -> Self + '_> {
        debug_assert!(self.is_normalized());
        (0..self.len_raw()).map(|i| {
            self.get_term(
                i.try_into()
                    .expect("degrees above u32::MAX are not supported"),
            )
            .expect("terms are within range len_raw, this should never fail")
        })
    }
}

impl<T: RealScalar + PartialOrd> Poly<T> {
    /// Monic polynomial from its complex roots.
    ///
    /// # Examples
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    /// use num::{Zero, One};
    ///
    /// let p = Poly::from_roots(&[Complex::new(-1.0, 0.0), Complex::zero(), Complex::one()]);
    /// assert_eq!(p, Poly::new(&[Complex::zero(), Complex::new(-1.0, 0.0), Complex::zero(), Complex::one()]))
    /// ```
    #[must_use]
    pub fn from_roots(roots: &[Complex<T>]) -> Self {
        if roots.is_empty() {
            return Self::one();
        }

        let mut roots = roots.to_owned();
        complex_sort_mut(roots.as_mut_slice());

        roots
            .into_iter()
            .map(|e| Self::line(c_neg(e), Complex::<T>::one()))
            .fold(Self::one(), |acc, x| acc * x)
            .normalize()
    }

    /// Compose two polynomials, returning a new polynomial.
    ///
    /// Substitute the given polynomial `x` into `self` and expand the
    /// result into a new polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use rust_poly::{Poly, poly};
    /// use num::{One, Complex};
    ///
    /// let f = poly![1.0, 2.0];
    /// let g = Poly::one();
    ///
    /// assert_eq!(f.clone().compose(g), f);
    #[must_use]
    pub fn compose(self, x: Self) -> Self {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());
        debug_assert!(x.is_normalized());

        if self.is_one() {
            return x;
        }

        if x.is_one() {
            return self;
        }

        (0..self.len_raw())
            .map(|i| Self::new(&[self.0[i].clone()]) * x.clone().pow_usize(i))
            .sum()
    }
}

impl<T: RealScalar> Poly<T> {
    /// Evaluate the polynomial for each entry of a slice.
    pub fn eval_multiple(&self, points: &[Complex<T>], out: &mut [Complex<T>]) {
        debug_assert!(self.is_normalized());

        // TODO: parallelize this loop
        for (y, x) in out.iter_mut().zip(points) {
            *y = self.eval(x.clone());
        }
    }

    /// Evaluate the polynomial at a single value of `x`.
    ///
    /// ```
    /// use rust_poly::Poly;
    /// use num::Complex;
    ///
    /// let p = Poly::new(&[Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0)]);
    /// let x = Complex::new(1.0, 0.0);
    /// assert_eq!(p.eval_point(x), Complex::new(6.0, 0.0));
    /// ```
    pub fn eval(&self, x: Complex<T>) -> Complex<T> {
        // use Horner's method: https://en.wikipedia.org/wiki/Horner%27s_method
        // in theory, Estrin's method is more parallelizable, but benchmarking
        // against fast_polynomial crate shows no significant difference, this
        // is close to optimal in terms of performance. You may get some small
        // speedups by dividing large polynomials into 4 or 8 evaluations and
        // computing them in parallel using SIMD, Rayon or GPU.
        debug_assert!(self.is_normalized());
        let mut eval = self.last();

        // inlined len_raw() to ensure safety conditions are not broken by
        // updating len_raw() implementation
        let n = self.0.len();

        for i in 1..n {
            // SAFETY: index n - i - 1 is always in bounds
            // PROOF:
            // 1. the safe bounds for the index are [0, n - 1]
            // 2. i is always in the range [1, n - 1]
            // 3. the range of the index n - i - 1 is given by:
            //    n - i - 1 => n - [1, n - 1] - 1
            //             <=> [n - 1 - 1, n - 1 - (n - 1)]
            //             <=> [n - 2, 0]
            //             <=> [0, n - 2]   reversing bounds is equivalent
            // 4. the range of the index [0, n - 2] is a subset of [0, n - 1],
            //    therefore the index is in bounds. QED
            let c = (unsafe { self.0.get_unchecked(n - i - 1) }).clone();
            eval = eval * x.clone() + c;
        }
        eval
    }
}

impl<T: RealScalar + PartialOrd> Poly<T> {
    /// Translate along x-axis (or x-plane) and y-axis (or y-plane).
    ///
    /// Using complex coordinates means you'll effectively be translating in
    /// 4D space.
    #[must_use]
    pub fn translate(mut self, x: Complex<T>, y: Complex<T>) -> Self {
        self = self.compose(Self::from_complex_slice(&[c_neg(x), Complex::<T>::one()]));
        self.0[0] += y;
        self
    }
}

impl<T: RealScalar> Poly<T> {
    /// Returns true if every coefficient in the polynomial is smaller than the
    /// tolerance (using complex norm).
    ///
    /// # Examples
    /// ```
    /// use rust_poly::{Poly, poly};
    ///
    /// assert!(poly![0.01, -0.01].almost_zero(&0.1));
    /// ```
    #[must_use]
    pub fn almost_zero(&self, tolerance: &T) -> bool {
        // invariant: polynomials are normalized
        debug_assert!(self.is_normalized());

        self.as_slice().iter().all(|c| c.norm_sqr() <= *tolerance)
    }
}

impl<T: RealScalar + Display> Display for Poly<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.iter().enumerate();
        if let Some((_, c)) = iter.next() {
            write!(f, "{}", complex_fmt(c))?;
        } else {
            return Ok(());
        }
        for (i, c) in iter {
            write!(f, " + {}*x^{}", complex_fmt(c), i)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn translate() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(
            p.translate(complex!(1.0), complex!(2.0)),
            poly![4.0, -4.0, 3.0]
        );
    }

    #[test]
    fn display() {
        let p = poly![(2.0, 0.0), (4.5, 0.0), (5.0, 1.0), (6.0, 1.5), (7.0, 2.0)];
        assert_eq!(
            p.to_string(),
            "2 + 4.5*x^1 + (5+i)*x^2 + (6+i1.5)*x^3 + (7+i2)*x^4".to_string()
        );
    }

    #[test]
    fn conj() {
        let p = poly![(1.0, -1.0), (2.0, 2.0), (3.0, -3.0)];
        let q = poly![(1.0, 1.0), (2.0, -2.0), (3.0, 3.0)];
        assert_eq!(p.conj(), q);
    }
}
