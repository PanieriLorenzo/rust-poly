//! The new API, will replace the [`poly`] module.

use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use num::{complex::ComplexFloat, FromPrimitive, One};

use crate::{
    num::Zero,
    scalar::ComplexScalar,
    storage::{BaseStore, MutStore, OwnedStore, OwnedUniStore, UniStore},
    util::doc_macros::panic_absurd_size,
};

use itertools::Itertools;

pub mod aliases;
pub mod poly_base;
pub mod roots;
pub mod roots_base;

/// A generic polynomial
///
/// This trait makes very few assumptions. May be univariate or multivariate,
/// integer, real or complex-valued, owned or borrowed, etc...
///
/// Must have a valid way to represent a zero polynomial, this is expressed by
/// implementing [`crate::num::Zero`]. A zero polynomial behaves like the neutral
/// element of addition over polynomials. A zero polynomial has degree -1 by
/// convention.
pub trait Poly<T>:
    ToOwned
    + Add<Self, Output = Self::Owned>
    + Mul<Self, Output = Self::Owned>
    + Sub<Self, Output = Self::Owned>
    + Div<Self, Output = Self::Owned>
    + Rem<Self, Output = Self::Owned>
where
    Self::Owned: OwnedPoly<T>,
    Self: Sized,
    // <<Self as poly2::Poly<T>>::BackingStorage as ToOwned>::Owned: OwnedStore<T>,
    <Self::BackingStorage as ToOwned>::Owned: OwnedStore<T>,
    <Self::Owned as Poly<T>>::BackingStorage: OwnedStore<T>,
{
    type BackingStorage: BaseStore<T>;

    /// Implementation detail
    #[doc(hidden)]
    fn __as_store(&self) -> &Self::BackingStorage;

    /// Implementation detail
    #[doc(hidden)]
    fn __from_store(store: Self::BackingStorage) -> Self;

    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    /// Return the degree as `usize`.
    ///
    /// Note that unlike [`Poly::degree`], this will saturate at 0 for zero
    /// polynomials. As the degree of zero polynomials is undefined.
    fn degree_usize(&self) -> usize;

    /// The total number of coefficients in the polynomial.
    fn size(&self) -> usize;

    /// The degree of a polynomial (the maximum exponent)
    ///
    /// Note that this will return `-1` for zero polynomials. The degree of
    /// zero polynomials is undefined, but we use the `-1` convention adopted
    /// by some authors \[which authors?\].
    ///
    /// # Panics
    #[doc = panic_absurd_size!()]
    fn degree(&self) -> i64 {
        if self.is_zero() {
            return -1;
        }
        self.degree_usize()
            .try_into()
            .expect("usize did not fit into i64")
    }

    /// Raises a polynomial to an integer power.
    ///
    /// # Caveats
    /// We adopt the convention that $0^0=1$, even though some authors leave this
    /// case undefined. We believe this to be more useful as it naturally arises
    /// in integer exponentiation when defined as repeated multiplication.
    #[must_use]
    fn pow(&self, pow: u32) -> Self::Owned {
        self.pow_usize(pow as usize)
    }

    /// Same as [`Poly::pow`], but takes a `usize` exponent.
    fn pow_usize(&self, pow: usize) -> Self::Owned;

    fn terms(&self) -> impl Iterator<Item = Self::Owned>;

    /// Return an iterator over the coefficients in ascending order of degree
    #[inline]
    #[must_use]
    fn coeffs<'a>(&'a self) -> impl Iterator<Item = &'a T>
    where
        T: 'a,
    {
        self.__as_store().iter()
    }

    /// Polynomial composition.
    fn compose(&self, other: Self::Owned) -> Self::Owned;

    fn eval(&self, x: T) -> T;

    /// Evaluate the polynomial for each entry of a slice. May be faster than
    /// repeatedly calling [`Poly::eval`].
    fn eval_multiple(&self, points: &[T], out: &mut [T])
    where
        T: Clone,
    {
        // TODO: parallelize this loop
        for (y, x) in out.iter_mut().zip(points) {
            *y = self.eval(x.clone());
        }
    }

    /// Compute the conjugate polynomial, that is a polynomial where every
    /// coefficient is conjugated.
    ///
    /// To evaluate a conjugate polynomial, you must evaluate it at the conjugate
    /// of the input, i.e. `poly.conj().eval(z.conj())`
    #[must_use]
    fn conj(&self) -> Self::Owned
    where
        T: Clone + ComplexFloat,
    {
        let data = self.__as_store();
        Self::Owned::__from_store(
            <Self::Owned as Poly<T>>::BackingStorage::from_iter(
                &data.shape(),
                data.iter().copied().map(num::complex::ComplexFloat::conj),
            )
            .unwrap(),
        )
    }

    #[inline]
    fn as_slice(&self) -> &[T] {
        self.__as_store().as_slice()
    }

    #[must_use]
    fn to_vec(&self) -> Vec<T>
    where
        T: Clone,
    {
        Vec::from(self.as_slice())
    }
}

/// Univariate polynomial
pub trait UniPoly<T>: Poly<T>
where
    Self::Owned: OwnedUniPoly<T>,
    Self::BackingStorage: UniStore<T>,
    // HACK: can't this diamond be flattened?
    <Self::BackingStorage as ToOwned>::Owned: OwnedUniStore<T>,
    <Self::Owned as Poly<T>>::BackingStorage: OwnedUniStore<T>,
{
    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0];
    /// assert_eq!(p.shift_up(2), poly![0.0, 0.0, 1.0, 2.0, 3.0]);
    /// ```
    #[must_use]
    fn shift_up(&self, n: usize) -> Self::Owned
    where
        T: Zero + Clone,
    {
        let mut v = <Self::Owned as Poly<T>>::BackingStorage::zeros(&[n]);
        v.extend(self.__as_store().iter().cloned());
        Self::Owned::__from_store(v)
    }

    /// # Examples
    /// ```
    /// # use rust_poly::{poly, Poly};
    /// let p = poly![1.0, 2.0, 3.0, 4.0];
    /// assert_eq!(p.shift_down(2), poly![3.0, 4.0]);
    /// ```
    #[must_use]
    fn shift_down(&self, n: usize) -> Self::Owned
    where
        T: Clone,
    {
        let s = self.__as_store();
        let v = <Self::Owned as Poly<T>>::BackingStorage::from_iter(
            &s.shape(),
            s.iter().skip(n).cloned(),
        )
        .unwrap();
        Self::Owned::__from_store(v)
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
    fn get_term(&self, degree: i64) -> Option<Self::Owned>
    where
        T: Clone + num::Zero,
    {
        Some(Self::Owned::term(
            self.__as_store().iter().nth(degree as usize)?.clone(),
            degree,
        ))
    }

    /// Translate along x-axis and y-axis.
    #[must_use]
    fn translate(&self, x: T, y: T) -> Self::Owned
    where
        T: Neg<Output = T> + One,
    {
        let s = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[2], [-x, T::one()]).unwrap();
        let x_trans = self.compose(Self::Owned::__from_store(s));
        x_trans + Self::Owned::constant(y)
    }

    /// Derivative
    ///
    /// # Panics
    /// On very large degree polynomials coefficients may overflow in `T`
    #[must_use]
    fn diff(self) -> Self::Owned
    where
        T: FromPrimitive + Clone + Mul<T, Output = T>,
    {
        // derivative of constant is zero
        if self.degree() == 0 {
            return Self::Owned::zero();
        }

        let coeffs = (0..self.size())
            // HACK: this douple try_from is because Complex<f32> or Complex<f64>
            //       does not implement TryFrom<usize>, but do implement TryFrom<u64>
            .map(|x| T::from_usize(x).expect("overflow"))
            // .map(|x| Complex::new(x, T::zero()))
            .zip(self.coeffs())
            .map(|(n, c)| n * c.clone())
            .skip(1); // shift degrees down
        let store = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[self.size() - 1], coeffs)
            .unwrap();
        Self::Owned::__from_store(store)
    }

    /// Antiderivative (with C=0)
    ///
    /// # Panics
    /// On very large degree polynomials coefficients may overflow in `T`
    #[must_use]
    fn integral(self) -> Self::Owned
    where
        T: Zero + FromPrimitive + Clone + Div<T, Output = T>,
    {
        let coeffs = itertools::chain(
            [T::zero()],
            (1..=self.size())
                .map(|x| T::from_usize(x).expect("overflow"))
                .zip(self.coeffs())
                .map(|(n, c)| c.clone() / n),
        );
        let store = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[self.size() + 1], coeffs)
            .unwrap();
        Self::Owned::__from_store(store)
    }
}

pub trait MutPoly<T>: Poly<T>
where
    Self::Owned: OwnedPoly<T>,
    Self::BackingStorage: MutStore<T>,
    <Self::BackingStorage as ToOwned>::Owned: OwnedStore<T>,
    <Self::Owned as Poly<T>>::BackingStorage: OwnedStore<T>,
{
    /// Implementation detail
    #[doc(hidden)]
    fn __as_store_mut(&mut self) -> &mut Self::BackingStorage;

    #[inline]
    fn as_mut_slice(&mut self) -> &mut [T] {
        self.__as_store_mut().as_mut_slice()
    }

    #[inline]
    #[must_use]
    fn coeffs_mut<'a>(&'a mut self) -> impl Iterator<Item = &'a mut T>
    where
        T: 'a,
    {
        self.as_mut_slice().iter_mut()
    }
}

pub trait OwnedPoly<T>: Poly<T, Owned = Self>
where
    Self::BackingStorage: OwnedStore<T>,
{
    fn zero() -> Self;

    fn one() -> Self;

    fn constant(val: T) -> Self {
        let s = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[1], [val]).unwrap();
        Self::__from_store(s)
    }
}

pub trait OwnedUniPoly<T>: OwnedPoly<T> + UniPoly<T>
where
    <Self::BackingStorage as ToOwned>::Owned: OwnedUniStore<T>,
    <Self::Owned as Poly<T>>::BackingStorage: OwnedUniStore<T>,
{
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
    fn line_from_points(p1: (T, T), p2: (T, T)) -> Self
    where
        T: Clone + Sub<Output = T> + Div<Output = T> + Mul<Output = T>,
    {
        let slope = (p2.1 - p1.1.clone()) / (p2.0 - p1.0.clone());
        let offset = p1.1.clone() - slope.clone() * p1.0;
        let v = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[2], [offset, slope]).unwrap();
        Self::__from_store(v)
    }

    fn term(coeff: T, degree: i64) -> Self
    where
        T: Clone + num::Zero,
    {
        let v = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[1], [coeff]).unwrap();
        Self::__from_store(v).shift_up(degree as usize)
    }

    fn from_iter(iter: impl IntoIterator<Item = T>) -> Self {
        let data = iter.into_iter().collect_vec();
        let v = <Self::Owned as Poly<T>>::BackingStorage::from_iter(&[data.len()], data).unwrap();
        Self::__from_store(v)
    }

    #[must_use]
    fn from_real_iterator(coeffs: impl IntoIterator<Item = T::ComponentScalar>) -> Self
    where
        T: ComplexScalar,
    {
        Self::from_iter(coeffs.into_iter().map(|x| T::from_real(x)))
    }

    /// Fit a polynomial to a set of points or constraints on the derivatives
    /// of the polynomial at these points.
    ///
    /// Points or constraints are given as a list of tuples, the first value
    /// is the x coordinate, the second is the y coordinate, the final value
    /// is the order of the derivative. For simple polynomial interpolation,
    /// without any constraints on the derivatives, the last value should be 0.
    fn fit(constraints: &[(T, T, usize)]) -> Self {
        todo!()
    }

    /// Get the nth [Chebyshev polynomial](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
    /// of the first kind.
    ///
    /// ```
    /// use rust_poly::{poly, Poly};
    ///
    /// assert_eq!(Poly::cheby(2), poly![-1.0, 0.0, 2.0]);
    /// assert_eq!(Poly::cheby(3), poly![0.0, -3.0, 0.0, 4.0]);
    /// assert_eq!(Poly::cheby(4), poly![1.0, 0.0, -8.0, 0.0, 8.0])
    /// ```
    // TODO: technically it can panic in some extreme cases, would need to
    //       do some boundary testing to write a proper doc comment
    #[allow(clippy::missing_panics_doc)]
    #[must_use]
    fn cheby1(n: usize) -> Self
    where
        T: Zero + One + FromPrimitive + Sub<T, Output = T>,
    {
        // TODO: make the first 32-ish explicit for performance
        match n {
            0 => Self::from_iter([T::one()]),
            1 => Self::from_iter([T::zero(), T::one()]),
            2 => Self::from_iter([
                T::zero() - T::one(),
                T::zero(),
                T::from_u8(2).expect("overflow"),
            ]),
            3 => Self::from_iter([
                T::zero(),
                T::zero() - T::from_u8(3).expect("overflow"),
                T::zero(),
                T::from_u8(4).expect("overflow"),
            ]),
            4 => Self::from_iter([
                T::one(),
                T::zero(),
                T::zero() - T::from_u8(8).expect("overflow"),
                T::zero(),
                T::from_u8(8).expect("overflow"),
            ]),
            _ => {
                Self::from_iter([T::zero(), T::from_u8(2).expect("overflow")]) * Self::cheby1(n - 1)
                    - Self::cheby1(n - 2)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::poly2::{Poly, UniPoly};

    #[test]
    fn diff() {
        let p = poly![1.0f64, 2.0, 3.0];
        assert_eq!(p.diff(), poly![2.0, 6.0]);
    }

    /// This was a bug
    #[test]
    fn diff1() {
        let one = poly![1.0];
        assert_eq!(one.diff().degree(), -1);
    }

    #[test]
    fn integral() {
        let p = poly![1.0, 2.0, 3.0];
        assert_eq!(p.integral(), poly![0.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn integral_diff() {
        let p = poly![1.0, 2.0, 3.0];
        let q = p.clone().integral().diff();
        assert_eq!(p, q);
    }
}
