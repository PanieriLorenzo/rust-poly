// internal utilities for dealing with Complex annoyiances

use std::cmp::Ordering;

use num::{Complex, FromPrimitive, One, ToPrimitive, Zero};

use itertools::Itertools;

use crate::RealScalar;

use f128::f128;

use super::vec::slice_mean;

/// cast complex to primitive
pub(crate) fn c_to_f64<T: RealScalar>(z: Complex<T>) -> Complex<f64> {
    Complex::new(
        z.re.to_f64().expect("overflow"),
        z.im.to_f64().expect("overflow"),
    )
}

/// cast complex from primitive
pub(crate) fn c_from_f64<T: RealScalar>(z: Complex<f64>) -> Complex<T> {
    Complex::new(
        T::from_f64(z.re).expect("overflow"),
        T::from_f64(z.im).expect("overflow"),
    )
}

pub(crate) fn c_to_f128<T: RealScalar>(z: Complex<T>) -> Complex<f128> {
    let z = c_to_f64(z);
    let re = f128::from(z.re);
    let im = f128::from(z.im);
    Complex::new(re, im)
}

pub(crate) fn c_from_f128<T: RealScalar>(z: Complex<f128>) -> Complex<T> {
    let re: f64 = z.re.into();
    let im: f64 = z.im.into();
    c_from_f64(Complex::new(re, im))
}

// neg operator for Complex, as it does not implement std::ops::Neg
pub(crate) fn c_neg<T: RealScalar>(x: Complex<T>) -> Complex<T> {
    Complex::<T>::zero() - x
}

// min based on norm1
pub(crate) fn c_min<T: RealScalar>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    if a.norm_sqr() < b.norm_sqr() {
        a
    } else {
        b
    }
}

/// arg using `ToPrimitive`
pub(crate) fn c_arg<T: RealScalar>(z: Complex<T>) -> T {
    let z = c_to_f64(z);
    T::from_f64(z.im.atan2(z.re)).expect("overflow")
}

/// exp using `ToPrimitive`
pub(crate) fn c_exp<T: RealScalar>(z: Complex<T>) -> Complex<T> {
    c_from_f64(c_to_f64(z).exp())
}

/// powf using `ToPrimitive`
pub(crate) fn c_powf<T: RealScalar>(z: Complex<T>, e: T) -> Complex<T> {
    c_from_f64(c_to_f64(z).powf(e.to_f64().expect("overflow")))
}

/// sqrt using `ToPrimitive`
pub(crate) fn c_sqrt<T: RealScalar>(z: Complex<T>) -> Complex<T> {
    c_from_f64(c_to_f64(z).sqrt())
}

// sort a vector of complex numbers lexicographically, using their real part first
pub(crate) fn complex_sort_mut<T: RealScalar>(v: &mut [Complex<T>]) {
    v.sort_by(|a, b| {
        let re_ord = a.re.partial_cmp(&b.re).unwrap_or(Ordering::Equal);
        if re_ord != Ordering::Equal {
            return re_ord;
        }
        a.im.partial_cmp(&b.im).unwrap_or(Ordering::Equal)
    });
}

/// formatting for Complex, because the implementation is not good enough for me
pub(crate) fn complex_fmt<T: std::fmt::Display + Zero + One + PartialEq>(c: &Complex<T>) -> String {
    let r = &c.re;
    let i = &c.im;
    if i.is_zero() {
        format!("{r}")
    } else if i.is_one() {
        format!("({r}+i)")
    } else {
        format!("({r}+i{i})")
    }
}

pub fn complex_mean<T: RealScalar>(v: &[Complex<T>]) -> Complex<T> {
    let (res, ims): (Vec<_>, Vec<_>) = v.iter().map(|z| (z.re.clone(), z.im.clone())).unzip();
    Complex::new(slice_mean(&res), slice_mean(&ims))
}

#[cfg(test)]
mod test {
    use num::Complex;

    #[test]
    fn c_to_f64() {
        let src: Complex<f32> = Complex::new(-12.34, 56.78);
        let dst = super::c_to_f64(src);
        assert_eq!(f64::from(src.re), dst.re);
        assert_eq!(f64::from(src.im), dst.im);
    }

    #[test]
    fn c_from_f64() {
        let src: Complex<f64> = Complex::new(-12.34, 56.78);
        let dst: Complex<f32> = super::c_from_f64(src);
        assert_eq!(src.re as f32, dst.re);
        assert_eq!(src.im as f32, dst.im);
    }

    #[test]
    fn c_arg() {
        let src = Complex::new(0.25, 0.75);
        assert_eq!(super::c_arg(src), src.arg());

        let src = Complex::new(-0.25, 0.75);
        assert_eq!(super::c_arg(src), src.arg());

        let src = Complex::new(0.25, -0.75);
        assert_eq!(super::c_arg(src), src.arg());

        let src = Complex::new(-0.25, -0.75);
        assert_eq!(super::c_arg(src), src.arg());
    }

    #[test]
    fn c_exp() {
        let src = Complex::new(0.25, -0.75);
        assert_eq!(super::c_exp(src), src.exp());
    }

    #[test]
    fn c_powf() {
        let src = Complex::new(0.25, -0.75);
        assert_eq!(super::c_powf(src, 1.23), src.powf(1.23));
    }

    #[test]
    fn c_sqrt() {
        let src = Complex::new(0.25, -0.75);
        assert_eq!(super::c_sqrt(src), src.sqrt());
    }
}
