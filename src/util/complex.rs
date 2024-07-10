// internal utilities for dealing with Complex annoyiances

use std::cmp::Ordering;

use num::{Complex, One, Zero};

use crate::RealScalar;

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
