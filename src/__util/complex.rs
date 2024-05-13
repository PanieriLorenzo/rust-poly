// internal utilities for dealing with Complex annoyiances

use std::cmp::Ordering;

use num::{traits::float::FloatCore, Complex, Float, One, Zero};

use crate::{Scalar, ScalarOps};

// neg operator for Complex, as it does not implement std::ops::Neg
pub(crate) fn c_neg<T: Scalar>(x: Complex<T>) -> Complex<T> {
    Complex::<T>::zero() - x
}

// min based on norm1
pub(crate) fn c_min<T: Scalar + PartialOrd>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    if a.norm_sqr() < b.norm_sqr() {
        a
    } else {
        b
    }
}

pub(crate) fn c_max<T: Scalar + PartialOrd>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    if a.norm_sqr() > b.norm_sqr() {
        a
    } else {
        b
    }
}

// sort a vector of complex numbers  in place by their real component first,
// then their imaginary component
pub(crate) fn complex_sort_mut<T: Scalar + PartialOrd>(v: &mut na::DVector<Complex<T>>) {
    v.as_mut_slice().sort_by(|a, b| {
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
