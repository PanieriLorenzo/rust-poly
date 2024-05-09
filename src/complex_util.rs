// internal utilities for dealing with Complex annoyiances

use std::cmp::Ordering;

use num::{Complex, One, Zero};

use crate::Scalar;

// neg operator for Complex, as it does not implement std::ops::Neg
pub fn c_neg<T: Scalar>(x: Complex<T>) -> Complex<T> {
    Complex::<T>::zero() - x
}

// sort a vector of complex numbers  in place by their real component first,
// then their imaginary component
pub fn complex_sort_mut<T: Scalar + PartialOrd>(v: &mut na::DVector<Complex<T>>) {
    v.as_mut_slice().sort_by(|a, b| {
        let re_ord = a.re.partial_cmp(&b.re).unwrap_or(Ordering::Equal);
        if re_ord != Ordering::Equal {
            return re_ord;
        }
        a.im.partial_cmp(&b.im).unwrap_or(Ordering::Equal)
    });
}

/// formatting for Complex, because the implementation is not good enough for me
pub fn complex_fmt<T: std::fmt::Display + Zero + One + PartialEq>(c: &Complex<T>) -> String {
    let r = &c.re;
    let i = &c.im;
    if i.is_zero() {
        format!("{}", r)
    } else if i.is_one() {
        format!("({}+i)", r)
    } else {
        format!("({}+i{})", r, i)
    }
}
