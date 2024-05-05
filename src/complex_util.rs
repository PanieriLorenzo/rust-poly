// internal utilities for dealing with Complex annoyiances

use std::cmp::Ordering;

use num::{Complex, Zero};

use crate::{FloatScalar, Scalar};

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
