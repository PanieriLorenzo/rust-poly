use std::ops::{Add, Div};

use num::{FromPrimitive, Zero};

pub fn slice_mean<T: Zero + Add<Output = T> + Div<Output = T> + FromPrimitive + Clone>(
    v: &[T],
) -> T {
    let num = v.len();
    v.iter().fold(T::zero(), |a, b| a + b.clone()) / T::from_usize(num).expect("overflow")
}

pub fn slice_is_sorted<T, K: PartialOrd>(v: &[T], key: impl Fn(&T) -> K) -> bool {
    if v.len() < 2 {
        return true;
    }
    v.iter()
        .zip(v.iter().skip(1))
        .all(|(a, b)| key(a) <= key(b))
}

pub fn vec_sorted_insert<T, K: PartialOrd>(v: &mut Vec<T>, elem: T, key: impl Fn(&T) -> K) {
    debug_assert!(
        slice_is_sorted(&v, &key),
        "can't perform sorted insert on empty vector"
    );

    // binary search
    let mut low_index = 0;
    let mut high_index = v.len();
    let val = key(&elem);
    while low_index < high_index {
        let mid_index = (low_index + high_index) / 2;
        if key(&v[mid_index]) < val {
            low_index = mid_index + 1;
        } else {
            high_index = mid_index;
        }
    }
    v.insert(low_index, elem);
}

#[cfg(test)]
mod test {
    #[test]
    fn is_sorted() {
        assert!(super::slice_is_sorted(&[1, 2, 3, 4], |n| n * 2));
        assert!(!super::slice_is_sorted(&[1, 3, 2, 4], |n| n * 2));
    }

    #[test]
    fn sorted_insert() {
        let mut sorted = vec![1, 2, 3, 5];
        super::vec_sorted_insert(&mut sorted, 4, |n| n * 2);
        assert_eq!(sorted, &[1, 2, 3, 4, 5]);
    }
}
