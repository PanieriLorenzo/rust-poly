// internal utilities for dealing with Array annoyances
// (mostly missing NumPy functions)
use ndarray::{stack, Array, Array1, Array2};
use num_traits::Zero;

pub(crate) fn np_diag<T: Clone + Zero>(diagonal: Array1<T>, offset: isize) -> Array2<T> {
    let diag_matrix: Array2<T> = Array2::<T>::from_diag(&diagonal);
    if offset == 0 {
        return diag_matrix;
    }
    let n = diagonal.len();
    if offset < 0 {
        let mut extended_matrix: Array2<T> = Array2::<T>::zeros([0, n]);
        // TODO: place this in function "h_stack"
        for _ in 0..offset.abs() {
            extended_matrix.push_row(Array1::<T>::zeros([n]).view());
        }

        // TODO: place this in function "v_stack"
    }
    todo!();
}
