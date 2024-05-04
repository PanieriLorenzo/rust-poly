use num::cast;

pub fn usize_to_i32(x: usize) -> i32 {
    x.try_into().unwrap_or(i32::MAX)
}

pub fn usize_to_u32(x: usize) -> u32 {
    x.try_into().unwrap_or(u32::MAX)
}

pub fn usize_to_f64(x: usize) -> f64 {
    cast(x).unwrap_or(f64::INFINITY)
}
