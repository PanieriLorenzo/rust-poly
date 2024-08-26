use num::Float;

use super::F128;

type F = F128;

const fn abs_diff(a: i32, b: i32) -> u32 {
    a.wrapping_sub(b).wrapping_abs() as u32
}

pub(crate) fn pow(a: F, b: i32) -> F {
    let mut a = a;
    let recip = b < 0;
    let mut pow = abs_diff(b, 0);
    let mut mul = F::ONE;
    loop {
        if (pow & 1) != 0 {
            mul = mul * a;
        }
        pow >>= 1;
        if pow == 0 {
            break;
        }
        a = a * a;
    }

    if recip {
        mul.recip()
    } else {
        mul
    }
}
