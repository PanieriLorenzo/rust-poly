#![allow(clippy::unreadable_literal)]
#![allow(unused)]

// TODO: clean up these, as some of the things Herbie does are completely superfluous
use num::complex::Complex64;

/// Original: `rb * rb - ib * ib - 4.0 * ra * rc + 4.0 * ia * ic`
///
/// Optimized with Herbie 2.1
/// - Assuming all variables are within [-1e100, 1e100]
/// - Accuracy 100% -> 100% (no change)
/// - Speed 1.3x
/// - `FPCore` output: `(fma (- rb ib) (+ ib rb) (* (fma (- ra) rc (* ia ic)) 4.0)))`
///
/// ## Reproduce
/// ```racket
/// herbie shell --seed 1
/// (FPCore (rb ib ra rc ia ic)
/// :name "rb * rb - ib * ib - 4.0 * ra * rc + 4.0 * ia * ic"
/// :precision binary64
/// :pre (and (and (and (and (and (and (<= -1e+100 rb) (<= rb 1e+100)) (and (<= -1e+100 ib) (<= ib 1e+100))) (and (<= -1e+100 ra) (<= ra 1e+100))) (and (<= -1e+100 rc) (<= rc 1e+100))) (and (<= -1e+100 ia) (<= ia 1e+100))) (and (<= -1e+100 ic) (<= ic 1e+100)))
/// (+ (- (- (* rb rb) (* ib ib)) (* (* 4.0 ra) rc)) (* (* 4.0 ia) ic)))
/// ```
#[inline]
fn quad_discriminant_re(ra: f64, ia: f64, rb: f64, ib: f64, rc: f64, ic: f64) -> f64 {
    f64::mul_add(rb - ib, ib + rb, f64::mul_add(-ra, rc, ia * ic) * 4.0f64)
}

/// Original: `2.0 * rb * ib - 4.0 * ra * ic - 4.0 * ia * rc`
///
/// Optimized with Herbie 2.1
/// - Assuming all variables are within [-1e100, 1e100]
/// - Accuracy 100% -> 100% (no change)
/// - Speed 1.3x
/// - `FPCore` output: `(fma (* 2.0 rb) ib (* (fma rc ia (* ra ic)) -4.0)))`
///
/// ## Reproduce
/// ```racket
/// herbie shell --seed 1
/// (FPCore (rb ib ra ic ia rc)
/// :name "2.0 * rb * ib - 4.0 * ra * ic - 4.0 * ia * rc"
/// :precision binary64
/// :pre (and (and (and (and (and (and (<= -1e+100 rb) (<= rb 1e+100)) (and (<= -1e+100 ib) (<= ib 1e+100))) (and (<= -1e+100 ra) (<= ra 1e+100))) (and (<= -1e+100 ic) (<= ic 1e+100))) (and (<= -1e+100 ia) (<= ia 1e+100))) (and (<= -1e+100 rc) (<= rc 1e+100)))
/// (- (- (* (* 2.0 rb) ib) (* (* 4.0 ra) ic)) (* (* 4.0 ia) rc)))
/// ```
fn quad_discriminant_im(ra: f64, ia: f64, rb: f64, ib: f64, rc: f64, ic: f64) -> f64 {
    f64::mul_add(2.0f64 * rb, ib, f64::mul_add(rc, ia, ra * ic) * -4.0f64)
}

/// See `fpcore/functions/re_quad_plus.txt`
fn quad_plus_re(ra: f64, ia: f64, rb: f64, ib: f64, re_sqrt_delta: f64, im_sqrt_delta: f64) -> f64 {
    let t_0: f64 = -0.5f64 * (f64::mul_add((rb - re_sqrt_delta) / ia, ra, ib - im_sqrt_delta) / ia);
    let tmp: f64;
    if ia <= -5900000.0f64 {
        tmp = t_0;
    } else if ia <= -6.6e-81f64 {
        tmp = ((ia * (ib - im_sqrt_delta)) + (ra * (rb - re_sqrt_delta)))
            / (-2.0f64 * ((ia * ia) + (ra * ra)));
    } else if ia <= 1.6e+102f64 {
        tmp = f64::mul_add(
            ((ib - im_sqrt_delta) * ia) / ra,
            -0.5f64,
            (re_sqrt_delta - rb) * 0.5f64,
        ) / ra;
    } else {
        tmp = t_0;
    }
    tmp
}

/// See: `fpcore/functions/quad_plus_im.txt`
fn quad_plus_im(ra: f64, ia: f64, rb: f64, ib: f64, re_sqrt_delta: f64, im_sqrt_delta: f64) -> f64 {
    let t_0: f64 = f64::mul_add(
        ib - im_sqrt_delta,
        (ra / ia) * -0.5f64,
        0.5f64 * (rb - re_sqrt_delta),
    ) / ia;
    let tmp: f64;
    if ia <= -1.66e+33f64 {
        tmp = t_0;
    } else if ia <= -2.65e-79f64 {
        tmp = ((ia * (rb - re_sqrt_delta)) - (ra * (ib - im_sqrt_delta)))
            / (2.0f64 * ((ia * ia) + (ra * ra)));
    } else if ia <= 9e+101f64 {
        tmp = f64::mul_add(
            ((rb - re_sqrt_delta) * ia) / ra,
            -0.5f64,
            (ib - im_sqrt_delta) * 0.5f64,
        ) / -ra;
    } else {
        tmp = t_0;
    }
    tmp
}

fn quad_minus_re(
    ra: f64,
    ia: f64,
    rb: f64,
    ib: f64,
    re_sqrt_delta: f64,
    im_sqrt_delta: f64,
) -> f64 {
    let tmp: f64;
    if ia <= -5900000.0f64 {
        tmp = -0.5f64 * (f64::mul_add((re_sqrt_delta + rb) / ia, ra, im_sqrt_delta + ib) / ia);
    } else if ia <= -6.6e-81f64 {
        tmp = ((ia * (ib + im_sqrt_delta)) + (ra * (rb + re_sqrt_delta)))
            / (-2.0f64 * ((ia * ia) + (ra * ra)));
    } else if ia <= 9.5e+102f64 {
        tmp = (-0.5f64 / ra) * f64::mul_add(ib + im_sqrt_delta, ia / ra, rb + re_sqrt_delta);
    } else {
        tmp = -0.5f64
            * f64::mul_add(
                ra / ia,
                (rb + re_sqrt_delta) / ia,
                (ib + im_sqrt_delta) / ia,
            );
    }
    tmp
}

fn quad_minus_im(
    ra: f64,
    ia: f64,
    rb: f64,
    ib: f64,
    re_sqrt_delta: f64,
    im_sqrt_delta: f64,
) -> f64 {
    let t_0: f64 = -0.5f64 * (im_sqrt_delta + ib);
    let t_1: f64 = f64::mul_add(re_sqrt_delta + rb, 0.5f64, t_0 * (ra / ia)) / ia;
    let t_2: f64 = f64::mul_add(ra, ra, ia * ia);
    let tmp: f64;
    if ia <= -1.02e+136f64 {
        tmp = t_1;
    } else if ia <= -4.8e-78f64 {
        tmp = f64::mul_add(
            (im_sqrt_delta + ib) / t_2,
            -ra * 0.5f64,
            ((re_sqrt_delta + rb) * 0.5f64) * (ia / t_2),
        );
    } else if ia <= 9.5e+102f64 {
        tmp = f64::mul_add(re_sqrt_delta + rb, (ia / ra) * 0.5f64, t_0) / ra;
    } else {
        tmp = t_1;
    }
    tmp
}

/// Quadratic formula, a is the highest degree term
pub fn quadratic_c64(a: Complex64, b: Complex64, c: Complex64) -> (Complex64, Complex64) {
    let (ra, ia) = (a.re, a.im);
    let (rb, ib) = (b.re, b.im);
    let (rc, ic) = (c.re, c.im);
    let re_delta = quad_discriminant_re(ra, ia, rb, ib, rc, ic);
    let im_delta = quad_discriminant_im(ra, ia, rb, ib, rc, ic);
    let sqrt_delta = Complex64::new(re_delta, im_delta).sqrt();
    let (re_sqrt_delta, im_sqrt_delta) = (sqrt_delta.re, sqrt_delta.im);
    let rx1 = quad_plus_re(ra, ia, rb, ib, re_sqrt_delta, im_sqrt_delta);
    let ix1 = quad_plus_im(ra, ia, rb, ib, re_sqrt_delta, im_sqrt_delta);
    let rx2 = quad_minus_re(ra, ia, rb, ib, re_sqrt_delta, im_sqrt_delta);
    let ix2 = quad_minus_im(ra, ia, rb, ib, re_sqrt_delta, im_sqrt_delta);
    (Complex64::new(rx1, ix1), Complex64::new(rx2, ix2))
}

#[cfg(test)]
mod test {
    use num::complex::Complex64;

    use super::quadratic_c64;

    #[test]
    fn test_quadratic_c64() {
        const TOL: f64 = 0.000000000000001;
        let a = Complex64::new(3.0, 0.0);
        let b = Complex64::new(6.0, 6.0);
        let c = Complex64::new(15.0, -30.0);
        let (r1, r2) = quadratic_c64(a, b, c);
        let expected1 = Complex64::new(1.0, 2.0);
        let expected2 = Complex64::new(-3.0, -4.0);
        // note that order may be flipped, so we have to test both cases
        dbg!(r1);
        dbg!(r2);
        assert!((r1 - expected1).norm() < TOL || (r2 - expected1).norm() < TOL);
        assert!((r2 - expected2).norm() < TOL || (r1 - expected2).norm() < TOL);
    }
}
