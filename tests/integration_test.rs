use std::f64::consts::PI;

use num::{complex::ComplexFloat, Complex, Zero};
use rand::prelude::*;
use rust_poly::{complex, poly, Poly, Poly64};

fn gen_rand_roots(rng: &mut impl Rng) -> Vec<Complex<f64>> {
    const MAX_LEN: u32 = 7;
    let len = rng.next_u32() % MAX_LEN;
    let mut v: Vec<Complex<f64>> = vec![];
    for _ in 0..len {
        let re = rng.gen::<f64>();
        let im = rng.gen::<f64>();
        v.push(complex!(re, im));
    }
    v
}

/// a tour of the most important features and operators, in random combinations
#[test]
fn stress_test_1() {
    const SEED: u64 = 0;
    const MAX_POW: u32 = 3;
    let mut rng = StdRng::seed_from_u64(SEED);
    const ITER: usize = 250;
    for _ in 0..ITER {
        let roots1 = Poly::from_roots(&gen_rand_roots(&mut rng));
        let roots2 = Poly::from_roots(&gen_rand_roots(&mut rng));
        let lhs = if rng.gen_bool(0.5) {
            roots1.clone()
        } else {
            roots2.clone()
        };
        let rhs = if rng.gen_bool(0.5) { roots2 } else { roots1 };
        let ops = rng.next_u32() % 5;
        dbg!(ops);
        let res = match ops {
            0 => lhs + rhs,
            1 => lhs - rhs,
            2 => lhs * rhs,
            3 => {
                lhs / if rhs.is_zero() {
                    continue;
                } else {
                    rhs
                }
            }
            4 => lhs.compose(rhs),
            _ => unreachable!(),
        };
        let res = res.pow(rng.next_u32() % MAX_POW);
        let _ = res.try_roots(0.1, 100, 10, None, None, None);
    }
}

#[test]
fn back_and_forth_1() {
    const EPSILON: f64 = 1E-14;

    let p = poly![2.0, -3.0, 4.0, 1.0];

    // because p is monic, we expect pp to be almost identical
    let pp = Poly::from_roots(
        p.try_roots(EPSILON, 100, 10, None, None, None)
            .unwrap()
            .as_slice(),
    );

    // assert almost equal
    for i in 0..p.len() {
        assert!((p[i] - pp[i]).abs() < EPSILON);
    }
}

#[test]
fn big_roots() {
    assert!(Poly64::bessel(85)
        .unwrap()
        .try_roots(1E-14, 1000, 10, None, None, None)
        .unwrap()
        .iter()
        .all(|z| !z.is_nan()));
    assert!(Poly64::reverse_bessel(50)
        .unwrap()
        .try_roots(1E-14, 1000, 10, None, None, None)
        .unwrap()
        .iter()
        .all(|z| !z.is_nan()));
    assert!(Poly64::legendre(100)
        .try_roots(1E-14, 1000, 10, None, None, None)
        .unwrap()
        .iter()
        .all(|z| !z.is_nan()));
}
