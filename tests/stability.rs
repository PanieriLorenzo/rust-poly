//! These tests are tuned to just barely succeede, to find what the stability
//! limits of different root finders are and to catch stability regressions.
//!
//! This should be updated if the stability is improved.
//!
//! The goal here is thus not to benchmark the rate of convergence, but see if
//! they can converge at all for a general polynomial.
//!
//! All methods provide an error considerably higher than epsilon, because repeated
//! division used in shrinking introduces considerable rounding error, so we
//! allow a validation error of up to 0.1, above which we consider it a failure.
//!
//! All other parameters are chosen case-by-case to maximize the degree that can
//! be handled, so the main parameter of interest here is the maximum degree
//! that can be handled, but also of interest is the number of iterations
//! necessary to achieve that, and to some extent what epsilon is necessary to
//! achieve at least 0.1 validation error.

use rust_poly::{
    __util::testing::{
        check_roots, test_case_conj_roots, test_case_multiple_roots, test_case_roots,
        RandStreamC64Cartesian, RandStreamC64Polar, RandStreamR64,
    },
    roots::{halley, naive, newton},
};

const LOG_LEVEL: log::Level = log::Level::Warn;

fn test_setup() {
    let _ = simple_logger::init_with_level(LOG_LEVEL);
}

#[test]
fn naive_real() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(1, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(2, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 6);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn newton_real() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(1, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(2, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 6);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(25), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn halley_real() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(1, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(2, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 5);
        let roots = halley(&mut poly.clone(), Some(1E-11), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn naive_real_multiplicity_1() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 6, 1);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn newton_real_multiplicity_1() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 6, 1);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn halley_real_multiplicity_1() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 5, 1);
        let roots = halley(&mut poly.clone(), Some(1E-12), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn naive_real_multiplicity_3() {
    test_setup();
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 6, 3);
        let roots = naive(&mut poly.clone(), Some(1E-15), Some(200), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

//NOTE: straight up does not work
// #[test]
// fn newton_real_multiplicity_3() {
//     test_setup();
//     let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
//     let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
//     for i in 0..1000 {
//         let (poly, expected_roots) =
//             test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 4, 3);
//         let roots = newton(&mut poly.clone(), Some(1E-14), Some(50), &[]).unwrap();
//         assert!(
//             check_roots(roots.clone(), expected_roots.clone(), 0.1),
//             "@ {i}: {roots:?} != {expected_roots:?}, poly: {poly}",
//         );
//     }
// }

#[test]
fn naive_complex() {
    test_setup();
    let mut roots_stream = RandStreamC64Cartesian::new(3, -2.0, 2.0, -2.0, 2.0);
    let mut scale_stream = RandStreamC64Polar::new(4, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 10);
        let roots = naive(&mut poly.clone(), Some(1E-12), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn newton_complex() {
    test_setup();
    let mut roots_stream = RandStreamC64Cartesian::new(3, -2.0, 2.0, -2.0, 2.0);
    let mut scale_stream = RandStreamC64Polar::new(4, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 10);
        let roots = newton(&mut poly.clone(), Some(1E-12), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

#[test]
fn halley_complex() {
    test_setup();
    let mut roots_stream = RandStreamC64Cartesian::new(3, -2.0, 2.0, -2.0, 2.0);
    let mut scale_stream = RandStreamC64Polar::new(4, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 9);
        let roots = halley(&mut poly.clone(), Some(1E-10), Some(50), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.1),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

// NOTE: this does not work for degree 3 or higher, so it isn't interesting
//       (degree 2 always uses quadratic formula)
// #[test]
// fn naive_conjugate() {
//     // roots are within the unit circle, this is common in signal processing
//     let mut roots_stream = RandStreamC64Polar::new(5, 0.0, 1.0, 0.0, 1.0);
//     let mut scale_stream = RandStreamC64Polar::new(6, 1.0, 10.0, 0.0, 1.0);
//     for i in 0..1000 {
//         let (poly, expected_roots) = test_case_conj_roots(&mut roots_stream, &mut scale_stream, 2);
//         let roots = naive(&mut poly.clone(), Some(1E-2), Some(10), &[]).unwrap();
//         assert!(
//             check_roots(roots.clone(), expected_roots.clone(), 0.1),
//             "@ {i}: {roots:?} != {expected_roots:?}",
//         );
//     }
// }

// NOTE: this does not work for degree 3 or higher, so it isn't interesting
// #[test]
// fn newton_conjugate() {
//     // roots are within the unit circle, this is common in signal processing
//     let mut roots_stream = RandStreamC64Polar::new(5, 0.0, 1.0, 0.0, 1.0);
//     let mut scale_stream = RandStreamC64Polar::new(6, 1.0, 10.0, 0.0, 1.0);
//     for i in 0..1000 {
//         let (poly, expected_roots) = test_case_conj_roots(&mut roots_stream, &mut scale_stream, 2);
//         let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
//         assert!(
//             check_roots(roots, expected_roots, 1E-13),
//             "@ {}: {}",
//             i,
//             poly
//         );
//     }
// }
