//! These tests are tuned to just barely succeede, to find what the stability
//! limits of different root finders are and to catch stability regressions.
//!
//! This should be updated if the stability is improved.

use rust_poly::{
    __util::testing::{
        check_roots, test_case_conj_roots, test_case_multiple_roots, test_case_roots,
        RandStreamC64Cartesian, RandStreamC64Polar, RandStreamR64,
    },
    roots::{naive, newton},
};

/// - max degree: 6
/// - worst-case error: 0.0385
#[test]
fn naive_real() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(1, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(2, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 6);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.0385),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

/// - max degree: 6
/// - worst-case error: 0.0385
#[test]
fn newton_real() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(1, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(2, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 6);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots.clone(), expected_roots.clone(), 0.0385),
            "@ {i}: {roots:?} != {expected_roots:?}",
        );
    }
}

/// - max degree: 6
/// - worst-case error: 0.097
#[test]
fn naive_real_multiplicity_1() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 6, 1);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.97E-1),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 6
/// - worst-case error: 0.097
#[test]
fn newton_real_multiplicity_1() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 6, 1);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.097),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 4
/// - worst-case error: 0.0015
#[test]
fn naive_real_multiplicity_3() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 4, 3);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.0015),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 5
/// - worst-case error: 0.035
#[test]
fn newton_real_multiplicity_3() {
    let _ = simple_logger::init_with_level(log::Level::Debug);
    let mut roots_stream = RandStreamR64::new(7, -2.0, 2.0);
    let mut scale_stream = RandStreamR64::new(8, 1.0, 10.0);
    for i in 0..1000 {
        let (poly, expected_roots) =
            test_case_multiple_roots(&mut roots_stream, &mut scale_stream, 5, 3);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.035),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 10
/// - worst-case error: 0.0046
#[test]
fn naive_complex() {
    let mut roots_stream = RandStreamC64Cartesian::new(3, -2.0, 2.0, -2.0, 2.0);
    let mut scale_stream = RandStreamC64Polar::new(4, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 10);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.46E-2),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 10
/// - worst-case error: 0.0046
#[test]
fn newton_complex() {
    let mut roots_stream = RandStreamC64Cartesian::new(3, -2.0, 2.0, -2.0, 2.0);
    let mut scale_stream = RandStreamC64Polar::new(4, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, 10);
        let roots = newton(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 0.46E-2),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 2
/// - worst-case error: 1E-12
#[test]
fn naive_conjugate() {
    // roots are within the unit circle, this is common in signal processing
    let mut roots_stream = RandStreamC64Polar::new(5, 0.0, 1.0, 0.0, 1.0);
    let mut scale_stream = RandStreamC64Polar::new(6, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_conj_roots(&mut roots_stream, &mut scale_stream, 2);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 1E-13),
            "@ {}: {}",
            i,
            poly
        );
    }
}

/// - max degree: 2
/// - worst-case error: 1E-12
#[test]
fn newton_conjugate() {
    // roots are within the unit circle, this is common in signal processing
    let mut roots_stream = RandStreamC64Polar::new(5, 0.0, 1.0, 0.0, 1.0);
    let mut scale_stream = RandStreamC64Polar::new(6, 1.0, 10.0, 0.0, 1.0);
    for i in 0..1000 {
        let (poly, expected_roots) = test_case_conj_roots(&mut roots_stream, &mut scale_stream, 2);
        let roots = naive(&mut poly.clone(), Some(1E-14), Some(100), &[]).unwrap();
        assert!(
            check_roots(roots, expected_roots, 1E-13),
            "@ {}: {}",
            i,
            poly
        );
    }
}
