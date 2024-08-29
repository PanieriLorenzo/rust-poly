//! Exploratory tests which use randomized test cases

use fastrand::Rng;
use rust_poly::__testing::{check_roots, test_case_roots, RandStreamR64};

const LOG_LEVEL: log::Level = log::Level::Trace;

fn test_setup() {
    let _ = simple_logger::init_with_level(LOG_LEVEL);
}

#[test]
fn test_uniform_real_roots() {
    test_setup();
    let case = |deg, seed, tolerance_min, tolerance_max| {
        let mut seed_stream = Rng::with_seed(seed);
        let mut roots_stream = RandStreamR64::new(seed_stream.u64(..), -10.0, 10.0);
        let mut scale_stream = RandStreamR64::new(seed_stream.u64(..), 0.1, 10.0);
        for i in 0..1 {
            let (poly, expected_roots) = test_case_roots(&mut roots_stream, &mut scale_stream, deg);
            let roots = poly.roots(tolerance_min, 1000).unwrap();
            assert!(
                check_roots(roots.clone(), expected_roots.clone(), tolerance_max),
                "{:?} != {:?} @ iter = {}",
                roots,
                expected_roots,
                i
            );
        }
    };

    // case(3, 1, 1E-22, 1E-10);
    // case(4, 2, 1E-20, 1E-8);
    // case(5, 3, 1E-20, 1E-9);
    // case(6, 4, 1E-17, 1E-7);
    // case(7, 5, 1E-17, 1E-5);
    // case(8, 6, 1E-16, 1E-6);
    // case(9, 7, 1E-13, 1E-5);
    case(10, 8, 1E-15, 1E-5);
}
