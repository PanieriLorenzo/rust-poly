//! Exploratory tests which use randomized test cases

use fastrand::Rng;
use rust_poly::__testing::{check_roots, test_case_roots, RandStreamR64};

#[test]
fn test_uniform_real_roots() {
    let case = |deg, seed, tolerance_min, tolerance_max| {
        let mut seed_stream = Rng::with_seed(seed);
        let mut roots_stream = RandStreamR64::new(seed_stream.u64(..), -10.0, 10.0);
        let mut scale_stream = RandStreamR64::new(seed_stream.u64(..), 0.1, 10.0);
        for i in 0..10000 {
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

    case(3, 1, 1E-19, 1E-9);
    case(4, 2, 1E-15, 1E-7);
    case(5, 3, 1E-15, 1E-7);
    case(6, 4, 1E-15, 1E-6);
}
