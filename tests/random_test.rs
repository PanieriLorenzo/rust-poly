//! Exploratory tests which use randomized test cases

use std::iter;

use fastrand::Rng;
use num::complex::Complex64;
use rust_poly::__util::testing::{
    check_roots, test_case_roots, RandStreamC64Cartesian, RandStreamC64Polar, RandStreamR64,
};

#[test]
fn test_uniform_real_roots() {
    let case = |deg, seed| {
        let mut seed_stream = Rng::with_seed(seed);
        let mut roots_stream = RandStreamR64::new(seed_stream.u64(..), -10.0, 10.0);
        let mut scale_stream = RandStreamR64::new(seed_stream.u64(..), 0.1, 10.0);
        for i in 0..1000 {
            let (poly, mut expected_roots) =
                test_case_roots(&mut roots_stream, &mut scale_stream, deg);
            let mut roots = poly.try_roots(1E-14, 1000, 10, None, None, None).unwrap();
            assert!(
                check_roots(&mut roots, &mut expected_roots, 0.1),
                "{:?} != {:?} @ iter = {}",
                roots,
                expected_roots,
                i
            );
        }
    };

    case(1, 1);
    case(2, 2);
    case(3, 3);
}
