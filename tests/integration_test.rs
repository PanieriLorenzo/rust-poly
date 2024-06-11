use rust_poly::{complex, num::CheckedDiv, poly, Poly64};

/// Check that zero polynomials behave as expected
#[test]
fn zero_polynomial_properties() {
    let zp: Poly64 = poly![];
    assert_eq!(zp, poly![0.0]);
    assert_eq!(zp, poly![(0.0, 0.0)]);

    // identity of addition
    assert_eq!(&zp + poly![1.0, 2.0, 3.0], poly![1.0, 2.0, 3.0]);
    assert_eq!(poly![1.0, 2.0, 3.0] + &zp, poly![1.0, 2.0, 3.0]);

    // absorbing element of multiplication
    assert_eq!(&zp * poly![1.0, 2.0, 3.0], zp);

    // division by zero should fail
    assert!(poly![1.0, 2.0, 3.0].checked_div(&zp).is_none());

    // should yield zero everywhere
    assert_eq!(zp.eval_point(complex!(100.0, -100.0)), complex!(0.0, 0.0));

    // composition with any polynomial should yield zero
    assert_eq!(zp.clone().compose(poly![1.0, 2.0, 3.0]), zp);

    // composition of any polynomial with zero should yield leading coefficient
    assert_eq!(poly![1.0, 2.0, 3.0].compose(zp.clone()), poly![1.0]);

    // derivative of a constant should be zero
    assert_eq!(poly![1.0].diff(), zp);

    // zero to any power should be zero
    assert_eq!(zp.clone().pow(10), zp);

    // zero to the zero power is defined to be equal to 1 by convention
    assert_eq!(zp.pow(0), poly![1.0])
}
