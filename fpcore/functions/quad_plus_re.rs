#![allow(unused_mut, unused_parens)]

fn ex0(im_sqrt_delta: f64, ib: f64, re_sqrt_delta: f64, ra: f64, ia: f64, rb: f64) -> f64 {
	let mut t_0: f64 = -0.5f64 * (f64::mul_add(((rb - re_sqrt_delta) / ia), ra, (ib - im_sqrt_delta)) / ia);
	let mut tmp: f64;
	if ia <= -5900000.0f64 {
	    tmp = t_0;
	} else if ia <= -6.6e-81f64 {
	    tmp = ((ia * (ib - im_sqrt_delta)) + (ra * (rb - re_sqrt_delta))) / (-2.0f64 * ((ia * ia) + (ra * ra)));
	} else if ia <= 1.6e+102f64 {
	    tmp = f64::mul_add((((ib - im_sqrt_delta) * ia) / ra), -0.5f64, ((re_sqrt_delta - rb) * 0.5f64)) / ra;
	} else {
	    tmp = t_0;
	}
	tmp
}

