#![allow(unused_mut, unused_parens)]

fn ex0(im_sqrt_delta: f64, ib: f64, re_sqrt_delta: f64, ra: f64, ia: f64, rb: f64) -> f64 {
	let mut t_0: f64 = -0.5f64 * (im_sqrt_delta + ib);
	let mut t_1: f64 = f64::mul_add((re_sqrt_delta + rb), 0.5f64, (t_0 * (ra / ia))) / ia;
	let mut t_2: f64 = f64::mul_add(ra, ra, (ia * ia));
	let mut tmp: f64;
	if ia <= -1.02e+136f64 {
	    tmp = t_1;
	} else if ia <= -4.8e-78f64 {
	    tmp = f64::mul_add(((im_sqrt_delta + ib) / t_2), (-ra * 0.5f64), (((re_sqrt_delta + rb) * 0.5f64) * (ia / t_2)));
	} else if ia <= 9.5e+102f64 {
	    tmp = f64::mul_add((re_sqrt_delta + rb), ((ia / ra) * 0.5f64), t_0) / ra;
	} else {
	    tmp = t_1;
	}
	tmp
}

