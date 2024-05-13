use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rust_poly::{poly, Poly, Poly64, __util::casting::usize_to_scalar};

criterion_main!(micro_benches, realistic_benches);
criterion_group!(
    micro_benches,
    bessel,
    reverse_bessel,
    legendre,
    bench_usize_to_scalar
);

pub fn bessel(c: &mut Criterion) {
    let mut group = c.benchmark_group("bessel");
    for n in [1, 2, 4, 8, 16, 32, 64, 128] {
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| black_box(Poly64::bessel(black_box(n))))
        });
    }
    group.finish();
}

pub fn reverse_bessel(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_bessel");
    for n in [1, 2, 4, 8, 16, 32, 64, 128] {
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| black_box(Poly64::reverse_bessel(black_box(n))))
        });
    }
    group.finish();
}

pub fn bench_usize_to_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("usize_to_scalar");
    for n in [1, 10, 20, 50, 100, 200, 500] {
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| black_box(usize_to_scalar::<f64>(black_box(n))))
        });
    }
    group.finish();
}

pub fn legendre(c: &mut Criterion) {
    let mut group = c.benchmark_group("legendre");
    for n in [1, 2, 4, 8, 16, 24, 30] {
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| black_box(Poly64::legendre(black_box(n))))
        });
    }
}

criterion_group!(realistic_benches, bessel_filter_design);

pub fn bessel_filter_design(c: &mut Criterion) {
    let mut group = c.benchmark_group("bessel filter design");
    for n in [1, 2, 4, 8, 16, 32, 64, 128] {
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| {
                black_box(
                    Poly64::reverse_bessel(black_box(n))
                        .unwrap()
                        .try_n_roots(n, None, 1E-14, 100, None),
                )
            })
        });
    }
    group.finish();
}
