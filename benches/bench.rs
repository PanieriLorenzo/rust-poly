use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use rust_poly::{
    Poly64,
    __testing::{PolyStream, RandStreamC64Polar},
    roots::newton_deflate,
};

criterion_main!(/*micro_benches, */ realistic_benches, solver_benches);
criterion_group!(micro_benches, bessel, reverse_bessel, legendre,);

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
                        .roots(1E-14, 100),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(solver_benches, newton_complex_uniform);

fn uniform_cases(seed: u64, degree: usize) -> PolyStream<f64> {
    let roots = RandStreamC64Polar::new(seed, 0.001, 1000.0, 0.0, 1.0);
    PolyStream::new(degree, roots)
}

pub fn newton_complex_uniform(c: &mut Criterion) {
    const SEED: u64 = 1;
    const EPSILON: f64 = 1E-14;
    const MAX_ITER: usize = 1000;
    let mut group = c.benchmark_group("newton_complex_uniform");
    for degree in [1, 2, 3, 7, 15, 31, 63, 127, 255] {
        group.bench_function(BenchmarkId::from_parameter(degree), |b| {
            let mut cases = uniform_cases(SEED + degree as u64, degree);
            b.iter_batched(
                || cases.next().unwrap().1,
                |mut poly| {
                    black_box(newton_deflate(
                        black_box(&mut poly),
                        Some(EPSILON),
                        Some(MAX_ITER),
                        &[],
                    ))
                },
                BatchSize::SmallInput,
            )
        });
    }
}
