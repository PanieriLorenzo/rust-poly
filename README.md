
# rust-poly
[<img alt="Crates.io" src="https://img.shields.io/crates/v/rust-poly">](https://crates.io/crates/rust-poly)

Numeric manipulation of real and complex polynomials.

> Note: this crate is still in development and might change in the future.

Basic Goals:
- [x] addition, subtraction and multiplication of complex and real univariate polynomials
- [x] long division of complex and real univariate polynomials
- [x] finding complex roots of polynomials
- [x] indexing, slicing and iterating
- [x] from/into traits
- [ ] cross-compatibility with `num` types and others
    - [x] primitive floats
    - [x] Complex
    - [ ] Ratio
    - [ ] BigFloat
    - [ ] easily implementable traits for custom types

Future Goals:
- [ ] Extremely fast polynomial evaluation
(using parallelism and SIMD)
- [ ] Fast fixed-size polynomials
- [ ] Generating important polynomial sequences
    - [x] Chebyshev type 1 polynomials
    - [ ] Chebyshev type 2 polynomials
    - [x] Bessel polynomials
    - [x] Reverse Bessel polynomials
    - [ ] Hermite polynomials
    - [ ] Lagrange polynomials
    - [ ] Bezier polynomials
    - [x] Legendre polynomials
    - more to come...
- [ ] Random integration
- [ ] Real polynomial type
    - [ ] Real polynomial factoring
- [ ] Rational functions
    - [ ] Simplification
- [ ] Multivariate polynomials
- [ ] Interpolation
- [ ] Integer polynomials
- [ ] Stabilize API
- [ ] `no_std` support
- [ ] Rayon support
- [ ] Fixed point support
- [ ] SIMD support
- [ ] Make it go fast (fastest polynomial root finder?)
    - [ ] use GCD method for determining multiplicity of roots

Non-Goals:
- [ ] Symbolic polynomial manipulation (use a symbolic algebra crate)
- [ ] LAPACK or BLAS integration (compiling shared libraries in a portable way is a pain and I don't want to do it. You're welcome to contribute)

## Contributing & Development

If you want to contribute to this project, please read [CONTRIBUTING.md](CONTRIBUTING.md).

Before opening a pull request, make sure you have read [DEVELOPMENT.md](DEVELOPMENT.md).

## Licensing

This library is covered by the MIT license, see [LICENSE](LICENSE).

Parts of the source code are based on the [NumPy](https://github.com/numpy/numpy) library for Python, used in accordance to the original license, see [licenses/numpy/LICENSE.txt](licenses/numpy/LICENSE.txt).

The `__utils/linalg.rs` module is based on the [Rulinalg](https://github.com/AtheMathmo/rulinalg) crate, used in accordance to the original license, see [licenses/rulinalg/LICENSE.md](licenses/rulinalg/LICENSE.md).

The `poly/eval.rs` module is based on the [fast_polynomial](https://crates.io/crates/fast_polynomial) crate, used in accordance to the original license, see [licenses/fast_polynomial/LICENSE](licenses/fast_polynomial/LICENSE).
