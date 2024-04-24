
# rust-poly
[<img alt="Crates.io" src="https://img.shields.io/crates/v/rust-poly">](https://crates.io/crates/rust-poly)

Numeric manipulation of real and complex polynomials.

> Note: this crate is still in development and might change in the future.

Basic Goals:
- [x] addition, subtraction and multiplication of complex and real univariate polynomials
- [x] long division of complex and real univariate polynomials
- [x] finding complex roots of polynomials
- [ ] indexing, slicing and iterating
    - [x] indexing by term
    - [ ] indexing by coefficient
    - [x] iterating by coefficient
    - [ ] iterating by term
    - [ ] slicing coefficients by range
    - [ ] slicing terms by range
    - [ ] IntoIterator traits
- [x] from/into traits

Future Goals:
- [ ] Generating important polynomial sequences
    - [x] Chebyshev type 1 polynomials
    - [ ] Chebyshev type 2 polynomials
    - [x] Bessel polynomials
    - [x] Reverse Bessel polynomials
    - [ ] Hermite polynomials
    - [ ] Lagrange polynomials
    - more to come...
- [ ] Serde integration
- [ ] Random integration
- [ ] Real polynomial type
    - [ ] Real polynomial factoring
- [ ] Rational functions
    - [ ] Simplification

## Licensing

This library is covered by the MIT license, see [LICENSE](LICENSE).

Parts of the source code are based on the [NumPy](https://github.com/numpy/numpy) library for Python, used in accordance to the original license, see [licenses/numpy/LICENSE.txt](licenses/numpy/LICENSE.txt).
