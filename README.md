
# Rust Poly
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
- [ ] Ultra-high degree factorization using Lindsey-Fox
- [ ] Fast fixed-size polynomials
- [ ] More robust eigenvalue-based root finder using LAPACK bindings
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

## Contributing & Development

If you want to contribute to this project, please read [CONTRIBUTING.md](CONTRIBUTING.md).

Before opening a pull request, make sure you have read [DEVELOPMENT.md](DEVELOPMENT.md).

## Licensing

This library is covered by the MIT license, see [LICENSE](LICENSE).

Parts of the source code are based on the [NumPy](https://github.com/numpy/numpy) library for Python, used in accordance to the original license, see [licenses/numpy/LICENSE.txt](licenses/numpy/LICENSE.txt).

## References
- \[[Vestermark 2023](http://dx.doi.org/10.13140/RG.2.2.30423.34728)\]: Vestermark, Henrik. (2023). _Practical Implementation of Polynomial Root Finders._ DOI: 10.13140/RG.2.2.30423.34728.
- \[[McNamee 2005](https://www.researchgate.net/publication/228745231_A_comparison_of_a_priori_bounds_on_real_or_complex_roots_of_polynomials)\] Mcnamee, J. and Olhovsky M. (2005) _A comparison of a priori bounds on (real or complex) roots of polynomials._
- \[[McNamee 2007 I](https://shop.elsevier.com/books/numerical-methods-for-roots-of-polynomials-part-i/mcnamee/978-0-444-52729-5)\]: McNamee, J.M. _Numerical Methods for Roots of Polynomials - Part I._ ISBN: 978-0-444-52729-5.
- \[[McNamee 2007 II](https://shop.elsevier.com/books/numerical-methods-for-roots-of-polynomials-part-ii/mcnamee/978-0-444-52730-1)\]: McNamee, J.M. _Numerical Methods for Roots of Polynomials - Part II._ ISBN: 978-0-444-52730-1.
- \[Lagouanelle 1966\]: Lagouanelle, J.L. (1966) _Sur Une Metode de Calcul de l’Ordre de Multiplicite des Zeros d’Un Polynome._ Comptes Rendus de l'Académie des Sciences, 262, 626-627.
- \[[Madsen 1973](https://doi.org/10.1007/BF01933524)\]: Madsen, K. (1973) _A root-finding algorithm based on Newton's method._ BIT 13, 71–75. DOI: 10.1007/BF01933524.
- \[[Bini 1996](https://doi.org/10.1007/BF02207694)\]: Bini, D. A. (1996) _Numerical computation of polynomial zeros by means of Aberth’s method._ Numer Algor, vol. 13, no. 2, pp. 179–200. DOI: 10.1007/BF02207694.

