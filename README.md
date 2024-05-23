
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

## Development

We don't do GitHub actions over here, because I'm a solo developer and actions are expensive.

Automation is performed using nushell scripts found in `devops/`

### Conventional Commits

Commit messages must adhere to the [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) specification and additionally to the [commitlint conventional config](https://github.com/conventional-changelog/commitlint/tree/master/%40commitlint/config-conventional).

**Breaking Changes**:
- Any other type that ends in a bang, like `feat!` or `fix!`. Note that `feat!` can also be used to mark removal of features.

**Minor Changes**:
- `feat`: any new functionality that is visible to the user

**Patch Changes**:
- `fix`: a bugfix
- `docs`: documentation changes
- `build`: changes to build scripts
- `ci`: changes to CI pipelines
- `test`: adding tests or benchmarks
- `perf`: changes that affect performance
- `refactor`: major refactoring
- `revert`: reverting a change
- `style`: stylistic changes
- `chore`: any changes that are mostly administrative, e.g. small refactors, code style, comments, semver adjustments, etc...

### Changelog

The changelog is generated automatically from commit messages. During the publishing process of a new release, the generated changelog can be manually edited to include additional information or rephrase the changes.

Right now, changelogs are generated with [git-cliff](https://github.com/orhun/git-cliff). Configurations for how this happens are in `cliff.toml`.

### Releasing a New Version

TODO

### Benchmarking

```bash
cargo flamegraph --bench bench -- <name-of-bench> --bench
```

## Licensing

This library is covered by the MIT license, see [LICENSE](LICENSE).

Parts of the source code are based on the [NumPy](https://github.com/numpy/numpy) library for Python, used in accordance to the original license, see [licenses/numpy/LICENSE.txt](licenses/numpy/LICENSE.txt).
