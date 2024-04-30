
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

## Licensing

This library is covered by the MIT license, see [LICENSE](LICENSE).

Parts of the source code are based on the [NumPy](https://github.com/numpy/numpy) library for Python, used in accordance to the original license, see [licenses/numpy/LICENSE.txt](licenses/numpy/LICENSE.txt).
