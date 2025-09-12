# Changelog

All notable changes to this project will be documented in this file.

## [0.4.3] - 2025-09-12

### ⚙️ Miscellaneous Tasks

- Updated dependencies
- Internal cleanup

## [0.4.2] - 2024-10-11

### 🚀 Features

- Added `roots_expert` method with more customization
- [**breaking**] Changed signature of `from_xxx_iterator` constructors so they don't take a length parameter
- [**breaking**] Removed `initial_guesses_random` use `initial_guesses_circle` instead
- Added multiples handling to `roots` and `roots_expert`
- `Poly::constant` method for conveniently casting a constant to a polynomial

### 🚜 Refactor

- Removed deprecated internal `trivial_roots`
- Removed unused vec utils

### 📚 Documentation

- Updated `roots` and `roots_expert` docs

### ⚙️ Miscellaneous Tasks

- Updated dependencies

## [0.4.1] - 2024-08-28

### 🚀 Features

- Added `conj` method to compute the conjugate polynomial

### 🐛 Bug Fixes

- [**breaking**] Parallel not re-using given initial guesses

## [0.4.0] - 2024-08-26

### 🚀 Features

- [**breaking**] Removed dependency on `nalgebra` and related methods
- [**breaking**] Relaxed trait bounds for `RealScalar` to allow non-floats
- [**breaking**] Use f128 in `roots` for higher precision

### 🚜 Refactor

- Use `f128` type internally for higher precision

## [0.3.0] - 2024-08-24

### 🚀 Features

- [**breaking**] Redid `Poly::roots` entirely, now it's much more stable and efficient
- [**breaking**] Made the `Poly::from` trait impls less ambiguous
- [**breaking**] Changed how zero polynomials are represented
- [**breaking**] `eval_point` renamed to `eval`, `eval` renamed to `eval_multiple`
- [**breaking**] Changed `eval_multiple` to have an output buffer
- [**breaking**] Use a slice for indexing instead of re-inventing the wheel (removed `Get` trait)
- [**breaking**] Custom error type for root finders
- [**breaking**] Assert that testing api is not used in production
- New procedural root finder API (in the `roots` module)
    - `halley` root finder
    - `newton` root finder
    - `naive` root finder (Newton without variable stepsize and local-minimum detection)
    - `aberth_ehrlich` root finder
    - `parallel` root finder (runs `newton` or `halley` in parallel)
    - `initial_guess_circle`
- Relax bounds on `T` for `Poly<T>`
- Removed dependency `duplicate`
- Added `diff` and `integral`
- `translate` method to move polynomial on x and y axes
- Added `legendre` polynomials
- Impl `std::fmt::Display` for `Poly`
- Added `shift_up` and `shift_down` methods for changing degrees of polynomial
- `SafeConstants` trait for better numeric precision
- Improve precision of `eval` (at the expense of some speed)
- Added `from_real_iterator` and `from_complex_iterator`

### 🐛 Bug Fixes

- Potential bug in `poly!` macro, making nested complex numbers by accident
- Removed debug print from `bessel`
- `Poly::div_rem` was completely wrong
- `normalize` converted constant 0 to zero polynomial
- `diff` produced a zero polynomial instead of a constant zero on constant polynomials
- Math error in `Poly::companion` (closes #3)
- Roots accidentally required `na::Scalar` instead of `crate::Scalar`
- Root finders use exact solutions for linear and quadratic polynomials
- `poly` macro now works without having to import `num::Zero`
- `div_rem` stuck in endless loop if divisor is 1
- Multiple bugs preventing trivial polynomial roots to be found
- Removed old test relying on `plotters`, now that we use `plotly` instead
- Over 60 lints and warnings
- Cargo clippy fixes

### 📚 Documentation

- Moved dev docs from `README.md` to `DEVELOPMENT.md`
- Clarified behavior of zero polynomial to the power of zero
- Many others

### ⚡ Performance

- Added criterion benchmarking
- Made `eval` faster
- Various micro optimizations

## [0.2.0] - 2024-05-04

### 🚀 Features

- Impl `IntoIterator` for `&Poly` and `&mut Poly`

### 🐛 Bug Fixes

- [**breaking**] Clippy warnings

### 📚 Documentation

- Added developer documentation

### ⚙️ Miscellaneous Tasks

- Update itertools to 0.12.1
- Updated TODOs
- Added commitlint to pre-commit hooks
- Added git-cliff dev tool
- Added old changes to CHANGELOG.md

## [0.1.13] - 2023-11-01

### Changes
- Add Bessel polynomials
- Add revere bessel polynomials

### Bug Fixes
- `Poly::term()` having incorrect coefficients

## [0.1.12] - 2023-09-08

### Changes

- Div_rem now returns anyhow::Result
- Switched from num_complex and num_traits to num crate
- Div_rem takes ownership of both arguments
- Reduced cloning in ops
- More pre-computed chebyshev polynomials
- Div for poly and complex
- More combinations of borrowed and owned variants of operators
- More combinations of borrowed and owned operators
- Almost_zero method
- Checked_div and checked_rem
- CheckedDiv trait
- All combinations of owned and borrowed ops
- Missing cargo.lock
- Checked and negative indexing
- More concise macro syntax for complex values with no imaginary component

<!-- generated by git-cliff -->
