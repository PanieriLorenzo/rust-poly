# Project

Project Description

<em>[TODO.md spec & Kanban Board](https://bit.ly/3fCwKfM)</em>

### Todo

- [ ] Utilize multishift QR algorithm (the one from `cgees` in LAPACK) from the Karen Braman, Ralph Byers and Roy Mathias papers (late 90s early 0s).

### In Progress

- [ ] Benchmarks

### Done ✓

- [x] `From` traits are ambiguous, because `Complex<T: Scalar>` could be any of `Complex<T>`, `Complex<Complex<T>>`, `Complex<Complex<Complex<T>>>` etc... Maybe make `From` more specific to concrete data-types? This also affects the `poly` macro.
- [x] **Bug** - `roots` hangs on some inputs (such as `legendre(7)`) due to [nalgebra#1291](https://github.com/dimforge/nalgebra/issues/1291).

### Backlog

- [ ] Conform to [Rust API guidelines](https://rust-lang.github.io/api-guidelines/checklist.html)
- [ ] Split further into `int` and `real` impls for better organization, when you get to integer polynomial algorithms.