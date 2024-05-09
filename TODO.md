# Project

Project Description

<em>[TODO.md spec & Kanban Board](https://bit.ly/3fCwKfM)</em>

### Todo

- [ ] Split further into `int` and `real` impls for better organization, when you get to integer polynomial algorithms.
- [ ] Benchmarks

### Bugs


### In Progress


### Done ✓

- [x] `From` traits are ambiguous, because `Complex<T: Scalar>` could be any of `Complex<T>`, `Complex<Complex<T>>`, `Complex<Complex<Complex<T>>>` etc... Maybe make `From` more specific to concrete data-types? This also affects the `poly` macro.
