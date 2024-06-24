#!/usr/bin/env nu

def main [--fix] {
    if $fix {
        cargo +nightly fmt
        cargo clippy --fix --allow-dirty --allow-staged --all-features
        cargo deny check
        cargo +nightly fmt
    } else {
        cargo clippy --all-features
        cargo deny check
    }
}
