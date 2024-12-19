#!/usr/bin/env nu

def main [--fix] {
    if $fix {
        # source fixes
        cargo +nightly fmt
        cargo clippy --fix --allow-dirty --allow-staged --all-features
        cargo +nightly fmt

        # cargo
        cargo deny check
        cargo machete --fix
        cargo audit --deny warnings
        cargo features prune
    } else {
        # source fixes
        cargo clippy --all-features

        # cargo
        cargo deny check
        cargo machete
        cargo audit --deny warnings
    }
}
