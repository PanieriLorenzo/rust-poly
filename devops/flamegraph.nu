#!/usr/bin/env nu

def main [example] {
    cargo clean
    CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --example $"($example)"
}