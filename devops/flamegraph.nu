#!/usr/bin/env nu

def main [] {
    let level = (sysctl kernel.perf_event_paranoid | parse --regex '(?<level>-?\d)').level.0
    sudo sysctl $"kernel.perf_event_paranoid=0"
    try {
        cargo clean
        # you have to ask cargo very very nicely to _actually_ run benchmarks as
        # benchmarks and not as tests...
        CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --bench bench -- --bench --profile-time 10
    }
    sudo sysctl $"kernel.perf_event_paranoid=($level)"
}