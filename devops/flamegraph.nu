#!/usr/bin/env nu

def main [--bench, --test] {
    # swap the current paranoid level with 0
    let level = (sysctl kernel.perf_event_paranoid | parse --regex '(?<level>-?\d)').level.0
    sudo sysctl $"kernel.perf_event_paranoid=0"

    if $bench {
        try {
            cargo clean
            # you have to ask cargo very very nicely to _actually_ run benchmarks as
            # benchmarks and not as tests...
            CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --bench bench -- --bench --profile-time 10
        }
    }

    if $test {
        try {
            cargo clean
            CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --test random
        }
    }

    # restore old paranoid value
    sudo sysctl $"kernel.perf_event_paranoid=($level)"
}