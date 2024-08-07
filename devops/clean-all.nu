#!/usr/bin/env nu

# cleanup temp file generated by tests
rm -rf temp

# cleanup node modules
rm -rf node_modules
npm cache verify

# cleanup husky files
rm -rf "--help" "--verbose" ".husky/_"

# cleanup cargo
cargo clean

# cleanup flamegraph
rm -f flamegraph.svg
rm -f perf.data
rm -f perf.data.old
