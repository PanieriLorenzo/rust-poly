#!/usr/bin/env nu
#
# build and publish a new version

use std log

def thin-context [] {
    git remote update
    mut context = {
        version: (git-cliff --bumped-version),
        origin: (git for-each-ref --format='%(upstream:short)' $"(git symbolic-ref -q HEAD)"),
    }
    $context.on_main = $context.origin == "origin/main" or $context.origin == "origin/master"
    let tags = git tag | lines
    let tags = $tags | where $it =~ '^v(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$'
    $context.version_old = ($tags | reverse).0
    $context.up_to_date = true              # a lie
    $context.dirty = false                  # a lie
    return $context
}

def gather-context [] {
    log info "gathering context"
    mut context = thin-context
    let local_commit = git rev-parse @
    let remote_commit = git rev-parse @{u}
    $context.up_to_date = $local_commit == $remote_commit
    let res = git diff HEAD --quiet | complete
    $context.dirty = $res.exit_code != 0
    return $context
}

def pre-flight-checks [context] {
    log info "running pre-flight checks"
}


def main [
    --force (-f)            # ignore most errors and continue with defaults
    --interactive (-i)      # prompt on errors instead of failing
    --dry (-d)              # do not apply any changes
] {
    log debug "validating arguments"
    if $force and $interactive {
        log error "cannot use --force and --interactive together"
        exit 1
    }

    let context = gather-context
}
