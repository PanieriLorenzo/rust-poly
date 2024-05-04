#!/usr/bin/env nu
#
# build and publish a new version

use std log
use _lib.nu yes-or-no

const COMMIT_TYPES = [__unknown, feat, fix, docs, build, ci, test, perf, refactor, revert, style, chore]
const FORCE_PROMPT = "Run with --force to continue anyway."

def pre-flight-checks [conf] {
    log info "running pre-flight checks"
    log debug "ensuring remote branch is main or master"
    let origin = git for-each-ref --format='%(upstream:short)' $"(git symbolic-ref -q HEAD)"
    if $origin != "origin/main" and $origin != "origin/master" {
        log critical "current branch is not main/master, cannot release from non-main/master branches"
        exit 1
    }

    git remote update
    let local_commit = git rev-parse @
    let remote_commit = git rev-parse @{u}
    if $local_commit != $remote_commit {
        log critical "local branch is not in sync with remote, can only release if local commit matches remote. Push/pull to resolve."
        exit 1
    }

    let res = git diff HEAD --quiet | complete
    if $res.exit_code != 0 {
        let err_msg = "There are local uncommitted changes that wouldn't be included."
        if not $conf.force and not $conf.interactive {
            log error $"($err_msg) ($FORCE_PROMPT)"
            exit 1
        }
        log warning $err_msg
        if $conf.interactive and not $conf.force {
            if not (yes-or-no $"($err_msg) Continue anyway?") {
                exit 0
            }
        }
    }

    log info "pre-flight checks OK"
}

def find-latest-tag [] {
    let tags = git tag | lines

    # filter for semver tags only
    let tags = $tags | where $it =~ '^v(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$'
    return ($tags | reverse).0
}

def parse-commits [conf, tag: string] {
    let commits = git log $"($tag)..HEAD" --pretty=%s | lines | reverse

    # parse conventional commits
    let conventional_commits = $commits | parse --regex '^(?P<type>\w+)(?:\((?P<scope>\w+)\))?(?P<breaking>!)?:\s*+(?P<description>.*+)$'
    let conventional_commits = $conventional_commits | insert is_conventional true

    if ($conventional_commits | length) != ($commits | length) {
        let err_msg = "Non-conventional commits were found in the git history."
        if not $conf.force and not $conf.interactive {
            log error $"($err_msg) ($FORCE_PROMPT)"
            exit 1
        }
        log warning $err_msg
        if $conf.interactive and not $conf.force {
            if not (yes-or-no $"($err_msg) They will be ignored. Continue anyway?") {
                exit 0
            }
        }
    }

    # combine and fill missing data
        print $commit
    let commits = $conventional_commits
    let commits = ($commits | update breaking {|it| $it.breaking == "!"})
    return $commits
}

def repr-semver [semver: record] -> string {
    return $"v($semver.major).($semver.minor).($semver.patch)"
}

def calculate-semvers [conf, latest_tag: string, commits] {
    # fully parse semver string
    mut latest_tag = ($latest_tag | parse --regex '^v(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<pre_release>[^+]*+))?(?:\+(?P<metadata>.*+))?$').0
    $latest_tag.major = ($latest_tag.major | into int)
    $latest_tag.minor = ($latest_tag.minor | into int)
    $latest_tag.patch = ($latest_tag.patch | into int)
    if $latest_tag.pre_release != "" {
        log warning "pre-release tags are not currently supported and are ignored."
        exit 1
    }
    $latest_tag.pre_release = ""
    $latest_tag.metadata = ""

    # detect type of update, note that updates containing "unknown" will never be patches
    let types = $commits | group-by type | columns
    if not ($types | all {|it| $it in $COMMIT_TYPES}) {
        let err_msg = $"Commit type not supported. Supported types are: ($COMMIT_TYPES | str join ', ')"
        if not $conf.force and not $conf.interactive {
            log error $"($err_msg) ($FORCE_PROMPT)"
            exit 1
        }
        log warning $err_msg
        if $conf.interactive and not $conf.force {
            if not (yes-or-no $"($err_msg) They will be ignored. Continue anyway?") {
                exit 0
            }
        }
    }
    let is_major = $commits | get breaking | any {|it| $it}
    let is_minor = not $is_major and ($types | any {|it| $it in ["__unknown", "feat"]})
    let is_patch = not $is_major and not $is_minor

    # generate version proposals
    let next_major = repr-semver ($latest_tag | update major {|it| $it.major + 1} | update minor 0 | update patch 0)
    let next_minor = repr-semver ($latest_tag | update minor {|it| $it.minor + 1} | update patch 0)
    let next_patch = repr-semver ($latest_tag | update patch {|it| $it.patch + 1})

    mut proposed_versions = {major: $next_major, minor: $next_minor, patch: $next_patch}
    if $is_major {
        $proposed_versions.suggested = $proposed_versions.major
    }
    if $is_minor {
        $proposed_versions.suggested = $proposed_versions.minor
    }
    if $is_patch {
        $proposed_versions.suggested = $proposed_versions.patch
    }
    return $proposed_versions
}

def prepare-changelog [conf, commits, next_version: string] -> string {
    let path = $"($conf.repo_root)/CHANGELOG.md"
    mut changelog = {
        changes: [],
        fixes: [],
        internal_changes: [],
    }
    if not ($path | path exists) {
        $changelog.previous_changes = ""
    } else {
        $changelog.previous_changes = (open $path --raw | decode utf-8 | parse --regex '<!-- ?changes ?-->(.*)<!-- ?/ ?changes ?-->').0
    }
    for commit in $commits {
        match $commit.type {
            "chore" => { $changelog.internal_changes = $changelog.internal_changes ++ {title: $commit.description, breaking: $commit.breaking} }
        }
    }

}

# run lints and prepare next release (without publishing)
def "main release" [
    --force (-f),                   # ignore most errors and continue with sensible defaults
    --no-interactive (-I),          # do not query for inputs, assume sensible defaults
    --skip-pre-flight-checks (-s)   # skip pre-flight checks (this automatically enables --dry)
    --dry (-d)                      # do not apply any changes, prints a summary of the changes that would be made
] {
    let conf = {
        repo_root: (git rev-parse --show-toplevel)
        force: $force
        interactive: (not $no_interactive)
    }

    if not $skip_pre_flight_checks {
        pre-flight-checks $conf
    } else {
        log info "skipping pre-flight checks"
    }
    let latest_tag = find-latest-tag
    let commits = parse-commits $conf $latest_tag
    let proposed_versions = calculate-semvers $conf $latest_tag $commits
    mut next_version = ""
    if $conf.interactive {
        let prompt = "Choose the next version"
        $next_version = ([
            [name value];
            [$"default \(($proposed_versions.suggested)\)" $proposed_versions.suggested]
            [$"major \(($proposed_versions.major)\)" $proposed_versions.major]
            [$"minor \(($proposed_versions.minor)\)" $proposed_versions.minor]
            [$"patch \(($proposed_versions.patch)\)" $proposed_versions.patch]
            ["custom", "custom"]
        ] | input list -d name $prompt).value
        if $next_version == custom {
            $next_version = (input "Custom version (please respect SemVer!, e.g. v1.2.3): ")
        }
    } else {
        $next_version = $proposed_versions.suggested
    }
    let changelog = prepare-changelog $conf $commits $next_version
    ed $"($conf.repo_root)/CHANGELOG.md"
    print "lmao"
}



def main [] {}
