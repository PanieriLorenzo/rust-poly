# Developer Guide


## DevOps

We don't do GitHub actions over here, because I'm a solo developer and actions are expensive.

Automation is performed using nushell scripts found in `devops/` and with pre-commit hooks.

Please install pre-commit hooks before making any commits by running:
```bash
npm ci
```

### Conventional Commits

Commit messages must adhere to the [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) specification and additionally to the [commitlint conventional config](https://github.com/conventional-changelog/commitlint/tree/master/%40commitlint/config-conventional).

**Breaking Changes**:
- Any other type that ends in a bang, like `feat!` or `fix!`. Note that `feat!` can also be used to mark removal of features.

**Minor Changes**:
- `feat`: any new functionality that is visible to the user

**Patch Changes**:
- `fix`: a bugfix
- `docs`: documentation changes
- `build`: changes to build scripts
- `ci`: changes to CI pipelines
- `test`: adding tests or benchmarks
- `perf`: changes that affect performance
- `refactor`: major refactoring
- `revert`: reverting a change
- `style`: stylistic changes
- `chore`: any changes that are mostly administrative, e.g. small refactors, code style, comments, semver adjustments, etc...

### Changelog

The changelog is generated automatically from commit messages. During the publishing process of a new release, the generated changelog can be manually edited to include additional information or rephrase the changes.

Right now, changelogs are generated with [git-cliff](https://github.com/orhun/git-cliff). Configurations for how this happens are in `cliff.toml`.

### Releasing a New Version

```bash
git cliff --tag v0.3.0 --unreleased --prepend CHANGELOG.md
```

### Benchmarking

```bash
cargo flamegraph --bench bench -- <name-of-bench> --bench
```
