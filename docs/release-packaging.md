# OpenQP Release Packaging

This repository owns the build and publication of the `openqp` Python package.
The user-facing install command is:

```bash
pip install openqp
```

The `.github/workflows/build_wheels.yml` workflow builds source distributions
and one Linux smoke wheel whenever release-sensitive files change in an ordinary
pull request. The smoke wheel compiles OpenQP source against a restored bundled
externals cache so source changes are checked without rebuilding third-party
libraries on every PR. Full Linux and macOS wheel builds run for pull requests
labeled `release`, manual workflow dispatch, and GitHub Releases. Only the
`release: published` event uploads to PyPI.

The workflow does not require `OQP_EXTERNALS_ROOT`. Linux wheel jobs set the
standard `XDG_CACHE_HOME` to a cached checkout-local directory so OpenQP's
existing cache auto-discovery uses `$XDG_CACHE_HOME/openqp/externals`. macOS
wheel jobs use the default `~/Library/Caches/openqp/externals` location.

## Release Checklist

1. Update `project.version` in `pyproject.toml`.
2. Synchronize the version in `CMakeLists.txt` and `pyoqp/setup.py`, then run
   `python -m unittest tests.test_release_metadata`.
3. Write detailed release notes as described below and include them in the
   release-preparation pull request for review.
4. Merge the release-preparation pull request only after its checks pass.
5. Create and push a matching tag, for example `v1.2.0`.
6. Publish a GitHub Release from that tag using the reviewed release notes.
7. Wait for the full wheel workflow to finish, then confirm that the GitHub
   Release assets and PyPI files include the expected wheels and source
   distribution.

The workflow verifies that the GitHub Release tag is exactly `v` plus the
`pyproject.toml` version. For example, `version = "1.2.0"` must be released as
`v1.2.0`.

Pushing a `v*` tag by itself does not start this workflow. Publish a GitHub
Release from the tag instead, so the full wheel build and PyPI upload happen
once from the release event.

## Release Notes Requirements

Every release must describe the complete change set since the previous release,
not just the version bump or a short list of highlights.

1. Generate GitHub's candidate notes between the previous and new tags to
   collect merged pull requests and contributors.
2. Audit `git log <previous-tag>..HEAD` as well, because direct commits and
   multi-commit pull requests may not appear as separate entries in generated
   notes.
3. Reconcile the two lists so every merged pull request and every user-visible
   direct commit is represented exactly once.
4. Organize the notes under these headings where applicable:
   - Highlights
   - New features and scientific capabilities
   - Performance improvements
   - Bug fixes and reliability
   - Packaging, platforms, and CI
   - Documentation and developer experience
5. Explain the user or developer impact of each item and link it to its pull
   request or commit. Call out changed defaults, compatibility constraints,
   optional dependencies, migration steps, and known limitations explicitly.
6. End with contributor acknowledgements and a full comparison link between
   the two release tags.

Use the following structure for the reviewed release-note document and the
GitHub Release body:

```markdown
# OpenQP vX.Y.Z

## Highlights

## New features and scientific capabilities

## Performance improvements

## Bug fixes and reliability

## Packaging, platforms, and CI

## Documentation and developer experience

## Contributors

**Full changelog:** https://github.com/Open-Quantum-Platform/openqp/compare/vPREVIOUS...vX.Y.Z
```

## PyPI Trusted Publishing

Configure the `openqp` project on PyPI with a GitHub Trusted Publisher instead
of storing a PyPI token in repository secrets.

Use these values:

- Owner: `Open-Quantum-Platform`
- Repository: `openqp`
- Workflow: `.github/workflows/build_wheels.yml`
- Environment: `pypi`

The workflow requests GitHub OIDC credentials only in the publishing job, and
that job runs only after a GitHub Release is published.

## Current Platform Plan

The first automated release workflow builds:

- Ordinary pull requests: source distribution and one Linux x86_64 CPython 3.11
  smoke wheel using the reusable bundled-externals cache
- Pull requests labeled `release`: full Linux and macOS wheel matrix
- Linux x86_64 manylinux wheels
- Linux aarch64 manylinux wheels
- macOS x86_64 wheels for macOS 15 or newer
- macOS arm64 wheels for macOS 15 or newer
- Source distribution

Windows wheels should be added after native Windows runtime layout is validated.
Until then, Windows users should use WSL2, Docker, or a future conda-forge
package.

## Installer Site Repository

A separate `Open-Quantum-Platform/openqp-install` repository provides the public
install page without owning the package build itself. It should link to:

- PyPI: `https://pypi.org/project/openqp/`
- GitHub Releases from this repository
- Docker image instructions
- Windows WSL2/Docker instructions until native Windows wheels are available

Keep source builds, wheel builds, PyPI publishing, and release assets in this
main `openqp` repository so every package is built from the exact source tag.
