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
labeled `release`, manual workflow dispatch, `v*` tag pushes, and GitHub
Releases. Only the `release: published` event uploads to PyPI.

## Release Checklist

1. Update `project.version` in `pyproject.toml`.
2. Create and push a matching tag, for example `v1.2.0`.
3. Publish a GitHub Release from that tag.
4. Wait for the full wheel workflow to finish.
5. Confirm that the release assets and PyPI files include the expected wheels.

The workflow verifies that the GitHub Release tag is exactly `v` plus the
`pyproject.toml` version. For example, `version = "1.2.0"` must be released as
`v1.2.0`.

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
