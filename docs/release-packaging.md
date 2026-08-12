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
labeled `release`, manual build dispatches, and published GitHub Releases.

A published GitHub Release builds and verifies candidate distributions but does
not make the wheel files public. It neither attaches them to the release nor
uploads them to PyPI. Publication is a separate, explicit manual action protected
by exact tag, commit, workflow-run, package-metadata, and GitHub-environment
checks.

Neither that two-stage workflow nor a source-only GitHub Release resolves
copyright ownership. Before a v1.3.0 branch is merged, tagged, or made public
under the OpenQP Research License, the release owner must clear the separate
source-distribution authorization gate below.

The workflow does not require `OQP_EXTERNALS_ROOT`. Linux wheel jobs set the
standard `XDG_CACHE_HOME` to a cached checkout-local directory so OpenQP's
existing cache auto-discovery uses `$XDG_CACHE_HOME/openqp/externals`. macOS
wheel jobs use the default `~/Library/Caches/openqp/externals` location.

## Release Checklist

1. Update `project.version` in `pyproject.toml`.
2. Synchronize the version in `CMakeLists.txt` and `pyoqp/setup.py`, then run
   `python -m unittest tests.test_release_metadata`.
3. Verify that `project.license`, `project.license-files`, `LICENSE`,
   `LICENSING.md`, `THIRD_PARTY_NOTICES.md`, and every referenced third-party
   license describe the exact source and binaries being built.
4. Write detailed release notes as described below and include them in the
   release-preparation pull request for review.
5. Complete the source-distribution authorization gate below. If it is not
   cleared, keep the release branch private and do not merge, tag, or publish a
   GitHub Release under the new license.
6. Merge the release-preparation pull request only after its checks and the
   authorization gate pass.
7. Create and push a matching tag, for example `v1.3.0`.
8. Publish a GitHub Release from that tag using the reviewed release notes. It
   is a source release at this point; do not manually attach unverified binaries.
9. Wait for the full distribution workflow to finish and record the successful
   release-event workflow run ID and exact 40-character tag commit SHA.
10. Complete the binary-distribution gate below. If any item is unresolved, stop
   here and leave the release source-only.
11. In the Actions UI select the exact release tag (for example `v1.3.0`) in
    **Use workflow from**, then manually dispatch `build_wheels.yml` with action
    `publish_pypi`, the same exact tag, SHA, verified release-run ID, and the
    requested confirmation text. A dispatch from `main` is rejected by the
    `pypi` environment's `v*` tag policy.
12. Review the protected `pypi` deployment, approve it only after checking the
    displayed tag and SHA, and then verify both GitHub release assets and PyPI.

The workflow verifies that the GitHub Release tag is exactly `v` plus the
`pyproject.toml` version. For example, `version = "1.2.0"` must be released as
`v1.2.0`.

Pushing a `v*` tag by itself does not start the distribution workflow. Publishing
a GitHub Release from the tag builds and verifies the full candidate matrix.
It still does not authorize public binary distribution.

## Source Distribution Authorization Gate

Before any v1.3.0 source is published under the Research License, obtain and
archive evidence that Open Quantum Inc. has sufficient authority to license
every surviving copyrightable contribution on both the research and commercial
terms. This includes founder rights, historical contributors, employer or
institution interests, university and sponsor rules, and imported third-party
material. A Git commit, pull-request checkbox, DCO, or newly adopted CLA is not
by itself retroactive permission for earlier contributions.

The signed founder/company assignment, contributor schedule and permissions,
university/funder review, and third-party provenance record must be reviewed by
qualified counsel. If any required authority is missing, do not publish the
new-license branch, tag, source archive, wheel, or container; either complete
the permissions, retain the applicable historical GPL route, or replace the
affected material with independently documented code.

## Binary Distribution Gate

Before approving either wheel, Docker-image, or other binary publication:

1. Confirm that every bundled LGPL library remains a separately replaceable
   shared library, including the DFT-D4 stack, and inspect the installed wheel or
   image rather than relying only on CMake options.
2. Confirm from both the dependency table and dynamic symbols that NLopt is
   absent and the native deterministic simplex-QP implementation is in use.
3. Confirm that license texts, copyright notices, modification notices, source
   locations, and any required written offer or corresponding source are
   shipped with the exact artifact.
4. Run platform-specific dependency inspection (`auditwheel`, `delocate`, and
   runtime load tests) on every wheel family.
5. Confirm that the company owns or has sufficient written permission for all
   OpenQP contributions being offered under the research/commercial terms.

`THIRD_PARTY_NOTICES.md` is necessary but is not, by itself, proof that a static
LGPL combination satisfies the relinking requirements. Never approve publication
solely because the notice file is present.

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

The workflow requests GitHub OIDC credentials only in the publishing job. That
job can run only from an explicit `publish_pypi` manual dispatch after preflight
checks and approval in the protected `pypi` environment. The environment allows
only `v*` tags, requires the configured reviewer, and does not permit admin
bypass.

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

Candidate files are retained as GitHub Actions artifacts until the explicit
publication gate is approved. Building an artifact is not authorization to
publish it.

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
Docker builds follow the same rule: pull-request and release validation may
build locally in CI, but the workflow must not push an image automatically.

## Docker Distribution Policy

Docker remains a supported OpenQP installation route. The public research image
and a commercially licensed deployment should normally use the same reviewed,
versioned bits; the OpenQP license grant, not a scientifically divergent binary,
determines the permitted use. Pulling a public image does not authorize use by
or for a commercial entity without a written commercial license.

Publish only an immutable version tag built from the exact matching release tag
and commit. The publication job must be a manual dispatch from that tag, pass
the same binary-distribution gate as wheels, and wait for approval in a
tag-restricted GitHub environment whose registry credentials are environment
secrets. Do not publish `latest` until the versioned digest has been verified.

The published image must:

1. embed `LICENSE`, `LICENSING.md`, `THIRD_PARTY_NOTICES.md`, and all applicable
   third-party license texts at a documented location;
2. retain or accompany the exact corresponding source and build changes required
   by every redistributed copyleft component;
3. demonstrate that the DFT-D4 stack is dynamically linked and replaceable and
   that NLopt is absent;
4. carry OCI source, version, revision, and license labels; and
5. publish an SBOM and build-provenance attestation together with the immutable
   image digest.

The current `docker-build.yml` is deliberately build-only while final compliance
and release authorization remain pending. Re-enable registry publication only
through the manual protected path described above.

### Current container candidate implementation

The candidate Dockerfile uses two digest-pinned inputs recorded in
`docker/base-images.lock.json`. The `openqp-buildenv:1` image is a builder stage
only; OpenQP is built as a wheel, installed from an offline wheelhouse, and
exercised with the full numerical/ABI/DFT-D4 wheel smoke test. The final stage
starts from a pinned Python slim runtime and receives only the installed wheel,
its required shared-library closure, legal materials, and runtime manifests.
It does not inherit the compiler, source checkout, external cache, or static
archives from the build environment.

The current pinned builder is Linux/amd64 with Python 3.12, so the candidate
runtime is also Python 3.12. A Python 3.11 container requires a new matching
builder; do not describe the current image as Python 3.11-compatible exact
bits. The proposed companion-repository migration is documented in
`docs/openqp-dev-buildenv-hardening-plan.md`.

The build produces an OCI archive with an embedded SPDX SBOM and max-mode
provenance, verifies the OCI layout and authenticated descriptors, validates
the required SPDX/SLSA predicate bodies and image subjects, and then lets the
ephemeral runner discard it. The workflow has no registry
credentials, package-write permission, login, push, release attachment,
retained image artifact, or exported build cache.

The container remains publication-ineligible for two explicit reasons:

1. build-system and runtime dependency wheels are resolved during a candidate
   build and hashed only afterward; publication requires a reviewed hash lock
   or reuse of already verified exact release-run wheelhouses; and
2. copied Ubuntu GNU runtime libraries have exact binary/source package
   versions and copyright files recorded, but the matching source archives are
   not yet bundled for components whose terms require corresponding source.

Both `/usr/share/licenses/openqp/wheelhouse-manifest.json` and
`runtime-library-manifest.json` record `publication_gate.ready: false`. The
wheelhouse manifest separately inventories the exact build-system wheels used
to create OpenQP and the runtime wheels available to the final installation.
final-stage smoke test also rejects NLopt files/dependencies/symbol strings,
build tools, caches, static archives, alternate DFT-D4 copies, missing
corresponding source, unresolved dependencies, and legal-file omissions.
