# `openqp-dev` Build-Environment Hardening Plan

This is a patch plan for `Open-Quantum-Platform/openqp-dev`; it does not modify
or publish that repository.  It is based on the current one-commit `main`
layout: `docker/buildenv.Dockerfile`, `.github/workflows/publish-buildenv.yml`,
and `README.md`.

## Current risk to remove

The current Dockerfile starts from mutable `ubuntu:24.04`, downloads OpenBLAS
without a checksum, installs moving Ubuntu/PyPI dependency sets, and publishes
the mutable `openqp/openqp-buildenv:1` tag automatically after a matching push
to `main` whenever Docker Hub credentials exist.  The consumer cannot determine
the builder bits from the tag alone.

Registry values resolved on 2026-08-09 (record these in a lock file and
re-resolve deliberately when updating) are:

- `ubuntu:24.04` index:
  `sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea`
- `ubuntu:24.04` Linux/amd64 manifest:
  `sha256:019e8eb29a85e74d64925745884f2ec79aa27e3feab36353d24656f4d6b89467`
- proposed `python:3.11.15-slim-trixie` index:
  `sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff`
- proposed Python Linux/amd64 manifest:
  `sha256:78b39ef14d8e2b4d71f8dc304f1328c37df95fe0ef99477c2ae6bd3d03784553`

These values are review inputs, not permission to publish.

## Patch 1: make the existing builder reproducible and inspectable

1. Add `docker/base-images.lock.json` with the Ubuntu reference, index digest,
   Linux/amd64 manifest digest, resolution date, expected Python/GCC versions,
   OpenBLAS version, and upstream OpenBLAS archive SHA-256.
2. Change the first Dockerfile line to
   `FROM ubuntu:24.04@sha256:561618... AS buildenv`; CI must compare the full
   line to the lock file before building.
3. Replace unverified `wget` extraction with a download followed by
   `sha256sum --check` against the reviewed lock value before `tar` runs.
4. Point APT at one reviewed Ubuntu snapshot and install exact package versions
   from a checked-in `docker/ubuntu-packages.lock`.  Fail when any version is no
   longer available; do not silently fall back to the moving archive.
5. Remove `pip3 install dftd4`.  OpenQP now uses its native package-local DFT-D4
   stack, and the main project installs its own Python build requirements.
   Either remove `cffi` as well or install Python tooling only from a
   hash-locked `docker/python-build-requirements.lock` using
   `pip --require-hashes`.
6. Add `scripts/verify_buildenv.py` and run it in the Dockerfile.  It must check
   the exact Python, GCC/G++/GFortran, CMake, Ninja, and OpenBLAS versions;
   call `openblas_get_config`; require ILP64, unsuffixed BLAS symbols, OpenMP,
   and `DYNAMIC_ARCH`; reject LP64 or NetLib fallback; and emit a small JSON
   component manifest with hashes but no transient absolute paths.
7. Preserve the OpenBLAS BSD notice and record every copied GNU runtime binary,
   its binary/source package version, copyright file, and verified matching
   source archive.  Stage that complete runtime kit under a versioned path such
   as `/opt/openqp-buildenv-runtime/` for the main repository to copy without
   rediscovering cross-distribution dependencies.

## Patch 2: publish only immutable, reviewed buildenv artifacts

Replace `.github/workflows/publish-buildenv.yml` with two jobs:

1. Resolve every third-party action release line to a reviewed full commit SHA
   and use `owner/action@<40-hex-sha> # vN`; do not execute moving major tags.
   Re-resolve and review those SHAs deliberately during dependency updates.
2. `pull_request` and ordinary `push` run a credential-free, build-only
   Linux/amd64 OCI validation.  Generate SBOM and max-mode provenance, validate
   them on the ephemeral runner, and do not upload the OCI archive or export a
   source-bearing cache.
3. A separate `workflow_dispatch` publication accepts the exact 40-character
   source SHA, immutable requested tag, and explicit confirmation text.  It
   must verify that checkout `HEAD` equals the input SHA and run only in a
   protected `docker-buildenv` GitHub environment with required review, no
   admin bypass, and environment-scoped Docker Hub credentials.
4. Use an immutable tag containing the compatibility tuple and source commit,
   for example
   `py311-gcc14-openblas0.3.30-<12-character-commit>`.  Never let the main
   OpenQP repository consume `:1`, `:latest`, or another moving alias.
5. Only the protected job may run `docker/login-action` and `push: true`.
   Request only the minimum workflow permissions.  Publish SBOM and provenance
   with the image, record the returned registry digest, and immediately verify
   the tag/digest/platform/attestation subjects with
   `docker buildx imagetools inspect`.
6. Open a reviewed main-repository change that updates
   `docker/base-images.lock.json` to the returned builder digest.  The OpenQP
   Docker workflow must continue to reference `tag@digest`, so a moved tag
   cannot change a release candidate.

## Patch 3: move the Docker line to Python 3.11 without ABI ambiguity

The current `openqp-buildenv:1` is Ubuntu 24.04/Python 3.12.  Do not install a
wheel produced there into a Python 3.11 runtime and call it equivalent.

1. Build a new v2 builder from the digest-pinned
   `python:3.11.15-slim-trixie` image above.
2. From one pinned Debian snapshot, install exact GCC 14/G++ 14/GFortran 14,
   CMake, Ninja, make, and download tools.  Build the reviewed OpenBLAS version
   (preferably align it with the release-wheel value, currently 0.3.30) as
   shared ILP64/OpenMP/DYNAMIC_ARCH with a verified archive SHA-256.
3. Run the buildenv verifier and stage the complete shared runtime closure,
   notices, and matching source packages.  Publish under a new immutable tag;
   do not overwrite the existing builder digest.
4. In the main OpenQP repository, update the builder lock and use the exact same
   pinned Python 3.11 slim base for the final stage.  Rebuild the OpenQP wheel,
   run the wheel and final-container numerical/D4/NLopt gates, and compare the
   OCI SBOM against the reviewed component inventory.
5. Retire the temporary Python 3.12 container path only after the Python 3.11
   candidate passes all release gates.  Historical digests remain immutable.

## Acceptance criteria before any public image

- Both builder and runtime `FROM` inputs match checked-in registry digests and
  Linux/amd64 manifests.
- Every network download and Python wheel is version/hash locked before build.
- OpenBLAS version/configuration is read from the actual binary; ILP64 and
  OpenMP are proven.
- The final image has no compiler, build source, cache, static archive, NLopt,
  or alternate DFT-D4 copy.
- Every copied shared object has an exact hash, notice, source-package record,
  and corresponding source where its license requires it.
- DFT-D4 is dynamically linked through the canonical SOVERSION files,
  individually replaceable, numerically validated, and accompanied by the
  complete patched corresponding source.
- The exact release tag/SHA/image digest, OCI labels, SPDX SBOM, and provenance
  are verified before a protected manual registry push.
- OpenQP chain-of-title and the repository-wide binary-distribution checklist
  have been cleared independently; container hardening alone is not a release
  authorization.
