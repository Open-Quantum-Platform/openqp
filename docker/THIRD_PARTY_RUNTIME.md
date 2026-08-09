# OpenQP Container Runtime Notice

This file supplements `THIRD_PARTY_NOTICES.md`, which primarily describes the
default wheel matrix.  It documents the separately built Linux container
candidate.  Third-party terms remain controlling for third-party material.

## Exact candidate inventory

The container is compiled by the digest-pinned `openqp-buildenv:1` image
recorded in `docker/base-images.lock.json`.  That builder is currently
Linux/amd64, Python 3.12, and is expected to contain OpenBLAS 0.3.28.  The build
does not trust that expectation alone: it calls `openblas_get_config` on the
actual linked shared library, rejects any other version, and records the
library SHA-256 in `/usr/share/licenses/openqp/runtime-library-manifest.json`.
This Docker-specific 0.3.28 inventory is distinct from wheel jobs that may
build and repair against a different reviewed OpenBLAS version.

The final image is based on the separately digest-pinned Python 3.12 slim
runtime in the same lock file.  The Python minor version is constrained by the
current builder; this image does not claim to contain Python 3.11 artifacts.
Exact dependency-wheel filenames, versions, compatibility tags, sizes, and
SHA-256 hashes for a candidate build are recorded in
`/usr/share/licenses/openqp/wheelhouse-manifest.json`.  Copyright files for
copied system shared libraries are installed under
`/usr/share/licenses/openqp/runtime-packages/`.

The DFT-D4 shared-library stack and its patched complete corresponding source
are installed by the OpenQP wheel under
`oqp/share/corresponding-source/dftd4-stack`.  The final-stage test verifies
the exact DFT-D4 dependency graph and loaded paths, removes each canonical
library in turn to prove there is no hidden fallback copy, and rejects static
archives or NLopt artifacts.

## Candidate-only publication gate

This container workflow is deliberately build-only.  Its OCI archive, SBOM,
and provenance are verified on the ephemeral CI runner and are neither pushed
to a registry nor retained as a release artifact.

Two inputs are intentionally recorded but are not yet sufficient for public
publication:

1. Python dependency wheels are resolved during the candidate build.  Their
   exact hashes are recorded afterward, but a public release requires either a
   reviewed checked-in hash lock or reuse of the already verified exact
   release-run wheelhouse.
2. The Ubuntu binary/source package names and versions for copied GCC/OpenMP/
   Fortran runtime libraries are recorded, but matching source archives are not
   bundled.  Public publication requires retrieving, verifying, and retaining
   the exact corresponding source for every component whose terms require it.

`runtime-library-manifest.json` and `wheelhouse-manifest.json` therefore carry
an explicit `publication_gate.ready: false` marker.  A future publication
workflow must fail unless those blockers have been replaced by reviewed,
verified materials and the complete binary-distribution checklist has passed.
