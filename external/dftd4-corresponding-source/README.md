# DFT-D4 stack corresponding source

This directory is installed with OpenQP binary packages and container images.
It contains the complete source trees used to build the replaceable shared
libraries in OpenQP's `lib` directory.

| Installed source tree | Upstream archive | SHA-256 | License |
| --- | --- | --- | --- |
| `mctc-lib-0.4.2/` | <https://github.com/grimme-lab/mctc-lib/archive/refs/tags/v0.4.2.tar.gz> | `c7aa45c0a3e6f96e3316e15fc6cdbe48b15234940d3773927a57bb7bfe9744ac` | Apache-2.0 |
| `multicharge-0.3.0/` | <https://github.com/grimme-lab/multicharge/archive/refs/tags/v0.3.0.tar.gz> | `2fcc1f80871f404f005e9db458ffaec95bb28a19516a0245278cd3175b63a6b2` | Apache-2.0 |
| `dftd4-3.7.0/` | <https://github.com/dftd4/dftd4/archive/refs/tags/v3.7.0.tar.gz> | `f00b244759eff2c4f54b80a40673440ce951b6ddfa5eee1f46124297e056f69c` | LGPL-3.0-or-later |

The license documents and copyright notices supplied by each upstream project
are retained inside its source tree.

## OpenQP build-time modifications

The installed trees are the post-patch sources that produced the binaries:

- `patches/mctc-lib-0.4.2-disable-tests.patch` comments out mctc-lib's
  top-level `test` subdirectory.
- `patches/dftd4-3.7.0-disable-tests-and-mstore.patch` comments out DFT-D4's
  top-level `test` subdirectory and the test-only, unpinned `mstore` lookup.
- OpenQP makes no source-tree modification to multicharge 0.3.0.

These changes affect only which tests and test dependency are configured; they
do not change library implementation source. The two installed files are the
exact unified patches applied by the OpenQP ExternalProjects with `patch -p1`;
their SHA-256 values are recorded in `BUILD-INFO.json`. The installed
`apply-patch.cmake` helper documents the idempotent clean/already-applied/error
state checks used for reusable source caches.

`BUILD-INFO.json` is generated during OpenQP configuration and records the
resolved build: archive pins and licenses, CMake and generator versions,
platform, compiler IDs/versions and executable basenames, forwarded C and
Fortran flags, OpenMP, resolved BLAS/LAPACK library names and integer width,
shared-library mode, and canonical runtime filenames. An exact Git commit is
included only when the build source has verifiable Git metadata and a clean
working tree; otherwise the field is JSON `null`. Transient source, build, and
cache path prefixes are replaced with descriptive placeholders, and resolved
libraries and compiler executables are recorded by name rather than by volatile
absolute path.

`openqp-external-build.cmake` is a reference copy of the OpenQP ExternalProject
recipe. It contains variable-based selection logic and the archive/patch pins,
but it is neither a standalone build script nor a record of the values resolved
for a particular binary. Use `BUILD-INFO.json` for those values and the complete
post-patch source trees plus `patches/` to reproduce or modify the libraries.

The three libraries are built with `BUILD_SHARED_LIBS=ON`. DFT-D4 and
multicharge use their `WITH_ILP64` option when OpenQP uses an eight-byte BLAS
interface; their ordinary Fortran integers remain four bytes. `WITH_OpenMP`
matches the OpenQP build. On macOS Accelerate ILP64 builds, the two BLAS-using
shared libraries and their command-line tools receive OpenQP's generated
alias-list link option so classic Fortran names resolve to the genuine
`$NEWLAPACK$ILP64` entry points rather than Accelerate's LP64 compatibility
interface.

To use a modified implementation, rebuild the ABI-compatible shared library
with the same Fortran compiler ABI and BLAS/OpenMP settings, then replace its
SOVERSION-named file in OpenQP's package-local `lib` directory before starting
a new process.
