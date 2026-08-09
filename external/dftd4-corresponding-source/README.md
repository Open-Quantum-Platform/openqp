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

- In mctc-lib 0.4.2, OpenQP comments out the top-level `test` subdirectory.
- In DFT-D4 3.7.0, OpenQP comments out the top-level `test` subdirectory and
  the test-only, unpinned `mstore` package lookup.
- OpenQP makes no source-tree modification to multicharge 0.3.0.

These changes affect only which tests and test dependency are configured; they
do not change library implementation source. `openqp-external-build.cmake` is a
copy of the OpenQP build recipe containing the exact archive pins, patch
commands, compiler settings, shared-library mode, and BLAS/OpenMP selection.

The three libraries are built with `BUILD_SHARED_LIBS=ON`. DFT-D4 and
multicharge use their `WITH_ILP64` option when OpenQP uses an eight-byte BLAS
interface; their ordinary Fortran integers remain four bytes. `WITH_OpenMP`
matches the OpenQP build.

To use a modified implementation, rebuild the ABI-compatible shared library
with the same Fortran compiler ABI and BLAS/OpenMP settings, then replace its
SOVERSION-named file in OpenQP's package-local `lib` directory before starting
a new process.
