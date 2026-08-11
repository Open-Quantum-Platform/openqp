# Third-Party Software Notices

OpenQP includes or links third-party software. The OpenQP Research License and
any separate commercial agreement apply only to OpenQP material for which the
Licensor has authority to grant those rights. Third-party software remains
subject to its own license terms.

This notice describes the default binary wheels built by
`.github/workflows/build_wheels.yml`. It does not replace the upstream license
texts in `licenses/third_party/`, nor does it cover every optional build
configuration.

## Components in the default OpenQP wheel

| Component | Version | How used in the default wheel | License | Upstream source |
| --- | --- | --- | --- | --- |
| Libxc | 7.0.0 | Statically linked into `liboqp` | MPL-2.0 | <https://gitlab.com/libxc/libxc/-/tree/7.0.0> |
| tagarray | 1.0.0 | Statically linked into `liboqp` | MIT; Copyright (c) 2023 Igor S. Gerasimov | <https://github.com/Open-Quantum-Platform/tagarray/tree/v1.0.0> |
| libecpint | 1.0.7 | Statically linked into `liboqp` | MIT; Copyright (c) 2021 Robert A. Shaw | <https://github.com/robashaw/libecpint/tree/v1.0.7> |
| Faddeeva implementation bundled by libecpint | libecpint 1.0.7 copy | Statically linked into `liboqp` | MIT; Copyright (c) 2012 Massachusetts Institute of Technology | <http://ab-initio.mit.edu/Faddeeva> |
| mctc-lib | 0.4.2 | Statically linked into `liboqp` as part of the DFT-D4 stack | Apache-2.0; Copyright 2020-2025 Sebastian Ehlert | <https://github.com/grimme-lab/mctc-lib/tree/v0.4.2> |
| multicharge | 0.3.0 | Statically linked into `liboqp` as part of the DFT-D4 stack | Apache-2.0; Copyright 2021 Sebastian Ehlert | <https://github.com/grimme-lab/multicharge/tree/v0.3.0> |
| DFT-D4 | 3.7.0 | Statically linked into `liboqp` | LGPL-3.0-or-later; Copyright 2017-2021 E. Caldeweyher, S. Ehlert, S. Grimme | <https://github.com/dftd4/dftd4/tree/v3.7.0> |
| OpenBLAS | 0.3.30 | Shared library bundled into Linux wheels by `auditwheel`; macOS wheels use the system Accelerate framework instead | BSD-3-Clause | <https://github.com/OpenMathLib/OpenBLAS/tree/v0.3.30> |
| GCC runtime libraries | Release-builder version | Linux wheels bundle `libgfortran`, `libgomp`, and `libquadmath`; macOS wheels may additionally bundle `libgcc_s` and `libstdc++` | The licenses marked in the corresponding GCC sources, including GPLv3 with GCC Runtime Library Exception 3.1 and LGPL-2.1-or-later components | <https://gcc.gnu.org/> |

The exact source archive URLs and SHA-256 checksums for the CMake-managed
libraries above are recorded in `external/CMakeLists.txt`. OpenBLAS is pinned
to tag `v0.3.30` by the wheel workflow. Wheel repair tools can add platform
runtime libraries after the OpenQP build; therefore each release artifact
should also be inspected after repair.

## License texts

- Apache-2.0: `licenses/third_party/apache-2.0.txt`
- BSD-3-Clause for OpenBLAS: `licenses/third_party/openblas-bsd-3-clause.txt`
- GCC Runtime Library Exception 3.1: `licenses/third_party/gcc-runtime-library-exception-3.1.txt`
- GNU GPL version 3: `licenses/third_party/gpl-3.0.txt`
- GNU LGPL version 2.1: `licenses/third_party/lgpl-2.1.txt`
- GNU LGPL version 3: `licenses/third_party/lgpl-3.0.txt`
- Libxc MPL-2.0: `licenses/third_party/libxc-mpl-2.0.txt`
- libecpint MIT license: `licenses/third_party/libecpint-mit.txt`
- Faddeeva MIT notice: `licenses/third_party/faddeeva-mit.txt`
- tagarray MIT license: `licenses/third_party/tagarray-mit.txt`

The Apache-2.0 text is shared by mctc-lib and multicharge. The GPLv3 and
LGPLv3 texts together constitute the license documents shipped upstream with
DFT-D4.

## OpenQP build changes to upstream source trees

The default build applies the following build-time changes:

- `external/fix_tagarray_linelen.py` wraps an overlong Fortran source line in
  tagarray so it compiles with current GNU Fortran compilers.
- `external/CMakeLists.txt` patches mctc-lib and DFT-D4 build files to omit
  their test subdirectories and DFT-D4's test-only `mstore` lookup. These
  patches do not alter the DFT-D4, multicharge, or mctc-lib library
  implementation sources.

The patch logic is distributed with OpenQP and is the preferred form for
reproducing these build changes.

## Important LGPL distribution requirement

Including notices and license texts is necessary but is not, by itself,
sufficient for distribution of a statically linked LGPL Combined Work. The
default build currently statically links DFT-D4 into `liboqp`. A distributor
must additionally satisfy the applicable LGPL source, installation-information,
and relinking requirements, including providing the materials needed to relink
the application with a modified version of the LGPL library where required.

NLopt has been removed from the OpenQP source, build, link, and binary
distribution paths and is therefore no longer a component of the default wheel.

The current OpenQP wheel build does not generate such a relinking kit. A public
binary release must therefore either change those dependencies to a compliant
dynamic-link arrangement, omit/replace the LGPL code, or add and verify the
required relinking materials before publication.

## Optional configurations

The default PyPI wheel disables Libint, ddX, OpenTrustRegion, MPI, and the
OpenQP-DFTB staging hook. Enabling any optional component can add further
third-party license and redistribution requirements. A distributor of a
non-default build must audit that artifact separately and include the
corresponding notices and license texts.
