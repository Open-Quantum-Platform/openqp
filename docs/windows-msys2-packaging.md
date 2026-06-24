# Windows MSYS2 Packaging

This track is separate from the Linux and macOS PyPI wheel release work. Its
goal is a native Windows command-line OpenQP package built in MSYS2, without GUI
or GPU packaging scope.

## Target

Build and validate OpenQP in the MSYS2 UCRT64 environment:

- `openqp` command-line launcher
- package-local OpenQP runtime files
- `liboqp.dll` plus the required MSYS2 runtime DLLs
- bundled basis sets and examples
- no CUDA, GPU, or desktop GUI dependency

The first distributable is a local MSYS2 pacman package plus a small clickable
installer wrapper. The package recipe can later be hardened for official MSYS2
repository submission once the remaining dependency gap is closed.

## Why UCRT64

MSYS2 provides several environments. UCRT64 is the best default for new native
Windows builds because it uses the Universal C Runtime and keeps the MinGW-w64
toolchain, Python, CMake, Ninja, and Fortran compiler in one coherent prefix.
Avoid mixing this with Python from python.org or Visual Studio compilers in the
first pass.

## Package Build Shape

Use the MSYS2 UCRT64 Python and compilers for everything. The package recipe is
`tools/windows_msys2/mingw-w64-openqp/PKGBUILD`; build it with:

```bash
pacman -Syu
git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
OQP_MSYS2_INSTALL_DEPS=1 tools/windows_msys2/build_pkg.sh
```

This creates a normal pacman package under `dist/msys2-ucrt64/packages`, for
example `mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst`.

Install that package in a UCRT64 shell with:

```bash
pacman -U dist/msys2-ucrt64/packages/mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst
python -m pip install --user basis_set_exchange geometric  # temporary dependency bridge
openqp examples/other/h2o_rhf_6-31g_hf.inp --omp 2
```

For a clickable preview installer, use the MSI produced by the workflow:

- `OpenQP-MSYS2-UCRT64-<version>.msi`

The MSI installs an OpenQP installer bundle under Program Files and creates an
`Install OpenQP for MSYS2` shortcut. That shortcut runs the bundled PowerShell
installer, detects MSYS2 at `C:\msys64`, installs the embedded `.pkg.tar.*`
artifact with `pacman -U`, installs the temporary PyPI dependency bridge for
`basis_set_exchange` and geomeTRIC, checks that `oqp` imports, and creates a
desktop `OpenQP UCRT64.cmd` launcher.

This MSI is still a preview artifact. CI builds the MSYS2 package, runs an
OpenQP smoke calculation, builds the MSI, silently installs the MSI wrapper, and
checks that the wrapper payload is present under Program Files. It does not yet
prove the full clean-machine end-user flow where a user installs MSYS2, runs the
MSI, clicks the Start Menu shortcut, and completes the embedded pacman install.

Uninstalling the MSI removes only the wrapper bundle and shortcuts. It does not
remove the package already installed inside MSYS2; remove that separately from a
UCRT64 shell with `pacman -R mingw-w64-ucrt-x86_64-openqp`.

For a no-MSI fallback, place these three files in the same Windows folder:

- `install_openqp_msys2.cmd`
- `install_openqp_msys2.ps1`
- `mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst`

Then double-click `install_openqp_msys2.cmd`. If MSYS2 is installed elsewhere,
run:

```powershell
.\install_openqp_msys2.cmd -Msys2Root D:\msys64
```

The script configures OpenQP with:

- `USE_LIBINT=OFF`
- `ENABLE_OPENMP=ON`
- `ENABLE_OPENTRAH=OFF`
- `LINALG_LIB=netlib`
- `LINALG_LIB_INT64=ON`

`LINALG_LIB=netlib` is intentional for the first Windows proof. It avoids
starting with OpenBLAS DLL and ILP64-discovery complexity. Once the native
package imports, runs a small job, and carries the right DLLs, OpenBLAS can be
added as a performance follow-up.

## CI Shape

Use a GitHub Actions workflow with `msys2/setup-msys2`:

- runner: `windows-2022` or newer
- MSYS2 environment: `UCRT64`
- packages: UCRT64 GCC/GFortran, Python, CMake, Ninja, and zip tools
- artifact: pacman `.pkg.tar.*` package plus the clickable installer files

The workflow supports manual preview runs with `workflow_dispatch` and automatic
release runs with `release: published`. It remains separate from the Linux/macOS
PyPI wheel matrix and does not upload a Windows wheel to PyPI.

The release workflow validates the OpenQP pacman package and the MSI wrapper
installation. Treat Windows MSI releases as preview until the full installer has
been exercised on a fresh Windows machine with a normal MSYS2 install.

## Installability Status

The supported package path is now a local MSYS2 pacman package. Build with
`build_pkg.sh`, distribute the produced `.pkg.tar.zst`, and install with
`pacman -U` from an MSYS2 UCRT64 shell.

There is one important upstream-packaging blocker: OpenQP imports
`basis_set_exchange` at runtime and supports geomeTRIC workflows, but those two
Python projects are not yet available as MSYS2 UCRT64 packages. The local package
builder can install them from PyPI for validation; an official MSYS2 submission
should first add package recipes for them or make the OpenQP imports fully
optional.

## Distribution

For preview testing, run the `Windows MSYS2 package` workflow manually. For
normal releases, publish the GitHub Release from the OpenQP tag. The same
workflow runs on `release: published`, rebuilds the MSYS2 package and MSI from
that tag, creates `SHA256SUMS.txt`, and uploads these files back to the release:

- `OpenQP-MSYS2-UCRT64-<version>.msi`
- `mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst`
- `install_openqp_msys2.cmd`
- `install_openqp_msys2.ps1`
- `SHA256SUMS.txt`

Users install MSYS2 once, run the MSI, then launch `Install OpenQP for MSYS2`
from the Start Menu or Desktop. The MSI still uses the real pacman package
underneath, so the same artifact can be tested manually with `pacman -U` from an
MSYS2 UCRT64 shell. Do not label it as an official MSYS2 repository package until
the missing Python dependencies are packaged for MSYS2 or made optional in
OpenQP.

## Validation Gates

Before calling the MSYS2 package useful:

- `python -c "import oqp; print(oqp.oqp_root)"` succeeds from the installed
  prefix, not from the source tree.
- the `openqp` launcher from the installed prefix completes
  `examples/other/h2o_rhf_6-31g_hf.inp --omp 2`.
- The installed prefix contains `oqp/include/oqp.h`, `oqp/lib/liboqp.dll`, and
  `oqp/share/basis_sets`.
- A fresh UCRT64 shell can run `openqp` from the pacman-installed package without
  pointing `OPENQP_ROOT` at the source tree.
- The produced MSI installs silently and places `install_openqp_msys2.cmd`,
  `install_openqp_msys2.ps1`, the embedded `.pkg.tar.zst`, and `README.txt`
  under Program Files.
- A human clean-machine test completes the MSI shortcut flow against a normal
  MSYS2 UCRT64 installation.
- A later official submission has MSYS2 packages for `basis_set_exchange` and
  geomeTRIC, or OpenQP treats them as optional runtime features.

## Non-goals

- No GUI application, GPU helper, NSIS installer, or MSIX installer in this
  pass. The MSI is only a Windows-installer wrapper around the MSYS2 pacman
  package and PowerShell installer.
- No GPU or CUDA support in this pass.
- No PyPI Windows wheel claim until MSYS2-built runtime behavior is proven and
  the Python ABI/toolchain story is decided.
- No changes to the Linux/macOS release wheel PR.
