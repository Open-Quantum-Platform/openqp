# Windows MSYS2 Packaging

This directory contains the first command-line packaging path for Windows using
MSYS2 UCRT64.

## Build The Pacman Package

Install MSYS2 from `https://www.msys2.org/`, open the **UCRT64** shell, update
MSYS2, then build the package:

```bash
pacman -Syu
# close/reopen UCRT64 if MSYS2 asks, then run pacman -Syu again

git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
OQP_MSYS2_INSTALL_DEPS=1 tools/windows_msys2/build_pkg.sh
```

This creates and installs a package like:

```bash
dist/msys2-ucrt64/packages/mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst
```

Users can then install that artifact in UCRT64 with:

```bash
pacman -U mingw-w64-ucrt-x86_64-openqp-<version>-1-any.pkg.tar.zst
python -m pip install --user basis_set_exchange geometric  # temporary dependency bridge
openqp examples/other/h2o_rhf_6-31g_hf.inp --omp 2
```

For the friendliest install path, use the MSI produced by the manual
`Windows MSYS2 package` workflow:

```text
OpenQP-MSYS2-UCRT64-<version>.msi
```

The MSI installs the OpenQP MSYS2 installer bundle and creates an
`Install OpenQP for MSYS2` shortcut. The shortcut installs the embedded pacman
package into an existing MSYS2 UCRT64 installation.

This MSI is a preview artifact. The workflow checks that the MSI wrapper builds,
installs silently, and contains the expected files, but the full clean-machine
click-through install still needs real Windows/MSYS2 validation before calling it
stable.

Uninstalling the MSI removes the wrapper bundle and shortcuts, not the OpenQP
package installed inside MSYS2. Remove that separately from UCRT64 with
`pacman -R mingw-w64-ucrt-x86_64-openqp`.

For a simpler no-MSI double-click install, put the package artifact next to
`install_openqp_msys2.cmd` and `install_openqp_msys2.ps1`, then double-click:

```text
install_openqp_msys2.cmd
```

The clickable installer expects MSYS2 UCRT64 under `C:\msys64` by default. If
MSYS2 is somewhere else, pass the root path:

```powershell
.\install_openqp_msys2.cmd -Msys2Root D:\msys64
```

It installs the pacman package, installs the temporary `basis_set_exchange` and
geomeTRIC Python dependency bridge, checks the OpenQP import, and creates a
desktop `OpenQP UCRT64.cmd` launcher.

To build the MSI locally on Windows after creating the pacman package, install
the WiX Toolset command-line tool:

```powershell
dotnet tool install --global wix
tools/windows_msys2/build_msi.ps1
```

## Install From Source Prefix

Install MSYS2 from `https://www.msys2.org/`, open the **UCRT64** shell, then run:

```bash
pacman -Syu
```

If MSYS2 asks you to close the terminal, close it, reopen **UCRT64**, and run the
same update command again. Then clone OpenQP and build:

```bash
git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
OQP_MSYS2_INSTALL_DEPS=1 tools/windows_msys2/build_msys2.sh
```

For an already-prepared UCRT64 shell, this shorter command is enough:

```bash
tools/windows_msys2/build_msys2.sh
```

The script builds a Windows wheel, installs that wheel into
`dist/msys2-ucrt64/prefix`, imports `oqp` from that installed prefix, and runs:

```bash
openqp examples/other/h2o_rhf_6-31g_hf.inp --omp 2
```

The output is a wheelhouse and an installed-prefix archive under
`dist/msys2-ucrt64`. The script validates the installed prefix by importing
`oqp` from that prefix and running the `openqp` launcher from that prefix. This
is not a GUI package and does not include GPU support.

To use the installed prefix in a later UCRT64 shell:

```bash
source dist/msys2-ucrt64/prefix/activate-openqp.sh
openqp examples/other/h2o_rhf_6-31g_hf.inp --omp 2
```

## Status

`mingw-w64-openqp/PKGBUILD` is the real MSYS2 package recipe. The remaining
upstream-packaging blocker is dependency coverage: `basis_set_exchange` and
geomeTRIC still need MSYS2 Python packages before OpenQP can be submitted cleanly
to the official MSYS2 repository.

`build_msys2.sh` is kept as a lower-level source-prefix validation helper. The
zip archive it creates is not a standalone installer.

## Release Automation

The `.github/workflows/windows-msys2.yml` workflow can be run manually for
preview testing. On a published GitHub Release, it also rebuilds the MSYS2
package and MSI from that release tag, creates `SHA256SUMS.txt`, and uploads the
Windows assets back to the same release.
