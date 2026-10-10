# Standalone and integrated packages

An OQP Studio release contains three related distributions. They use the same
OpenQP input and output formats, but they serve different purposes.

| Distribution | Contains | Use it when |
| --- | --- | --- |
| OQP Studio standard package | Desktop interface and Studio service | OpenQP is already installed locally, or the engine will be installed separately |
| Standalone OpenQP engine | Command-line `openqp` program, libraries, basis data, and examples | Calculations will run without the Studio interface, or one engine installation will be shared with scripts |
| OQP Studio `with-engine` package | Desktop interface, Studio service, and a bundled OpenQP engine | A complete offline desktop installation is preferred |

The word **standalone** in an asset such as
`openqp-1.3.1-macos-arm64.zip` refers to the command-line OpenQP engine. The
standard OQP Studio installer is a separate, smaller download. The
`with-engine` installer combines Studio and an engine in one package.

## Which package should I choose?

Choose **with-engine** for the simplest first installation. It works without a
second engine download, and **OpenQP (bundled)** records the engine version used
for a calculation.

Choose the **standard Studio package** when a local development build of OpenQP
must be tested, disk space is limited, or engine updates should be managed
separately. In **Execution**, select **OpenQP (local)** and verify the displayed
path and version before starting a job.

Choose the **standalone engine** when no graphical interface is required. It can
also supply the local runner for the standard Studio package after its
`openqp` executable is placed on the application's search path.

## Release asset names

Replace `0.2.3` and `1.3.1` below with the versions shown on the
[OpenQP releases page](https://github.com/Open-Quantum-Platform/openqp/releases).

### macOS

| Machine | Standard Studio | Integrated Studio | Standalone engine |
| --- | --- | --- | --- |
| Apple Silicon (`arm64`) | `OQP-Studio-0.2.3-macos-apple-silicon.dmg` | `OQP-Studio-0.2.3-macos-apple-silicon-with-engine.dmg` | `openqp-1.3.1-macos-arm64.zip` |
| Intel (`x86_64`) | `OQP-Studio-0.2.3-macos-intel.dmg` | `OQP-Studio-0.2.3-macos-intel-with-engine.dmg` | `openqp-1.3.1-macos-x86_64.zip` |

Each macOS Studio build also has an `.app.tar.gz` alternative for terminal
installation.

### Windows and Linux

| Platform | Standard Studio | Integrated Studio | Standalone engine |
| --- | --- | --- | --- |
| Windows x64 | `OQP-Studio-0.2.3-windows-x64-setup.exe` | `OQP-Studio-0.2.3-windows-x64-with-engine-setup.exe` | `openqp-1.3.1-windows-x86_64.zip` |
| Debian/Ubuntu x86_64 | `OQP-Studio-0.2.3-linux-x86_64.deb` | `OQP-Studio-0.2.3-linux-x86_64-with-engine.deb` | `openqp-1.3.1-linux-x86_64.tar.gz` |
| Other x86_64 Linux | `OQP-Studio-0.2.3-linux-x86_64.AppImage` | Use the standard AppImage with a standalone engine | `openqp-1.3.1-linux-x86_64.tar.gz` |
| Linux AArch64 | No Studio installer | No Studio installer | `openqp-1.3.1-linux-aarch64.tar.gz` |

The integrated Linux package is currently distributed as a `.deb`; there is no
`with-engine` AppImage. The integrated Windows package uses the interactive
`.exe` installer; the standard package additionally provides an `.msi` for
managed deployment.

## Verify a download

Download `SHA256SUMS` from the same release. On macOS, verify one selected file
with:

```bash
FILE=OQP-Studio-0.2.3-macos-apple-silicon-with-engine.dmg
grep "  $FILE$" SHA256SUMS | shasum -a 256 -c -
```

On Linux:

```bash
FILE=OQP-Studio-0.2.3-linux-x86_64-with-engine.deb
grep "  $FILE$" SHA256SUMS | sha256sum -c -
```

On Windows PowerShell:

```powershell
Get-FileHash .\OQP-Studio-0.2.3-windows-x64-with-engine-setup.exe -Algorithm SHA256
```

Compare the displayed Windows hash with the corresponding line in
`SHA256SUMS`.

## Use the standalone engine

The standalone archive is self-contained; a separate Python or OpenQP
installation is not required.

### macOS

```bash
mkdir -p "$HOME/Applications/openqp-1.3.1"
cd "$HOME/Applications/openqp-1.3.1"
unzip "$HOME/Downloads/openqp-1.3.1-macos-arm64.zip"
OMP_NUM_THREADS=4 ./openqp calculation.oqp
```

Use the `macos-x86_64` archive on an Intel Mac. If macOS quarantines the
extracted files after a browser download, verify the archive first and then run:

```bash
xattr -dr com.apple.quarantine "$HOME/Applications/openqp-1.3.1"
```

### Linux

```bash
mkdir -p "$HOME/opt"
tar xzf openqp-1.3.1-linux-x86_64.tar.gz -C "$HOME/opt"
OMP_NUM_THREADS=4 "$HOME/opt/openqp/openqp" calculation.oqp
```

### Windows PowerShell

```powershell
New-Item -ItemType Directory -Force "$HOME\OpenQP\1.3.1" | Out-Null
Expand-Archive .\openqp-1.3.1-windows-x86_64.zip "$HOME\OpenQP\1.3.1"
$env:OMP_NUM_THREADS = "4"
& "$HOME\OpenQP\1.3.1\openqp.exe" .\calculation.oqp
```

The calculation writes its `.log`, `.json`, `.molden`, trajectory, and other
requested output files beside the input unless the input or caller specifies a
different working directory.

## Use the integrated package

Install the `with-engine` asset using the platform instructions in
[Installing OQP Studio](installation.md). Then:

1. Open **Execution**.
2. Select **OpenQP (bundled)**.
3. Confirm that the bundled runner has a check mark and a version.
4. Select the number of threads and results directory.
5. Run the generated `.oqp` input.

The bundled engine is private to that application installation. Updating a
system Python package or another `openqp` executable does not replace it. To
update both Studio and the bundled engine, install a newer `with-engine`
package. Existing projects remain in the selected results directory.

## Use the standard Studio package with a local engine

Install OpenQP separately with the package instructions in the main
[Installation](../installation.md) chapter, or extract a standalone engine and
make its executable available to Studio. Restart Studio after changing the
search path. In **Execution**, select **OpenQP (local)** only after checking the
displayed executable path and version.

The standard package can also obtain a managed engine through
**Execution > Install compute engine** when that command is available. This
requires network access. It is not needed by a `with-engine` installation.

## Updating and removing packages

- Replacing Studio does not delete projects stored outside the application.
- Replacing a standalone engine directory does not update an engine embedded in
  a `with-engine` application.
- Removing Studio does not remove a separately installed local or standalone
  OpenQP engine.
- Keep the old engine when a calculation must be reproduced with its original
  version; select the corresponding local executable in Execution.

Record the Studio version, runner name, OpenQP version, input, and output files
with any calculation that must be reproduced.
