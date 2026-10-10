# Installing OQP Studio

Installers are published on the
[OpenQP releases page](https://github.com/Open-Quantum-Platform/openqp/releases).
Each compatible OpenQP release carries the standalone engine, standard Studio
application, and integrated Studio application together.
For an application that can calculate without downloading OpenQP later, choose
an asset whose name contains `with-engine`.

Before downloading, see [Standalone and integrated packages](packages.md) for
the distinction between the standard Studio application, the standalone
command-line OpenQP engine, and the `with-engine` integrated application.

## macOS

Choose the architecture that matches `uname -m`:

| `uname -m` | Asset label |
| --- | --- |
| `arm64` | `macos-apple-silicon` |
| `x86_64` | `macos-intel` |

The `.dmg` is the normal graphical installer. Open it and drag **OQP Studio** to
**Applications**. Studio 0.2.4 requires macOS 15 or later. Its ad-hoc signature
is verified during packaging, but the app is not Apple-notarized. If the first
launch is blocked, verify the download checksum, then use **System Settings →
Privacy & Security → Open Anyway** for this app.

The `.app.tar.gz` asset can instead be installed entirely from the terminal.
Files downloaded by `curl` do not receive the browser quarantine attribute:

```bash
cd ~/Downloads
curl -L -o oqp-studio.tar.gz "<asset URL>"
tar xzf oqp-studio.tar.gz -C /Applications
```

On the first calculation, macOS may ask whether OQP Studio may access the
Documents folder. Grant access if results should use the default
`~/Documents/OQP Studio/jobs` location, or choose another folder from
**File > Results folder**.

## Windows

Choose the `windows-x64-with-engine-setup.exe` installer for the complete
package. Windows may display a SmartScreen warning while builds remain
unsigned. Verify that the installer came from the official release before
choosing **More info > Run anyway**.

Studio can use the bundled native Windows engine, a local OpenQP on `PATH`, or
OpenQP in WSL. Native bundled execution is the simplest default.

Native Windows OpenQP builds currently do not support ddX. Use an OpenQP build
inside WSL for PCM/ddX calculations; other workflows can use the native bundled
engine.

## Linux

The 0.2.4 desktop packages target Ubuntu 24.04 or a compatible newer Linux
distribution (glibc 2.39+, GTK 3 and WebKitGTK 4.1). Debian 13 is a compatible
base; older RHEL/Rocky systems can use the standalone engine or a remote engine
without running the Studio desktop locally. A `.deb` can be
installed with:

```bash
VERSION=0.2.4  # replace when installing a later release
sudo apt install "./OQP-Studio-${VERSION}-linux-x86_64-with-engine.deb"
```

For an AppImage, make the downloaded file executable before opening it:

```bash
VERSION=0.2.4  # replace when installing a later release
chmod +x "OQP-Studio-${VERSION}-linux-x86_64.AppImage"
"./OQP-Studio-${VERSION}-linux-x86_64.AppImage"
```

## Engine choices

**OpenQP (bundled)** uses the engine shipped with a `with-engine` installer.
It is isolated from changes to the shell environment and is the reproducible
choice for a new installation.

**OpenQP (local)** uses the `openqp` command found on the application's search
path. Select it when testing a separately installed or development OpenQP
build. Studio displays the detected path and version in Execution.

**WSL** is available on Windows when a usable OpenQP installation is detected
inside WSL. An unavailable runner remains visible as unavailable and cannot be
selected.

The plain installer obtains an engine on demand through **Execution > Install
compute engine**. A `with-engine` installer avoids that network step.

For complete standalone extraction commands, platform-specific asset names,
checksum verification, and update behavior, see
[Standalone and integrated packages](packages.md).

## Python package

The server and interface can also be installed from PyPI:

```bash
python -m pip install oqp-studio
oqp-studio
```

Use `python -m pip install "oqp-studio[desktop,chem]"` to add the native
pywebview window and RDKit-based 2D-sketch-to-3D conversion. The Tauri installer
is recommended for ordinary desktop use.
