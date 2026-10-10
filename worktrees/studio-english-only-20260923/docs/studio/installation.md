# Installing OQP Studio

Get Studio installers from the [Studio download page](download.md) or the
[independent Studio releases](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases).
Studio 0.2.4 and later are distributed separately from OpenQP engine releases.
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
**Applications**, then eject the disk image. Launch the installed copy from
Applications, rather than the copy inside the disk image. Studio 0.2.4 requires
macOS 15 or later; both standard and `with-engine` installers use the same
first-launch procedure below.

### macOS security: first launch

Studio 0.2.4 has an ad-hoc signature checked during packaging, but it does not
have an Apple Developer ID signature or Apple notarization. macOS may therefore
block the first launch with an **unidentified developer** or **cannot check for
malicious software** warning. An ad-hoc signature is not an Apple approval.

1. Download from the official Studio release linked above and
   [verify the download against SHA256SUMS](packages.md#verify-a-download).
2. Double-click **Applications → OQP Studio** once. If macOS blocks it, dismiss
   the warning with **Done** or **Cancel**; keep the app installed.
3. Open **Apple menu → System Settings → Privacy & Security**,
   then scroll to **Security**.
4. Find the message about **OQP Studio** and click **Open Anyway**.
   Confirm **Open** in the next dialog and authenticate
   with your Mac login password or Touch ID if asked.
5. Once approved, open Studio normally from Applications on subsequent launches.

If **Open Anyway** is missing, try opening the installed app again and return to
Settings immediately; Apple documents a roughly one-hour availability window.
On a managed Mac, ask the administrator if organization policy prevents an
exception. See Apple's [first-launch guidance](https://support.apple.com/en-us/102445)
and [unknown-developer instructions](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unknown-developer-mh40616/mac).

### If the app is still blocked or reported as damaged

A **damaged** message can also mean an incomplete download or altered app.
Download a fresh official copy, check its SHA256SUMS entry, and copy that app to
Applications again. Do not treat a failed checksum or signature check as a
quarantine problem.

For a verified official copy that remains blocked on your own Mac, the following
Terminal fallback removes the download quarantine attribute **only from OQP
Studio**, after checking the installed bundle's signature. Quit Studio first.
This is a user-selected exception for this app; it does not notarize the app.

```bash
APP="/Applications/OQP Studio.app"
codesign --verify --deep --strict --verbose=2 "$APP" &&
  xattr -dr com.apple.quarantine "$APP" &&
  open "$APP"
```

The `&&` operators stop the sequence if either check or attribute removal fails.
If you installed in your personal Applications folder, change the first line to
`APP="$HOME/Applications/OQP Studio.app"`. If permission is denied, install the
verified app there using Finder and retry with that path. Keep the app path
quoted because its name contains a space.

If `codesign` reports an invalid signature, use a fresh verified download or
report the error; do not re-sign the installed app to hide the failure. For an
explicit malware or revoked-authorization warning, stop and contact the
maintainer. Keep Gatekeeper and System Integrity Protection enabled; no
system-wide security change is required by these instructions.

The `.app.tar.gz` release asset is an alternative to the DMG. Verify its checksum
before extracting it with Archive Utility and moving **OQP Studio.app** to
Applications, then follow the same first-launch instructions.

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
