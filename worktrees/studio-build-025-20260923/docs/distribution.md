# Distribution channels

See [Studio releases and engine compatibility](studio-releases.md) for the
current repository policy and independent distribution proposal.

Studio source is private. The existing public desktop installers are attached
to [OpenQP v1.3.1](https://github.com/Open-Quantum-Platform/openqp/releases/tag/v1.3.1).
GitHub Actions are disabled. The old tag-triggered wheel/sdist/PyPI publisher
predates the source-privacy decision and must not be used for public promotion.

## Historical platform constraints


**Google Play** distributes Android apps only. OQP Studio is a desktop
application built on Tauri (a native webview plus a Python backend that spawns
OpenQP), so there is nothing to submit. An Android build would mean a
different product: a thin client talking to a remote OpenQP server, since
phones cannot run the Fortran core.

**Mac App Store** requires the App Sandbox. The Studio's whole purpose is to
execute a separately installed `openqp` binary — from a conda environment,
`/usr/local/bin`, or WSL — and to read and write job directories the user
chooses. A sandboxed app may not launch executables outside its own bundle,
and Apple rejects apps that depend on separately installed command-line
tools. Passing review would mean removing the ability to run calculations.

**Microsoft Store** is possible: registration is now free, and unpackaged
`.exe`/`.msi` installers have been accepted since 2021. But submissions must
be Authenticode-signed by a certificate chaining to a Microsoft-trusted CA,
so the Store does not avoid the signing cost (see
[code-signing.md](code-signing.md)). It is worth doing only alongside a
Windows certificate.
