# Download OQP Studio

Choose your operating system below. The **with-engine** installer includes
OpenQP and works without a separate engine download. A **standard** installer
uses a separately installed local or remote engine.

## Stable release: Studio 0.2.4

[Read release notes and compare versions](whats-new.md).

Studio **0.2.4** is released independently of OpenQP. The bundled engine is an
OpenQP **1.3.1 snapshot**, commit
[`1db79089b285`](https://github.com/Open-Quantum-Platform/openqp/commit/1db79089b285406338bd28a255b967e4011f8df7),
including NMR shielding and ACID/AICD current-density maps and vectors.
Both version numbers and the exact engine commit appear in **About**.

| Platform | Complete installation | Standard Studio |
| --- | --- | --- |
| macOS, Apple Silicon | [With engine (.dmg)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-macos-apple-silicon-with-engine.dmg) | [Standard (.dmg)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-macos-apple-silicon.dmg) |
| macOS, Intel | [With engine (.dmg)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-macos-intel-with-engine.dmg) | [Standard (.dmg)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-macos-intel.dmg) |
| Windows, x64 | [With engine (.exe)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-windows-x64-with-engine-setup.exe) | [Standard (.exe)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-windows-x64-setup.exe) |
| Linux, x86_64 (Ubuntu 24.04+/Debian 13+) | [With engine (.deb)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-linux-x86_64-with-engine.deb) | [Standard (.deb)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-linux-x86_64.deb) |
| Linux, x86_64 (AppImage) | Use the standard package with a local engine | [Standard (.AppImage)](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/OQP-Studio-0.2.4-linux-x86_64.AppImage) |

[Download SHA256SUMS](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/SHA256SUMS)
and follow the [verification instructions](packages.md#verify-a-download).
See [installation instructions](installation.md) and
[all release assets](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/tag/v0.2.4)
for alternative formats.

## Updates and historical releases

Stable is the default update channel. Preview is an explicit opt-in in **About**;
pre-release builds do not silently replace Stable. Every engine update receives
a new Studio release identifier and exact engine commit.

**Upgrading from 0.2.3:** manually download and install 0.2.4 once. The old updater
points to the former private repository. Keep the same operating system,
architecture and installer variant.

Historical [Studio 0.2.3 downloads](https://github.com/Open-Quantum-Platform/openqp/releases/tag/v1.3.1)
remain available with unchanged files and URLs.
See [engine compatibility](engine-compatibility.md) before selecting a different
local or remote engine.
