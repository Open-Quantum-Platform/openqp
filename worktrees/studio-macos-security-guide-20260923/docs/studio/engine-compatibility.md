# Studio and engine compatibility

Studio and OpenQP have independent version numbers. A newer Studio interface
may require engine functionality that is newer than the latest formal OpenQP
release. A tested engine snapshot can therefore be bundled with Studio without
waiting for matching OpenQP and Studio version numbers.

For each new Studio release, the download information will distinguish:

| Identifier | Meaning |
| --- | --- |
| Studio version | The desktop application release |
| Bundled OpenQP version | The engine's package version |
| OpenQP commit | The exact engine source used to build every platform |
| Stable or Preview | Whether the pair is the recommended release or an opt-in candidate |
| Platform and architecture | The operating system and CPU supported by that installer |

An engine built after its formal release must be identified as a snapshot,
with its exact commit, rather than described as the unchanged formal release.
The source identity must be the same across all platform builds of a candidate.

## Studio 0.2.4 Stable

| Component | Qualified identity |
| --- | --- |
| Studio | 0.2.4; exact source commit in the manifest and About |
| Bundled OpenQP | 1.3.1 snapshot |
| Engine commit | `1db79089b285406338bd28a255b967e4011f8df7` |
| Magnetic properties | NMR shielding, ACID/AICD maps and induced-current vectors |
| macOS | Apple Silicon or Intel; macOS 15 or later |
| Windows | x64 |
| Linux desktop | x86_64; Ubuntu 24.04+/compatible glibc 2.39+, GTK 3, WebKitGTK 4.1 |
| Standalone Linux engine | x86_64; built on glibc 2.28 |

[Download Studio](download.md). The package version 1.3.1 does not imply that
this later snapshot equals the formal OpenQP v1.3.1 release. Every platform uses
the same engine commit above.

Historical Studio 0.2.3 remains available on the OpenQP v1.3.1 release page. Its
location alone does not establish the exact engine commit inside an installer.

## Selecting another engine

When using **OpenQP (local)** or a remote engine, the engine actually selected
for a calculation can differ from the one bundled with Studio. Check its
version and required capabilities before running a new workflow. Updating the
interface alone does not add capabilities to an older engine.

A refresh of the bundled engine will be distributed under a new Studio release
identifier. Previously published installers will not be replaced under the
same version.
