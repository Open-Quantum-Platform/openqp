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

The current historical Studio 0.2.3 downloads are listed on the
[download page](download.md). Their location under the OpenQP 1.3.1 release
does not establish the exact engine commit inside an installer; consult the
artifact's own metadata. New releases must supply that provenance explicitly.

## Selecting another engine

When using **OpenQP (local)** or a remote engine, the engine actually selected
for a calculation can differ from the one bundled with Studio. Check its
version and required capabilities before running a new workflow. Updating the
interface alone does not add capabilities to an older engine.

A refresh of the bundled engine will be distributed under a new Studio release
identifier. Previously published installers will not be replaced under the
same version.
