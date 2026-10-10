# Release notes and version comparison

[Download Studio](download.md){ .md-button .md-button--primary }

Studio releases are numbered independently of OpenQP engine releases. Each
entry below states whether it is available to download or still in development.
The current public desktop release is **Studio 0.2.5**.

## Studio 0.2.5 — Stable

**Released:** 23 September 2026. **Previous version:** 0.2.4.
**Engine:** the same qualified OpenQP 1.3.1 snapshot used in Studio 0.2.4,
commit [`1db79089b285`](https://github.com/Open-Quantum-Platform/openqp/commit/1db79089b285406338bd28a255b967e4011f8df7).
Exact application and packaging commits and installer checksums are recorded in
the [release manifest](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.5/manifest.json).

### Compared with Studio 0.2.4

| Area | Studio 0.2.4 | Studio 0.2.5 |
| --- | --- | --- |
| NMR setup | Gauge selection can appear unset | GIAO is explicitly selected by default |
| ACID discovery | Current-density tools are harder to find | Dedicated ACID workflow, searchable as ACID or AICD |
| Atomic-property maps | Overlapping molecular representations and strong colors | One molecular representation, thin neutral bonds and a subdued blue–ivory–terracotta scale with a matching legend |
| Labels and camera | Updating a property can reset the view | Optional three-decimal labels and preserved camera position |
| Viewer controls | An initialization probe control is visible | Unnecessary probe control removed |

The NMR atomic-property view displays **atomic shielding values**, not a
continuous shielding field. The bundled engine has not changed in this release;
upgrading Studio does not add new capabilities to another local or remote engine.

Native standard and with-engine installers are available for macOS Apple Silicon,
macOS Intel, Windows x64 and Linux x86_64. All four targets passed installer and
backend checks and NMR/ACID engine qualification. Mac packages require macOS 15
or later and are ad-hoc signed, not Apple-notarized; see the
[first-launch security instructions](installation.md#macos-security-first-launch).
Windows installers are unsigned.

Upgrade using the same operating system, architecture and installer variant.
Studio 0.2.3 users must install the latest package manually once because that
version points to the former update repository. Published 0.2.4 installers
remain available unchanged.

## Studio 0.2.4 — previous stable release

**Previous version:** 0.2.3. **Engine:** OpenQP 1.3.1 snapshot,
commit [`1db79089b285`](https://github.com/Open-Quantum-Platform/openqp/commit/1db79089b285406338bd28a255b967e4011f8df7).
Exact source commits and installer checksums are recorded in the
[release manifest](https://github.com/Open-Quantum-Platform/oqp-studio-releases/releases/download/v0.2.4/manifest.json).

### Compared with Studio 0.2.3

| Area | Studio 0.2.3 | Studio 0.2.4 |
| --- | --- | --- |
| Startup and engine discovery | A slow backend or engine probe can appear unresponsive | Rotating clock and elapsed time during startup, discovery and submission; explicit failure messages |
| Interface responsiveness | Waiting for a backend response can block desktop interaction; a slow probe can delay other requests | Backend waits run off the desktop main thread; requests can progress independently |
| Repeated job submission | Submission can be clicked again while resource checks are pending | Submission remains disabled through the complete check-and-submit operation |
| Art rendering | Ray tracing starts automatically and can keep a slow GPU busy | Interactive preview by default; ray tracing starts explicitly and uses bounded resolution, samples and duration |
| Slow or unsupported graphics hardware | Ray tracing can cause long stalls | Preview-only mode on WebKit and unsupported configurations; slow rendering returns to preview where possible |
| Scene replacement | Expensive preparation or rejected input can leave stale artwork | Cancellable preparation; tracing and export wait for a new cube; rejected oversized scenes clear retained artwork |
| Navigation | Numbered steps such as Builder 1–5 | Feature names without step numbers |
| Application identity | Previous desktop icon | Pink-and-blue orbital Q icon in the desktop packages and header |
| Downloads | Studio installers are attached to an OpenQP engine release | Dedicated public Studio download page and independent Studio release numbering |
| Engine identification | An engine package version does not uniquely identify a later development snapshot | Release information distinguishes Studio version, engine version and exact engine commit |

### Supported distributions

Native installers cover Windows x64, Linux x86_64, and macOS on Apple Silicon
and Intel, with standard and with-engine variants. No measured installer
startup speedup is claimed.

Ray-tracing limits reduce sustained load, but they cannot interrupt a GPU driver
that has stopped responding. Preview remains the supported fallback. WebKit
uses preview-only rendering under the new policy.

### ACID/AICD current-density maps

The results viewer now includes ACID current-density vectors with arrow-density
and arrow-length controls. Direct cube opening and the surface panel both load
the associated current components. These source changes have been reviewed and
integrated after backend tests and frontend checks on Windows, macOS and Linux.

The bundled engine is qualified with NMR shielding and ACID calculations,
including the ACID map and all three induced-current component cube grids.
Older or unverified local/remote engines cannot run ACID through Studio until
their capabilities are identified. This is not included in the 0.2.3 installer.

### Upgrade notes

- Existing 0.2.3 installers remain downloadable during the transition.
- Upgrading from 0.2.3 requires a manual download: an already
  installed updater that points to the previous repository cannot be repaired
  by changing this website alone.
- Select the same installer variant when upgrading: **standard** or
  **with-engine**, with the matching operating system and architecture.
- A different local or remote OpenQP installation does not gain new engine
  capabilities when Studio is upgraded. Check
  [engine compatibility](engine-compatibility.md) before using a new workflow.

## Studio 0.2.3 — available

**Released:** 23 August 2026. **Previous version:** 0.2.2.

### Changes from 0.2.2

| Area | Added or corrected in 0.2.3 |
| --- | --- |
| Calculation setup | One canonical calculation driver and portable end-to-end workflow recipes |
| Molecular exploration | Bond scans, interactive reaction paths and MRSF excited-state maps |
| Comparing results | Comparison of calculations and display of atomic properties |
| Surface analysis | Volumetric surface tools and molecular symmetry analysis |
| Art | Configurable rendering, molecular and volumetric ray tracing |
| Rendering reliability | Hide coarse intermediate cube surfaces; fix WebKit tracing that stays at zero samples; keep appearance changes from repeatedly resetting tracing |
| Distribution | Harden packaged distribution and align licensing with OpenQP |
| Bundled engine limitation | ddX is disabled in the 0.2.3 engine bundles |

The [historical public download release](https://github.com/Open-Quantum-Platform/openqp/releases/tag/v1.3.1)
hosts Studio 0.2.3 installers alongside OpenQP 1.3.1 assets. These are two
different product version numbers. See the [download page](download.md) for
platform-specific installer links and checksums.

## Earlier versions

Studio **0.2.2** was released on **22 August 2026**, and **0.2.1** on
**21 August 2026**. Their release records do not contain detailed release
notes; no additional feature differences are inferred here. The 0.2.3
comparison above is based on the source changes between its two release tags.

## What each future release note will contain

Every independent release will state its Studio version and date, previous
Studio version, Stable or Preview channel, and bundled OpenQP version and exact
commit. Its notes will include a **previous version versus new version** table,
new features, fixes, known limitations, supported installer variants and upgrade
instructions. Only changes included in the published artifacts will move from
Unreleased into that version's entry.
