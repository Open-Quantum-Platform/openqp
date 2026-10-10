# Release notes and version comparison

[Download Studio](download.md){ .md-button .md-button--primary }

Studio releases are numbered independently of OpenQP engine releases. Each
entry below states whether it is available to download or still in development.
The current public desktop release is **Studio 0.2.3**, released on
**23 August 2026**.

## Next Studio release — unreleased

The following changes are implemented in the development source reviewed for
this page on **23 September 2026**. They are **not included in the current
0.2.3 downloads**. The next version number, release date and bundled engine
commit will be assigned after installer qualification.

### Compared with Studio 0.2.3

| Area | Studio 0.2.3 | Next release, in development |
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
| Engine identification | An engine package version does not uniquely identify a later development snapshot | Release information will distinguish Studio version, engine version and exact engine commit |

### Qualification still required

Windows, macOS and Linux installer startup, upgrade and bundled-engine checks
must pass before these changes are announced as a release. Source-level tests
and a frontend build alone do not qualify an installer. No measured installer
startup speedup is claimed yet.

Ray-tracing limits reduce sustained load, but they cannot interrupt a GPU driver
that has stopped responding. Preview remains the supported fallback. WebKit
uses preview-only rendering under the new policy.

### ACID current-density maps — integrated, unreleased

The results viewer now includes ACID current-density vectors with arrow-density
and arrow-length controls. Direct cube opening and the surface panel both load
the associated current components. These source changes have been reviewed and
integrated after backend tests and frontend checks on Windows, macOS and Linux.

Installer startup and bundled-engine compatibility still need qualification
before release. This is not a capability promised for the 0.2.3 installer.

### Upgrade notes for the next release

- Existing 0.2.3 installers remain downloadable during the transition.
- The first independent release may require a manual download: an already
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
