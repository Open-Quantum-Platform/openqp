# Studio responsiveness audit — 2026-09-22

Scope: desktop request transport, backend startup/discovery, job submission,
Art scene preparation/rendering, and existing backend/frontend checks. This
is not a claim that every scientific workflow or every GPU has been tested.

## Findings and changes

| Finding | User-visible effect | Change |
| --- | --- | --- |
| Synchronous Tauri command waits up to 180 seconds for the sidecar | Desktop main thread cannot paint or respond during backend startup/API calls | Async command with blocking wait on a worker; sidecar termination rejects pending requests |
| Stdio backend creates a new event loop and handles one request at a time | A slow runner probe delays health and cancel; FastAPI lifespan never runs | One lifespan/event loop, up to eight concurrent requests, replies matched by ID |
| No startup feedback or engine discovery error display | Slow launch looks hung; backend version mislabeled as engine version | Rotating clock and elapsed seconds for startup/discovery/submission; explicit failure state; separate Studio and runner version labels |
| Frontend probes runners twice; failed WSL probes raise | Duplicate cold WSL startup; one timeout can hide all runners | One detail request; WSL timeout means unavailable, preserving native runners |
| Engine version detection reads an entire executable if README is absent | Large binaries decoded/allocated just to check a shebang | Read at most 4096 bytes |
| Art starts ray tracing automatically at device pixel ratio up to 2, without limits | High GPU load persists and can make slower computers unresponsive | Raster preview by default; explicit ray opt-in; capped pixel count, bounces, samples and duration |
| Workaround forces synchronous WebKit shader compilation | A driver compile can block without letting a timer/cancel button run | Remove workaround; use raster-only safe mode on WebKit, software/unknown GPUs, <=2 logical CPUs, or missing parallel shader compilation |
| Synchronous BVH builds, repeated slider rebuilds, rendering in CSS-hidden iframe | Input changes and navigation keep doing expensive work | Worker BVH; cancellable preparation; return to preview on scene/material/quality changes; parent sends active state and defers hidden scenes |
| Submit button remains usable during resource preflight | Double clicks can start duplicate jobs and overload the computer | Disable submission for the full preflight/submit operation, recover on failure |

Loading a new cube disables ray tracing and export until that scene is ready.
Rejected oversized scenes clear retained geometry and cannot reappear after
appearance changes.

## Rendering bounds

- Preview: at most approximately 1 million physical pixels; pixel ratio <=1.25;
  redraw only on change, at most 20 frames/s.
- Ray tracing: at most approximately 260,000 pixels, <=4 bounces, 8x8 tiles,
  explicit start only. Stop at 64 samples or 20 seconds; Resume starts another
  bounded interval. Preparation timeout 15 seconds; Cancel restores preview.
- A CPU-side render call taking >180 ms restores preview. This is an additional
  guard, not a GPU-driver watchdog. JavaScript cannot interrupt a stuck driver.
- Art rejects >300 atoms and caps surface sampling at 60,000 cells. Large
  structures remain usable in Analysis. Cube grid validation remains bounded;
  oversized cube text (>24 MB) is rejected for Art.
- Hidden/inactive Art stops tracing. Inactive scene updates are held until the
  user returns. Software GPUs and macOS Tauri/WebKit remain in interactive
  rendering; ray-traced reflections/depth of field require a supported device.

## Verification and limits

Local checks include the complete backend suite, frontend TypeScript/Vite
build, render-budget/cancellation tests, live stdio health vs runner discovery,
macOS `cargo check --locked`, and browser interaction under simulated slow CPU
with a software GPU. Native checking uses an inert sidecar placeholder only
for compile-time resource validation; it does not build or execute an installer.

The existing macOS package uses an unpacked sidecar. Windows/Linux still use
PyInstaller onefile extraction, which can add cold-start delay. This audit
makes waiting visible and nonblocking; actual installed cold-start timing on
those platforms and hardware ray-tracing qualification remain release checks.
The large Ketcher bundle is a separate lazy-loaded sketcher entry; its first
open should be profiled separately before changing chemistry functionality.

Source checks do not claim a native installer startup speedup in seconds.
No scientific engine build, calculation, or public installer release is part
of this change. Public release remains subject to the previously agreed NMR
PR/MR prerequisite.

References: [Tauri command threading](https://v2.tauri.app/develop/calling-rust/)
and [path-tracer async BVH and rendering controls](https://github.com/gkjohnson/three-gpu-pathtracer).

## Optional browser regression harness

With the isolated backend on 127.0.0.1:18814 and Vite on 127.0.0.1:15173
(proxying that backend), run `frontend/tests/browser-responsiveness.cjs`
using Node. Install Playwright in a test environment or point
`PLAYWRIGHT_MODULE` to an existing Playwright package; `CHROME_EXECUTABLE`
selects an existing Chrome executable, and `STUDIO_TEST_URL` overrides Vite's
URL. Run from a disposable output directory: the harness writes two screenshots.
It forces a software GPU and 6x CPU throttling, delays the health response,
and verifies startup feedback, navigation, safe rendering, and rejection of
oversized structures without redisplaying a previously loaded orbital surface.
