# Handoff

> **Status note (2026-08-22):** This is a historical v0.2.1 handoff, not a
> description of current `main`. Release v0.2.2 is now published. Open problem
> 1 was addressed by `d9d51e5` (an unpacked macOS sidecar), and open problem 2
> by `6092931` (a pinned OpenQP commit in every engine build). Open problem 3
> remains an upstream follow-up.

Written 2026-08-21 (KST), handing this repository from a remote agent to one
working on a local Mac. It says where the project stands, what is unfinished,
and the few conventions that were learned the hard way.

The reason for the handoff is worth stating, because it should shape how the
next stretch of work is done: verifying a change here meant dispatching CI and
waiting the better part of an hour, so every wrong guess cost an hour. A local
Mac can measure in a minute what took a day to argue about. Measure first.

## Where things stand

- `main` at `2f07579`.
- One published release: **v0.2.1**. The fourteen before it were deleted (tags
  kept) once v0.2.1 was verified.
- v0.2.1 carries 21 assets: 8 slim installers, 6 all-in-one installers, 5
  engine bundles, wheel + sdist.
- The engine in it is OpenQP `v1.3.1` = commit `b039efc`.
- `2f07579` (the basis-set name fix, below) is **not in any release yet**.

The application: a Tauri 2 shell around a PyInstaller-frozen FastAPI sidecar.
One origin serves the UI, the viewer and the API. The compute engine is
OpenQP, either bundled inside the installer (`-with-engine`) or fetched on
first use.

## Conventions

**Never use a release to test a change.** Dispatching `release.yml` with an
empty `release_tag` runs all fifteen jobs and publishes nothing. Iterate on
that until it is clean; only then publish. This is not a style preference — a
no-tag run has already caught two defects that would otherwise have shipped.

**Results belong in the user's Documents folder**, at `~/Documents/OQP
Studio/jobs`, never inside the app bundle. `File > Results folder…` moves it.
The path is chosen lazily (`workspace.preferred()` touches no filesystem;
`workspace.ensure()` creates and write-probes, falling back down a chain).
Nothing may create that directory at import time.

**A version bump touches five files**, and missing one produces an installer
that disagrees with itself:

    backend/pyproject.toml
    backend/oqp_studio/__init__.py
    shell/src-tauri/Cargo.toml
    shell/src-tauri/Cargo.lock        # the `name = "oqp-studio"` entry
    shell/src-tauri/tauri.conf.json

**Pruning releases** is `prune-releases.yml`, which requires `keep=<tag>`,
`confirm=delete`, and optionally `delete_tags`. It refuses to run unless the
tag to keep exists.

## Open problem 1 — macOS startup takes twenty seconds

This is the first thing to work on, and it is unfinished because it was never
measured properly.

In v0.2.0 the backend took **14 s and 20.181 s** on real hardware to answer
`/api/health`, against a 30 s limit — close enough to the limit that a slow
disk read as "OQP Studio backend failed to start". v0.2.1 responded to that:
the login-shell PATH probe and the RDKit warm-up came out of the port-opening
path, a splash with a running counter went in, the limit went to 180 s, and
`main.py` now prints a timing trace at each startup step.

What is still unknown is where the twenty seconds go. Two candidates:

1. **PyInstaller `--onefile` extraction.** `backend/build_binary.py` passes
   `--onefile`, so every launch unpacks an interpreter plus NumPy, SciPy and
   the whole of RDKit into a temporary directory before Python runs at all.
   If this is it, the fix is `--onedir`, shipped as a Tauri resource:
   `build_binary.py`, the "Build backend sidecar binary" and "Stage sidecar
   with target triple" steps in `.github/workflows/installers.yml`, and
   `externalBin`/`resources` in `shell/src-tauri/tauri.conf.json`.
2. **Import cost** inside the frozen app, in which case the expensive imports
   move behind lazy loading.

### How to tell which

The decisive number is **how long before the first trace line appears**. That
interval is time spent before Python executes a single line, so it is
extraction, not imports.

Run the frozen sidecar straight from a terminal — this isolates it from Tauri
and is the cleanest measurement available:

    time "/Applications/OQP Studio.app/Contents/MacOS/oqp-studio-backend" --port 8899

It prints, to stderr:

    startup   1.23s  imports
    startup   1.24s  PATH
    startup   1.31s  TLS
    startup   1.35s  app built

Read it as: the delay *before* `startup ... imports` appears is extraction plus
interpreter start; the gaps *between* lines are import and setup cost. A large
first number and small gaps means `--onedir`. Small first number and a large
gap means the step named by the following line is the culprit.

**Do not look for these lines in `backend.log`.** That file is only written
when startup raises — `server_main.py` catches `BaseException` and dumps the
traceback there. The trace lines go to stderr, and the shell sidecar's output
receiver is dropped in `main.rs` (`let (_rx, child) = …`), so when the app is
launched from Finder nobody is listening. Worth fixing on its own: forwarding
that receiver into `backend.log` would make startup observable in the field
instead of only on a developer's machine.

On macOS `backend.log` lives at
`~/Library/Application Support/OQP Studio/backend.log`.

### A warning from the last attempt

The first diagnosis of this symptom was wrong, and confidently so: the blame
went to an import-time `mkdir` in `~/Documents` blocking on a macOS TCC
permission prompt. A manual test on real hardware showed no dialog, the
directory created normally, and the backend simply taking twenty seconds. The
lazy-Documents change that came out of it is still correct, but it did not fix
the symptom, and an hour of CI went into shipping it. Measure before changing.

## Open problem 2 — the engine reference moves

The three engine workflows resolve OpenQP's version at build time:

    tag=$(gh release view --repo Open-Quantum-Platform/openqp --json tagName -q .tagName)

OpenQP's `v1.3.1` tag moved three times in one day (`9d9cc53` → `35ee885` →
`b039efc`), so which source a bundle was built from depended on what hour the
build ran. v0.2.1 happened to catch the final commit; that was luck, and it
cost a round of forensics to establish.

`linux-bundle.yml`, `macos-bundle.yml` and `windows-bundle.yml` already accept
an `openqp_ref` input; `release.yml` never passes one. Pin it to a commit SHA
there before the next release.

To check afterwards what a bundle actually contains: each engine archive has a
`README.txt` recording `OpenQP commit : <sha>`, and the "Add the data tree and
a README" build step echoes its last five lines.

## Engine requirement — external ddX until native PCM lands

The current OpenQP source supports PCM through external ddX.  The native ddX
implementation exists but is not merged, so all three standalone-engine
workflows must keep `ENABLE_DDX=ON`, add the resulting ddX shared library to
the PyInstaller payload, and run the DDPCM water example in the clean-room
smoke test.  The app installer inherits this engine archive, so this is the
place to enforce it.

Do not remove that external library merely because native ddX appears in a
branch.  Remove it only after the native implementation is merged into the
OpenQP revision selected for the engine, and after the same PCM smoke test
passes using that revision.  This requirement was recorded after a local PCM
job failed with `PCM (ddX) cavity build failed: OpenQP was built without
OQP_ENABLE_DDX`.

## Open problem 3 — report the patch.exe crash upstream

The Windows engine failed twice in the same place, so not a flake:

    DFT-D4 source matches neither side of patch: mctc-lib-0.4.2-disable-tests
      Assertation failed!
      Program: C:\Strawberry\c\bin\patch.exe
      File: ...\patch-2.5.9-src\patch.c, Line 354

The patch is fine. Strawberry Perl's GNU patch 2.5.9 *crashes* on it, and
OpenQP reads the crash as the source not matching — the one thing it is
designed to refuse. The runner image also carries Git for Windows' patch
2.7.6, and which one `find_program(patch)` returns depends on PATH, which is
why the same commit built at 07:00 and not at 08:06.

Our workaround, in `windows-bundle.yml`, names the working one explicitly via
`OQP_DFTD4_PATCH_EXECUTABLE` and runs it once before the hour-long build
starts. OpenQP should pin this itself, and its error message should not point
at its own patch file when the tool is what died. Not yet reported.

## Recently fixed, for context

**Basis set names** (`2f07579`). The first entry in the basis dropdown — the
default every job started from — was `6-31g(d)`, a name the Basis Set Exchange
does not have, so a default submit always died at basis setup with
`KeyError: 'Basis set 6-31g(d) does not exist'`. Because that lands after
pyoqp's settings dump, the job read as a calculation that broke rather than as
a name that never existed. Verified against `basis_set_exchange` directly:
`6-31g*`, `6-31g**`, `6-31+g*`, `6-311g*`, `6-311g**`, `6-311+g**` all
resolve, while `6-31g(d)`, `6-311g(d)` and `6-311+g(d,p)` all raise —
`6-31g(d,p)` and `6-311g(d,p)` happen to exist as separate entries, which is
why only some parenthetical names failed. The dropdown now uses the asterisk
spelling and `normalizeBasis()` in `frontend/src/main.ts` translates a name
typed into the custom box.

**Windows engine build** (`43b694e`) — the patch.exe pin described above.

**PowerShell exit-code masking** — a failing `pip install` in the Windows
engine build was reported green, and surfaced later as a missing DLL. Every
`pwsh` block that runs native commands now sets `$ErrorActionPreference` and
`$PSNativeCommandUseErrorActionPreference`, and the build asserts `import oqp`
before going on.

**Release tag fallback** — `inputs.release_tag || github.ref_name` falls back
to the *branch* on a dispatch with no tag, so build-only runs tried to publish
to a release called "main". The tag is now computed once, in a `context` job
every other job reads.

## Map

| Path | What it holds |
| --- | --- |
| `backend/oqp_studio/main.py` | API, import/analysis endpoints, `/api/workspace`, startup traces |
| `backend/oqp_studio/jobs.py` | Job execution, `JOBS_ROOT`, `adopt`/`rebase`/`set_root` |
| `backend/oqp_studio/workspace.py` | Results-folder choice, write probe, fallback chain |
| `backend/oqp_studio/engine.py` | Bundled-engine discovery, macOS quarantine removal |
| `backend/oqp_studio/server_main.py` | Sidecar entry point, `backend.log` |
| `backend/build_binary.py` | PyInstaller invocation (`--onefile` today) |
| `shell/src-tauri/src/main.rs` | Splash, port wait (180 s), failure screen, sidecar spawn |
| `frontend/src/main.ts` | Input generation, opening and auto-displaying results |
| `frontend/index.html` | The form, including the basis dropdown |
| `.github/workflows/release.yml` | Orchestrates everything; `context` computes the tag once |
| `.github/workflows/installers.yml` | One installer recipe, built slim and with-engine |
| `.github/workflows/{linux,macos,windows}-bundle.yml` | Engine bundles |
| `.github/workflows/prune-releases.yml` | Deletes every release but one |

## Working on it

The loop to live in is the browser one: `uvicorn oqp_studio.main:app --reload
--port 8814` alongside `npm run dev`, open <http://localhost:5173>. Vite
proxies `/api`, so an edit to either side is live and neither Rust nor a
frozen backend is involved.

Only shell behaviour needs `cargo tauri dev`, and that needs a frozen backend
staged in `shell/src-tauri/binaries/` first — the directory ships empty, and
`externalBin` will not let the dev build start without it. The README's
"Working on the desktop shell" has the one-time build, and the note that the
shell reuses a same-version backend already on 8814, which is what keeps
Python edits live inside the window.

## Checks

    cd backend  && python -m pytest tests/ -q      # 18 tests
    cd frontend && npm ci && npx tsc --noEmit && npm run build

## Loose ends

- Thirteen tags remain from the deleted releases. Removing them was never
  decided either way; `prune-releases.yml` takes `delete_tags`.
- The `-with-engine` AppImage is not built: linuxdeploy dies on the 700 MB
  engine tree, so the Linux all-in-one is a `.deb` only. Slim Linux still gets
  both.
