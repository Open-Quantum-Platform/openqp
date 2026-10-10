# OQP Studio

Cross-platform (Windows / macOS / Linux) desktop GUI for the
[Open Quantum Platform](https://github.com/Open-Quantum-Platform/openqp):
build molecules, prepare validated `.oqp` inputs, run OpenQP locally or
remotely, and visualize results (geometries, molecular orbitals, NTOs,
Dyson orbitals, vibrations, spectra) with publication-quality graphics.

Design document: [OQP Studio — Design Proposal](https://github.com/Open-Quantum-Platform/openqp-docs/blob/main/docs/developers/oqp-studio-proposal.md)

## Repository layout

```
frontend/   TypeScript + Vite UI; Mol*-based 3D viewer, Ketcher sketcher,
            input forms, job monitor (also deployable as a website)
backend/    Python FastAPI local server; runs jobs through local or bundled
            OpenQP commands, execution adapters (WSL / SSH), Molden→cube grid engine
shell/      Tauri 2 desktop shell; produces MSI (Windows), DMG (macOS),
            AppImage/deb (Linux) installers
docs/       architecture notes and development guides
```

## Architecture (Phase 0)

```
┌────────────────────────────── OQP Studio ────────────────────────────────┐
│  Desktop shell: Tauri 2                                                  │
│  ┌──────────────────────────┐      ┌───────────────────────────────────┐ │
│  │ Frontend (TypeScript)    │ HTTP │ Local backend (Python, FastAPI)   │ │
│  │  • Builder (Ketcher 2D,  │◄────►│  • RDKit: SMILES→3D, MMFF pre-opt │ │
│  │    3D editor)            │  WS  │  • OpenQP: local or bundled       │ │
│  │  • Template & DB browser │      │  • Job queue                      │ │
│  │  • Input form/editor     │      │  • Grid engine: molden→MO cubes   │ │
│  │  • Mol*-based 3D viewer  │      │  • Execution adapters:            │ │
│  │  • Spectra/plots         │      │    local · bundled · WSL · SSH    │ │
│  └──────────────────────────┘      └───────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────────────┘
```

The compute engine is native on all three platforms, Windows included. Two
installers are published for each: the plain one downloads the engine when you
first ask for it, and the `-with-engine` one carries it, so a single download
installs an application that computes with no network afterwards. The WSL and
SSH adapters remain, for running against a cluster or an OpenQP you installed
yourself.

Bond-distance scans are available from the Workflow menu. A rigid scan runs
single-point calculations at the requested distances; a relaxed scan uses
OpenQP's native `freeze=distance(i,j)` constraint while optimizing every point.
Analysis plots relative energies in kcal/mol, and selecting a point opens that
calculation and structure. Optimization, IRC, and NEB histories use the same
interactive reaction-path plot when their step energies are present.

Analysis can compare any two completed projects without merging their result
identities. It reports the current-minus-reference energy, aligned Cartesian
structure RMSD, dipole-magnitude change, and matched excited-state energy
changes. Atomic results exported by OpenQP are also available as labeled 3D
maps: Mulliken, Lowdin, and RESP partial charges, plus coupled NMR shielding.

Workflow names and descriptions are searchable. End-to-end recipes extend a
workflow with Method, Execution, Analysis, and Art settings and can be exported
as `oqp-studio-recipe/v1` JSON, moved to another computer, and imported there.
Built-in recipes cover vertical absorption with NTO analysis, excited-state
relaxation with emission/ESA analysis, EKT IP/EA with Dyson orbitals,
vibrational characterization, and NMR shielding maps. NICS is listed but
disabled because OpenQP currently evaluates shielding only at real nuclei and
does not yet accept the ghost/probe centers required for a NICS value.

A recipe may contain optional Python for post-calculation analysis. Imported
Python is never trusted automatically: its source must be reviewed and trusted
on that computer before it can run. Trust is stored separately from exported
recipe JSON. Approved code runs only after a successful OpenQP calculation,
with a timeout and stdout/stderr captured in `postprocess.log`. This is a trust
boundary, not a sandbox; approved Python can access local files and the network.
See [docs/recipes.md](docs/recipes.md) for the portable format and trust model.

Gaussian cube results can be displayed directly or combined point by point as
a sum or primary-minus-secondary difference when their grids and atom headers
match. The surface controls expose the contour value and positive, negative,
or both signs; incompatible grids are rejected rather than resampled silently.

The Builder's symmetry analysis reports a likely molecular point group, the
accepted operation count, maximum coordinate residual, and symmetry-equivalent
atom groups at the displayed tolerance. Principal-axis alignment is a separate
explicit action, so analysis never changes the coordinates unless requested.

Until OpenQP's native continuum solver is merged, the bundled engine also
contains the external ddX library required for PCM calculations.

## Install with pip

```bash
pip install oqp-studio        # ships the server, UI, and results viewer
oqp-studio                    # opens in its own window
```

Add the extras you need: `pip install "oqp-studio[desktop,chem]"` for the
native window (pywebview) and 2D-sketch-to-3D conversion (RDKit).

## Run it as a desktop app (no browser)

After building the frontend once (see below), the backend can open its own
The published desktop application is the Tauri shell in `shell/`. It bundles
the frontend and starts the Python backend as a sidecar over standard input
and output; the installed application does not open a local HTTP port.

## Development quick start

Backend API development (requires Python ≥ 3.10; OpenQP optional — mock mode without it):

```bash
cd backend
python -m venv .venv && . .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
uvicorn oqp_studio.main:app --reload --port 8814
```

Frontend (requires Node ≥ 24.14.1, matching Ketcher's supported runtime):

```bash
cd frontend
npm install
npm run dev
```

Then open <http://localhost:5173>. This optional browser-only loop uses the
Vite proxy; it is not used by the desktop application.

### Working on the desktop shell

The desktop shell loads the built frontend directly and sends `/api` requests
to its Python sidecar over stdio. It has one prerequisite: `externalBin` in
`tauri.conf.json` declares the frozen backend, so stage one before running it.
Build it once:

```bash
cd frontend && npm run build          # build_binary.py embeds frontend/dist
cd ../backend && pip install -e . pyinstaller && python build_binary.py
mkdir -p ../shell/src-tauri/binaries
triple=$(rustc -Vv | sed -n 's/^host: //p')
if [ "$(uname)" = "Darwin" ]; then
  cp dist/oqp-studio-backend/oqp-studio-backend \
     "../shell/src-tauri/binaries/oqp-studio-backend-$triple"
  cp -R dist/oqp-studio-backend/_internal \
     ../shell/src-tauri/binaries/oqp-studio-backend-runtime
else
  cp dist/oqp-studio-backend \
     "../shell/src-tauri/binaries/oqp-studio-backend-$triple"
fi
```

```bash
cd shell/src-tauri
cargo run
```

For a frontend edit, rerun `npm run build` before `cargo run`. For a backend
edit, rerun `python build_binary.py`, copy the sidecar again, then rerun
`cargo run`. This checks the same no-network architecture as the packaged
application without creating a release.

### Measuring startup

The frozen backend prints a timing trace to stderr, which nothing captures
when the app is launched from Finder (`backend.log` is only written when
startup *raises*). Run it directly instead — this also isolates it from Tauri:

```bash
time "/Applications/OQP Studio.app/Contents/MacOS/oqp-studio-backend" --port 8899
```

The delay before the first `startup …s` line is time spent before Python runs
at all; the gaps between lines are import and setup cost. See
[docs/handoff.md](docs/handoff.md) for what each answer implies.

## Engine upstream and releasing

Studio's engine tracks the external-publication gateway
[`open-quantum-platform/openqp`](https://qchemlab.knu.ac.kr/open-quantum-platform/openqp),
not the private development superset. The build resolves gateway `main` to one
exact commit shared by all platform jobs. See
[Studio releases and engine compatibility](docs/studio-releases.md) for the
checkout command, source provenance, and proposed independent release channels.

Studio source remains private. Existing public installers are attached to
[OpenQP v1.3.1](https://github.com/Open-Quantum-Platform/openqp/releases/tag/v1.3.1).
The historical GitHub release/PyPI workflows are disabled and are not the new
publication path. A Studio tag alone does not authorize a public source release.

## Code signing

Installers are currently unsigned, so macOS and Windows show a first-run
warning. On macOS, clear the quarantine flag once after installing:

```bash
xattr -cr "/Applications/OQP Studio.app"
```

To skip the warning without any certificate, install the `.app.tar.gz` build
from the terminal instead — `curl` does not set the quarantine flag:

```bash
curl -L -o oqp-studio.tar.gz <tarball URL>
tar xzf oqp-studio.tar.gz -C /Applications
```

The release workflow signs and notarizes automatically once certificates are
added as repository secrets. Accredited educational institutions can obtain
the Apple Developer Program at no cost through Apple's fee waiver — see
[docs/code-signing.md](docs/code-signing.md).

## Status

Working end to end: build a molecule (2D sketcher, PubChem, samples, or raw
coordinates), pick a workflow, generate the `.oqp` input, run it through a
local OpenQP (native or WSL), and inspect the results —
orbitals, normal modes, geometries — rendered with Mol*. Installers are built
for Windows, macOS (Apple Silicon and Intel), and Linux.

Studio accepts and generates OpenQP's canonical concise `.oqp` input only; it
is not a Studio-specific format. Legacy sectioned input is not supported.
Geometry workflows use OpenQP's native optimizer automatically. Studio emits
native controls through the canonical `oqp(...)` call and never emits a legacy
`lib` option. For NACME, provide both the current structure and the
previous-geometry `.xyz` file, which Studio writes as `geom2` in the input.

Remaining: SSH/SLURM submission, style presets and high-resolution export,
code signing (see above). See the design proposal for the full roadmap.

Picking the work up: [docs/handoff.md](docs/handoff.md) says where the project
stands, which problems are open and what has already been ruled out on each.

## License

OQP Studio is source-available under the same dual-licensing model as OpenQP:
the [OpenQP Research License 1.0](LICENSE) permits qualified non-commercial
research use, while commercial use requires a separate written license. This
research license is not an open-source license. Third-party components remain
under their respective licenses.
