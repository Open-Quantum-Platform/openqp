from __future__ import annotations

import os
import re
import ssl
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from collections import OrderedDict
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import __version__, engine, environment, host, network, recipes, releases, scans
from .jobs import JobInfo, JobRequest, JobStatus, manager
from .runners import available_runners

# How long each piece of starting up took, on stderr. The shell gives the
# backend a fixed number of seconds to open its port and this took fourteen to
# twenty on a healthy machine, which is close enough to the limit that a
# slower one fails outright. Guessing at which part is slow is what produced
# the wrong fix last time, so it is measured.
_STARTED = time.monotonic()


def _trace(step: str) -> None:
    print(f"startup {time.monotonic() - _STARTED:6.2f}s  {step}", file=sys.stderr, flush=True)


def _warm_up_rdkit() -> None:
    """Import RDKit ahead of the first request, but not ahead of the server.

    The import costs a second or more in a frozen build, and it would
    otherwise be charged to whoever first converts a 2D sketch to 3D. It runs
    in a thread, which does not make it free: loading that many C extensions
    holds the GIL in long stretches, so doing it now competes with the work
    that opens the port. It waits until the server has had time to start.
    """

    def load() -> None:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import AllChem  # noqa: F401
        except ImportError:
            return
        _trace("rdkit ready")

    timer = threading.Timer(15.0, load)
    timer.daemon = True
    timer.start()


_trace("imports")
environment.enrich_path()
_trace("PATH")
network.activate()
_trace("TLS")
_warm_up_rdkit()

def _warm_runner_versions() -> None:
    """Read engine versions once, after the server is ready for requests."""
    local_openqp = environment.locate(os.environ.get("OQP_COMMAND", "openqp"))
    engine.version(local_openqp)
    engine.version(engine.bundled_or_downloaded())


@asynccontextmanager
async def _lifespan(_app: FastAPI):
    # Version probing may start a Python interpreter for a pip-installed
    # command. It must never add latency to the sidecar becoming reachable.
    threading.Thread(target=_warm_runner_versions, daemon=True).start()
    yield


MAX_SYMMETRY_REQUEST_BYTES = 1_100_000


class _SymmetryBodyLimit:
    """Reject oversized symmetry JSON before FastAPI buffers or binds it."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope.get("path") != "/api/symmetry":
            await self.app(scope, receive, send)
            return
        headers = dict(scope.get("headers", []))
        content_length = headers.get(b"content-length")
        if content_length is not None:
            try:
                if int(content_length) > MAX_SYMMETRY_REQUEST_BYTES:
                    await JSONResponse(
                        {"detail": "symmetry request body is too large"}, status_code=413,
                    )(scope, receive, send)
                    return
            except ValueError:
                pass

        chunks: list[bytes] = []
        total = 0
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                return
            chunk = message.get("body", b"")
            total += len(chunk)
            if total > MAX_SYMMETRY_REQUEST_BYTES:
                await JSONResponse(
                    {"detail": "symmetry request body is too large"}, status_code=413,
                )(scope, receive, send)
                return
            chunks.append(chunk)
            if not message.get("more_body", False):
                break

        body = b"".join(chunks)
        delivered = False

        async def replay():
            nonlocal delivered
            if delivered:
                return {"type": "http.request", "body": b"", "more_body": False}
            delivered = True
            return {"type": "http.request", "body": body, "more_body": False}

        await self.app(scope, replay, send)


app = FastAPI(title="OQP Studio backend", version=__version__, lifespan=_lifespan)
app.add_middleware(_SymmetryBodyLimit)
_trace("app built")

# FastAPI reads this as the file field; kept module level so it is not a
# function call in a default argument.
_UPLOAD = File(...)
_UPLOADS = File(...)
_IMPORT_NAME = Form("")


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "version": __version__, "build": releases.build_identity()}


@app.get("/api/runners")
def runners() -> dict:
    return available_runners()


@app.get("/api/runners/detail")
def runner_detail() -> dict:
    """Which runners work and, for the native one, the binary that was found."""
    from .jobs import JOBS_ROOT

    local_openqp = environment.locate(os.environ.get("OQP_COMMAND", "openqp"))
    bundled_openqp = engine.bundled_or_downloaded()
    return {
        "available": available_runners(),
        "jobs_root": str(JOBS_ROOT),
        **environment.describe(),
        "openqp": local_openqp,
        "bundled_openqp": bundled_openqp,
        "identities": {"local": engine.identity(local_openqp),
                       "bundled": engine.identity(bundled_openqp)},
        "versions": {
            "local": engine.version(local_openqp),
            "bundled": engine.version(bundled_openqp),
        },
    }


class ResourceRequest(BaseModel):
    input_text: str
    threads: int = 1


class SymmetryRequest(BaseModel):
    xyz: str
    tolerance: float = 0.05


@app.get("/api/host")
def host_status() -> dict:
    return host.snapshot()


@app.post("/api/host/admission")
def host_admission(req: ResourceRequest) -> dict:
    return host.admission(req.input_text, req.threads)


@app.post("/api/symmetry")
def molecular_symmetry(req: SymmetryRequest) -> dict:
    """Likely point group, accepted operations, and a principal-axis frame."""
    from . import symmetry

    if not 1.0e-4 <= req.tolerance <= 0.5:
        raise HTTPException(status_code=422, detail="tolerance must be between 0.0001 and 0.5 A")
    try:
        return symmetry.analyze(req.xyz, req.tolerance)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/api/jobs")
def submit_job(req: JobRequest) -> JobInfo:
    try:
        return manager.submit(req)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except OSError as exc:
        # A bare 500 says nothing; the thing that goes wrong here is the job
        # directory being somewhere the user cannot write, so name it.
        from .jobs import JOBS_ROOT

        raise HTTPException(
            status_code=500,
            detail=f"could not create a job directory under {JOBS_ROOT}: {exc}",
        ) from exc


@app.get("/api/recipes")
def list_recipes() -> list[recipes.RecipeRecord]:
    return recipes.list_records()


@app.get("/api/recipes/{recipe_id}")
def get_recipe(recipe_id: str) -> recipes.RecipeRecord:
    record = recipes.get(recipe_id)
    if record is None:
        raise HTTPException(status_code=404, detail="recipe not found")
    return record


@app.post("/api/recipes")
def save_recipe(req: recipes.SaveRecipeRequest) -> recipes.RecipeRecord:
    try:
        return recipes.save(req)
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.delete("/api/recipes/{recipe_id}")
def delete_recipe(recipe_id: str) -> dict:
    try:
        recipes.delete(recipe_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="recipe not found") from exc
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"deleted": recipe_id}


@app.post("/api/scans")
def submit_bond_scan(req: scans.BondScanRequest) -> dict:
    """Create a bond-distance scan whose calculation points run serially."""
    try:
        group_id, requests, values = scans.build(req)
        jobs = manager.submit_batch(
            requests, group_id=group_id, values=values, unit="A",
            state=scans.target_state(req.input_text),
        )
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"group_id": group_id, "jobs": jobs}


@app.get("/api/scans/{group_id}")
def bond_scan(group_id: str) -> dict:
    """Current energies and statuses for one persisted scan group."""
    points = []
    for info in manager.list():
        if info.group_id != group_id:
            continue
        scan_energy = _scan_point_energy(info) if info.status == "done" else None
        points.append({
            "job_id": info.id,
            "name": info.name,
            "status": info.status,
            "value": info.scan_value,
            "unit": info.scan_unit,
            "energy": scan_energy,
        })
    points.sort(key=lambda point: float(point["value"] or 0.0))
    if not points:
        raise HTTPException(status_code=404, detail="scan not found")
    return {"group_id": group_id, "points": points}


@app.get("/api/jobs")
def list_jobs() -> list[JobInfo]:
    return manager.list()


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> JobInfo:
    info = manager.get(job_id)
    if info is None:
        raise HTTPException(status_code=404, detail="job not found")
    return info


@app.post("/api/jobs/{job_id}/restart")
def restart_job(job_id: str) -> dict:
    try:
        return manager.restart_input(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str) -> JobInfo:
    try:
        manager.cancel(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    info = manager.get(job_id)
    assert info is not None
    return info


@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str) -> dict:
    """Permanently remove a completed project and all of its files."""
    try:
        manager.delete(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"could not delete project: {exc}") from exc
    _summary_cache.pop(job_id, None)
    _scan_energy_cache.pop(job_id, None)
    return {"deleted": job_id}


@app.get("/api/jobs/{job_id}/log")
def job_log(job_id: str) -> dict:
    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    return {"log": manager.log_tail(job_id)}


@app.get("/api/jobs/{job_id}/files")
def job_files(job_id: str) -> list[dict]:
    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    return manager.files(job_id)


@app.get("/api/jobs/{job_id}/files/{name}")
def job_file(job_id: str, name: str) -> FileResponse:
    path = manager.file_path(job_id, name)
    if path is None:
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(path, media_type="text/plain", filename=name)


@app.get("/api/jobs/{job_id}/cube-combine")
def combine_cube_files(job_id: str, left: str, right: str,
                       operation: str = "difference") -> PlainTextResponse:
    """Pointwise sum or difference of two compatible Gaussian cube files."""
    from . import cube

    left_path = manager.file_path(job_id, left)
    right_path = manager.file_path(job_id, right)
    if left_path is None or right_path is None:
        raise HTTPException(status_code=404, detail="cube file not found")
    if left_path.suffix.lower() not in {".cube", ".cub"} or right_path.suffix.lower() not in {
        ".cube", ".cub",
    }:
        raise HTTPException(status_code=422, detail="cube arithmetic requires .cube or .cub files")
    try:
        maximum_bytes = 64 * 1024 * 1024
        if left_path.stat().st_size + right_path.stat().st_size > maximum_bytes:
            raise HTTPException(
                status_code=413,
                detail=f"cube arithmetic is limited to {maximum_bytes // (1024 * 1024)} MiB combined",
            )
        result = cube.combine(
            left_path.read_text(errors="replace"),
            right_path.read_text(errors="replace"),
            operation,
        )
    except HTTPException:
        raise
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return PlainTextResponse(result, media_type="text/plain")


@app.get("/api/jobs/{job_id}/cube-geometry")
def cube_file_geometry(job_id: str, name: str) -> dict:
    """Atomic geometry embedded in a Gaussian cube, converted to angstrom."""
    from . import cube

    path = manager.file_path(job_id, name)
    if path is None:
        raise HTTPException(status_code=404, detail="cube file not found")
    if path.suffix.lower() not in {".cube", ".cub"}:
        raise HTTPException(status_code=422, detail="cube geometry requires a .cube or .cub file")
    try:
        with path.open(errors="replace") as stream:
            return {"xyz": cube.geometry_xyz(cube.parse_header(stream))}
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


class WorkspaceRequest(BaseModel):
    jobs_dir: str


@app.get("/api/workspace")
def workspace_status() -> dict:
    """Where results are written, and where the user asked for them."""
    from . import jobs as jobs_module
    from . import workspace

    return workspace.status(jobs_module.JOBS_ROOT)


@app.post("/api/workspace")
def set_workspace(req: WorkspaceRequest) -> dict:
    """Choose where results are written from now on."""
    from pathlib import Path

    from . import jobs as jobs_module
    from . import workspace

    if manager.busy():
        raise HTTPException(
            status_code=409,
            detail="a job is still running; let it finish before moving the folder")
    wanted = req.jobs_dir.strip()
    if wanted:
        try:
            target = workspace.usable(Path(wanted).expanduser())
        except OSError as exc:
            raise HTTPException(status_code=400,
                                detail=f"cannot write to {wanted}: {exc}") from exc
        workspace.save(str(target))
        jobs_module.set_root(target.resolve())
    else:
        # An empty box means "go back to choosing for me".
        workspace.save("")
        jobs_module.set_root(workspace.resolve())
    return workspace.status(jobs_module.JOBS_ROOT)


# Result files worth adopting when a whole folder is imported. Everything
# else in a run directory -- scratch, restart files, submission scripts --
# would only pad the file list.
RESULT_SUFFIXES = (
    ".log", ".out", ".json", ".molden", ".xyz", ".trj", ".cube", ".cub",
    ".inp", ".oqp", ".dat", ".txt",
)

# Ignore anything larger than this on a folder import: a stray checkpoint is
# not analysable and copying it would only cost the user disk.
MAX_IMPORT_BYTES = 512 * 1024 * 1024
MAX_COMPARISON_TEXT_BYTES = 16 * 1024 * 1024
MAX_COMPARISON_LINE = 64 * 1024
MAX_COMPARISON_ATOMS = 100_000


class ImportPathRequest(BaseModel):
    path: str
    name: str = ""


def _import_candidates(root: Path) -> list[Path]:
    if root.is_file():
        return [root]
    return sorted(
        path for path in root.iterdir()
        if path.is_file() and path.suffix.lower() in RESULT_SUFFIXES
        and path.stat().st_size <= MAX_IMPORT_BYTES
    )


@app.post("/api/jobs/import")
async def import_results(files: list[UploadFile] = _UPLOADS,
                         name: str = _IMPORT_NAME) -> JobInfo:
    """Adopt result files produced outside this app.

    A run made on a cluster, by the standalone command-line engine, or by an
    earlier session becomes a job here, so every analysis the app offers --
    summary, spectra, orbitals, normal modes, property maps -- applies to it
    exactly as it does to a run this app started.
    """
    payload: list[tuple[str, bytes]] = []
    for upload in files:
        data = await upload.read()
        if len(data) > MAX_IMPORT_BYTES:
            raise HTTPException(status_code=413,
                                detail=f"{upload.filename} is too large to import")
        payload.append((upload.filename or "file", data))
    try:
        return manager.adopt(name or "imported", payload)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/jobs/import-path")
def import_results_from_path(req: ImportPathRequest) -> JobInfo:
    """Adopt a result file, or a whole run directory, already on this machine.

    The desktop app and the backend share a filesystem, so pointing at a
    directory beats uploading it: a molden file of a few hundred megabytes is
    copied locally instead of pushed through the browser.
    """
    root = Path(req.path).expanduser()
    if not root.exists():
        raise HTTPException(status_code=404, detail=f"no such path: {root}")
    candidates = _import_candidates(root)
    if not candidates:
        raise HTTPException(status_code=400,
                            detail="no result files were found in that folder")
    try:
        return manager.adopt(req.name or root.stem or "imported",
                             [(path.name, path.read_bytes()) for path in candidates])
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


PUBCHEM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF?record_type=3d"
)


def _parse_sdf_atoms(sdf: str) -> list[list]:
    """Atom list [[symbol, x, y, z], ...] from a V2000 SDF counts+atom block."""
    lines = sdf.splitlines()
    if len(lines) < 4:
        raise ValueError("SDF too short")
    natoms = int(lines[3][0:3])
    atoms = []
    for line in lines[4 : 4 + natoms]:
        x, y, z, symbol = float(line[0:10]), float(line[10:20]), float(line[20:30]), line[31:34]
        atoms.append([symbol.strip(), x, y, z])
    return atoms


@app.get("/api/pubchem/{name}")
def pubchem_lookup(name: str) -> dict:
    """Fetch a 3D structure from PubChem by compound name."""
    url = PUBCHEM_URL.format(name=urllib.parse.quote(name))
    try:
        with urllib.request.urlopen(url, timeout=15,
                                    context=network.context()) as response:
            sdf = response.read().decode()
    except urllib.error.HTTPError as exc:
        detail = f"PubChem has no 3D record for '{name}'" if exc.code == 404 else str(exc)
        raise HTTPException(status_code=404, detail=detail) from exc
    except urllib.error.URLError as exc:
        # Match on the message as well as the type: the certificate error can
        # arrive wrapped differently depending on how the request failed.
        certificate_problem = (isinstance(exc.reason, ssl.SSLCertVerificationError)
                               or "CERTIFICATE_VERIFY" in str(exc))
        if certificate_problem:
            raise HTTPException(
                status_code=502,
                detail=(
                    "PubChem's certificate could not be verified — a network that "
                    "inspects TLS traffic is the usual cause. Open File \u2192 "
                    "Network settings to point OQP Studio at your network's root "
                    "certificate. Sketching or typing coordinates needs no network."
                ),
            ) from exc
        raise HTTPException(status_code=502, detail=f"PubChem unreachable: {exc}") from exc
    except OSError as exc:
        raise HTTPException(status_code=502, detail=f"PubChem unreachable: {exc}") from exc
    try:
        atoms = _parse_sdf_atoms(sdf)
    except (ValueError, IndexError) as exc:
        raise HTTPException(status_code=502, detail="could not parse PubChem SDF") from exc
    return {"name": name, "atoms": atoms}


RCSB_URL = "https://files.rcsb.org/download/{code}.pdb"


@app.get("/api/pdb/{code}")
def rcsb_lookup(code: str) -> dict:
    """Fetch an experimental structure from the RCSB Protein Data Bank.

    PubChem holds small molecules; proteins, nucleic acids and complexes come
    from the PDB, addressed by their four-character entry ID.
    """
    entry = code.strip().upper()
    if not re.fullmatch(r"[0-9A-Z]{4}", entry):
        raise HTTPException(
            status_code=400,
            detail=f"'{code}' is not a PDB ID — those are four characters, such as 1CRN",
        )
    try:
        with urllib.request.urlopen(RCSB_URL.format(code=entry), timeout=30,
                                    context=network.context()) as response:
            text = response.read().decode(errors="replace")
    except urllib.error.HTTPError as exc:
        detail = (f"the PDB has no entry {entry}" if exc.code == 404 else str(exc))
        raise HTTPException(status_code=404, detail=detail) from exc
    except urllib.error.URLError as exc:
        if isinstance(exc.reason, ssl.SSLCertVerificationError) or "CERTIFICATE_VERIFY" in str(exc):
            raise HTTPException(
                status_code=502,
                detail=(
                    "The PDB's certificate could not be verified — open File \u2192 "
                    "Network settings to trust your network's root certificate."
                ),
            ) from exc
        raise HTTPException(status_code=502, detail=f"RCSB unreachable: {exc}") from exc
    return {"code": entry, "pdb": text}


_molden_cache: OrderedDict[str, object] = OrderedDict()


def _load_molden(job_id: str, name: str):
    from . import molden

    key = f"{job_id}/{name}"
    if key not in _molden_cache:
        path = manager.file_path(job_id, name)
        if path is None:
            raise HTTPException(status_code=404, detail="file not found")
        _molden_cache[key] = molden.parse_molden(path.read_text(errors="replace"))
        while len(_molden_cache) > 8:
            _molden_cache.popitem(last=False)
    return _molden_cache[key]


def _dyson_source_state(job_id: str) -> str | None:
    """Return the physical neutral state used for an IP/EA calculation."""
    pattern = re.compile(r"PyOQP Physical target state\(s\):\s*(S\d+)", re.IGNORECASE)
    for entry in manager.files(job_id):
        name = entry["name"]
        if not name.lower().endswith((".log", ".out")):
            continue
        path = manager.file_path(job_id, name)
        if path is None:
            continue
        match = pattern.search(path.read_text(errors="replace"))
        if match:
            return match.group(1).upper()
    return None


@app.get("/api/jobs/{job_id}/molden/{name}/orbitals")
def molden_orbitals(job_id: str, name: str) -> dict:
    data = _load_molden(job_id, name)
    if not data.supported:
        raise HTTPException(
            status_code=422,
            detail=f"spherical shells not yet supported: {', '.join(data.unsupported)}",
        )
    dyson_source_state = _dyson_source_state(job_id)

    def orbital_info(orbital):
        """Expose Dyson metadata without treating its strength as occupancy."""
        result = {
            "index": orbital.index,
            "energy": orbital.energy,
            "spin": orbital.spin,
            "occupancy": orbital.occupancy,
            "symmetry": orbital.symmetry,
            "kind": "scf",
        }
        match = re.fullmatch(
            r"Dyson-(IP|EA)-state-(\d+)", orbital.symmetry or "", re.IGNORECASE,
        )
        if match:
            strength = orbital.occupancy
            result.update({
                "kind": "dyson",
                "dyson_kind": match.group(1).upper(),
                "state_index": int(match.group(2)),
                "source_state": dyson_source_state,
                "strength": strength,
                # OpenQP writes the Dyson pole strength in Molden's Occup field.
                "occupation": 2.0 * strength if strength is not None else None,
            })
        return result

    return {
        "atoms": len(data.atoms),
        "orbitals": [orbital_info(orbital) for orbital in data.orbitals],
    }


@app.get("/api/jobs/{job_id}/molden/{name}/geom.xyz")
def molden_geometry(job_id: str, name: str) -> PlainTextResponse:
    from . import molden

    return PlainTextResponse(molden.atoms_to_xyz(_load_molden(job_id, name), name))


@app.get("/api/jobs/{job_id}/molden/{name}/cube")
def molden_cube(job_id: str, name: str, mo: int) -> PlainTextResponse:
    from . import molden

    data = _load_molden(job_id, name)
    if not data.supported:
        raise HTTPException(status_code=422, detail="spherical shells not yet supported")
    if not 1 <= mo <= len(data.orbitals):
        raise HTTPException(status_code=404, detail="orbital index out of range")
    return PlainTextResponse(molden.orbital_cube(data, mo))


_vib_cache: OrderedDict[str, object] = OrderedDict()


def _load_vibrations(job_id: str, name: str):
    from . import molden

    key = f"{job_id}/{name}"
    if key not in _vib_cache:
        path = manager.file_path(job_id, name)
        if path is None:
            raise HTTPException(status_code=404, detail="file not found")
        _vib_cache[key] = molden.parse_vibrations(path.read_text(errors="replace"))
        while len(_vib_cache) > 8:
            _vib_cache.popitem(last=False)
    return _vib_cache[key]


@app.get("/api/jobs/{job_id}/molden/{name}/modes")
def molden_modes(job_id: str, name: str) -> dict:
    vib = _load_vibrations(job_id, name)
    return {
        "atoms": len(vib.atoms),
        "modes": [
            {
                "index": m.index,
                "frequency": m.frequency,
                "intensity": vib.intensities[i] if i < len(vib.intensities) else None,
            }
            for i, m in enumerate(vib.modes)
        ],
    }


@app.get("/api/jobs/{job_id}/molden/{name}/mode.xyz")
def molden_mode_trajectory(job_id: str, name: str, mode: int,
                           amplitude: float = 0.6) -> PlainTextResponse:
    from . import molden

    vib = _load_vibrations(job_id, name)
    if not 1 <= mode <= len(vib.modes):
        raise HTTPException(status_code=404, detail="mode index out of range")
    return PlainTextResponse(
        molden.mode_trajectory(vib, mode, amplitude=max(0.05, min(amplitude, 2.0)))
    )


@app.get("/api/jobs/{job_id}/molden/{name}/mode")
def molden_mode_vectors(job_id: str, name: str, mode: int) -> dict:
    """Equilibrium coordinates and one normal-mode displacement in Angstrom."""
    from . import molden

    vib = _load_vibrations(job_id, name)
    if not 1 <= mode <= len(vib.modes):
        raise HTTPException(status_code=404, detail="mode index out of range")
    normal_mode = vib.modes[mode - 1]
    displacements = normal_mode.displacements * molden.BOHR_TO_ANGSTROM
    return {
        "index": normal_mode.index,
        "frequency": normal_mode.frequency,
        "atoms": [
            {
                "element": element,
                "position": [x, y, z],
                "displacement": displacement.tolist(),
            }
            for (element, x, y, z), displacement in zip(vib.atoms, displacements)
        ],
    }


_summary_cache: OrderedDict[str, dict] = OrderedDict()
_scan_energy_cache: dict[str, float | None] = {}


def _job_summary(job_id: str) -> dict:
    from . import analysis

    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    cached = _summary_cache.get(job_id)
    if cached is None:
        paths = [
            path for entry in manager.files(job_id)
            if (path := manager.file_path(job_id, entry["name"])) is not None
        ]
        cached = analysis.summarize(paths)
        _summary_cache[job_id] = cached
        while len(_summary_cache) > 8:
            _summary_cache.popitem(last=False)
    return cached


def _scan_point_energy(info: JobInfo) -> float | None:
    """Read one terminal scan energy once, using the final optimized root."""
    if info.id in _scan_energy_cache:
        return _scan_energy_cache[info.id]
    report = _job_summary(info.id)
    energy = report.get("energy", {})
    value = energy.get("total", energy.get("components", {}).get("total"))
    if info.scan_state is not None:
        from . import analysis

        history = analysis.optimization_history(_job_paths(info.id)).get("steps", [])
        final_states = next(
            (step.get("states", []) for step in reversed(history) if step.get("states")),
            [],
        )
        states = final_states or report.get("states", [])
        state = next((row for row in states if row.get("index") == info.scan_state), None)
        value = state.get("total") if state else None
    _scan_energy_cache[info.id] = value
    return value


def _job_paths(job_id: str) -> list[Path]:
    return [
        path for entry in manager.files(job_id)
        if (path := manager.file_path(job_id, entry["name"])) is not None
    ]


def _summary_energy(summary: dict) -> float | None:
    energy = summary.get("energy", {})
    optimized_state = summary.get("excited_state_optimized")
    final_states = energy.get("final_states", {})
    if optimized_state is not None:
        value = final_states.get(optimized_state, final_states.get(str(optimized_state)))
        if value is not None:
            return float(value)
    value = energy.get("total", energy.get("components", {}).get("total"))
    if value is None:
        value = summary.get("scf", {}).get("energy")
    return float(value) if value is not None else None


def _comparison_frame(paths: list[Path]):
    from . import structure_io

    def rank(path: Path) -> tuple[int, str] | None:
        name = path.name.lower()
        if name in {"opt.xyz", "opt_geom.xyz"} or name.endswith((".namd.trj", ".trj")):
            return (0, name)
        if name.endswith(".xyz") and any(word in name for word in ("result", "final", "optimized")):
            return (1, name)
        if name.endswith((".log", ".out")):
            return (2, name)
        if name.endswith((".molden", ".freq.molden")):
            return (3, name)
        if name.endswith(".xyz"):
            return (4, name)
        if name.endswith(".json") and any(word in name for word in ("result", "output", "hess")):
            return (5, name)
        if name.endswith(".json"):
            return (6, name)
        if name.endswith((
            ".oqp", ".inp", ".pdb", ".ent", ".mol", ".sdf", ".sd", ".mol2",
            ".cdxml", ".cdx", ".smi", ".smiles",
        )):
            return (7, name)
        if name.endswith(".txt"):
            return (8, name)
        if name.endswith((".cube", ".cub")):
            return None
        return (9, name)

    candidates = [(priority, path) for path in paths if (priority := rank(path)) is not None]
    def parse_path(path: Path):
        payload = b"" if path.name.lower().endswith((".namd.trj", ".trj")) else path.read_bytes()
        return structure_io.parse(path.name, payload, path=str(path))

    def usable_frame(structure):
        if structure.frames and structure.frames[-1].atoms:
            return structure.frames[-1]
        return None

    def openqp_text_frame(path: Path):
        last = None
        step = 0

        def bounded_lines(stream):
            while line := stream.readline(MAX_COMPARISON_LINE + 1):
                if len(line) <= MAX_COMPARISON_LINE:
                    yield line
                    continue
                while line and not line.endswith("\n"):
                    line = stream.readline(MAX_COMPARISON_LINE + 1)
                yield ""

        with path.open(errors="replace") as stream:
            rows = iter(bounded_lines(stream))
            for line in rows:
                marker = re.search(r"Geometry Optimization Step\s+(\d+)", line)
                if marker:
                    step = int(marker.group(1))
                if "Cartesian Coordinate in Angstrom" not in line:
                    continue
                for line in rows:
                    if re.match(r"\s*-{5,}", line):
                        break
                atoms = []
                oversized = False
                for line in rows:
                    row = re.match(
                        r"\s*\d+\s+(\d+(?:\.\d+)?)\s+([-+0-9.Ee]+)\s+"
                        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
                        line,
                    )
                    if not row:
                        break
                    if oversized or len(atoms) >= MAX_COMPARISON_ATOMS:
                        atoms = []
                        oversized = True
                        continue
                    atomic_number = int(float(row.group(1)))
                    if 0 < atomic_number < len(structure_io.SYMBOLS):
                        atoms.append(
                            (
                                structure_io.SYMBOLS[atomic_number],
                                float(row.group(2)), float(row.group(3)), float(row.group(4)),
                            )
                        )
                if atoms and not oversized:
                    last = structure_io.Frame(atoms, f"step {step or 1}")
        return last

    def generic_text_frame(path: Path):
        with path.open("rb") as stream:
            data = stream.read(MAX_COMPARISON_TEXT_BYTES + 1)
        if len(data) > MAX_COMPARISON_TEXT_BYTES:
            return None
        text = data.decode(errors="replace")
        for reader in (structure_io.parse_xyz, structure_io.parse_pdb):
            frames = reader(text)
            if frames and frames[-1].atoms:
                return frames[-1]
        return None

    for priority in (0, 1):
        for _rank, path in sorted(item for item in candidates if item[0][0] == priority):
            try:
                structure = parse_path(path)
            except (OSError, ValueError):
                continue
            if frame := usable_frame(structure):
                return frame

    # A .txt file is promoted beside .log/.out only after its contents identify
    # it as an OpenQP log. Do this only after explicit result geometry tiers fail.
    log_candidates = [item for item in candidates if item[0][0] == 2]
    text_candidates = [item for item in candidates if item[0][0] == 8]
    for _rank, path in sorted([*log_candidates, *text_candidates], key=lambda item: item[1].name):
        try:
            frame = openqp_text_frame(path)
        except (OSError, ValueError):
            continue
        if frame:
            return frame

    for priority in range(3, 8):
        for _rank, path in sorted(item for item in candidates if item[0][0] == priority):
            try:
                structure = parse_path(path)
            except (OSError, ValueError):
                continue
            if frame := usable_frame(structure):
                return frame
    for _rank, path in sorted(item for item in candidates if item[0][0] == 8):
        try:
            frame = generic_text_frame(path)
        except OSError:
            continue
        if frame:
            return frame
    for _rank, path in sorted(item for item in candidates if item[0][0] == 9):
        try:
            frame = openqp_text_frame(path)
        except (OSError, ValueError):
            continue
        if frame:
            return frame
        try:
            frame = generic_text_frame(path)
        except OSError:
            continue
        if frame:
            return frame
    return None


def _geometry_comparison(left_paths: list[Path], right_paths: list[Path]) -> dict:
    import numpy as np

    left = _comparison_frame(left_paths)
    right = _comparison_frame(right_paths)
    if left is None or right is None:
        return {"available": False, "reason": "a comparable structure is missing"}
    left_symbols = [atom[0] for atom in left.atoms]
    right_symbols = [atom[0] for atom in right.atoms]
    if left_symbols != right_symbols:
        return {
            "available": False,
            "reason": "atom counts or atom ordering differ",
            "left_atoms": len(left_symbols), "right_atoms": len(right_symbols),
        }
    left_xyz = np.asarray([atom[1:] for atom in left.atoms], dtype=float)
    right_xyz = np.asarray([atom[1:] for atom in right.atoms], dtype=float)
    left_centered = left_xyz - left_xyz.mean(axis=0)
    right_centered = right_xyz - right_xyz.mean(axis=0)
    u, _, vt = np.linalg.svd(left_centered.T @ right_centered)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    difference = left_centered @ rotation - right_centered
    return {
        "available": True,
        "atoms": len(left_symbols),
        "rmsd_angstrom": float(np.sqrt(np.mean(np.sum(difference * difference, axis=1)))),
    }


@app.get("/api/comparison")
def compare_jobs(left: str, right: str) -> dict:
    """Compare two completed calculations without conflating their methods."""
    left_info = manager.get(left)
    right_info = manager.get(right)
    if left_info is None or right_info is None:
        raise HTTPException(status_code=404, detail="comparison project not found")
    terminal = {JobStatus.done, JobStatus.not_converged}
    if left_info.status not in terminal or right_info.status not in terminal:
        raise HTTPException(status_code=409, detail="comparison requires completed projects")
    _summary_cache.pop(left, None)
    _summary_cache.pop(right, None)
    left_paths = _job_paths(left)
    right_paths = _job_paths(right)
    left_summary = _job_summary(left)
    right_summary = _job_summary(right)
    left_energy = _summary_energy(left_summary)
    right_energy = _summary_energy(right_summary)
    state_left = {row["index"]: row for row in left_summary.get("states", [])}
    state_right = {row["index"]: row for row in right_summary.get("states", [])}
    states = []
    for index in sorted(state_left.keys() & state_right.keys()):
        left_ev = state_left[index].get("excitation_ev")
        right_ev = state_right[index].get("excitation_ev")
        if not isinstance(left_ev, (int, float)) or not isinstance(right_ev, (int, float)):
            continue
        states.append({
            "index": index, "left_ev": left_ev, "right_ev": right_ev,
            "delta_ev": right_ev - left_ev,
        })
    left_dipole = left_summary.get("dipole") or {}
    right_dipole = right_summary.get("dipole") or {}
    return {
        "left": {"id": left, "name": left_info.name, "energy": left_energy},
        "right": {"id": right, "name": right_info.name, "energy": right_energy},
        "energy_delta_hartree": (
            right_energy - left_energy
            if left_energy is not None and right_energy is not None else None
        ),
        "energy_delta_kcal_mol": (
            (right_energy - left_energy) * 627.509474
            if left_energy is not None and right_energy is not None else None
        ),
        "geometry": _geometry_comparison(left_paths, right_paths),
        "states": states,
        "dipole_delta_debye": (
            right_dipole["total_debye"] - left_dipole["total_debye"]
            if "total_debye" in left_dipole and "total_debye" in right_dipole else None
        ),
    }


@app.get("/api/jobs/{job_id}/summary")
def job_summary(job_id: str, refresh: bool = False) -> dict:
    """Energies, states, frequencies, thermochemistry and properties."""
    if refresh:
        _summary_cache.pop(job_id, None)
        _scan_energy_cache.pop(job_id, None)
    return _job_summary(job_id)


@app.get("/api/jobs/{job_id}/optimization")
def job_optimization(job_id: str) -> dict:
    """Geometry, convergence measures and response data for each opt step."""
    from . import analysis

    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    return analysis.optimization_history(_job_paths(job_id))


@app.get("/api/jobs/{job_id}/excited-analysis")
def job_excited_analysis(job_id: str, ref: int = 0, target: int = 1) -> dict:
    """Physical-root MRSF NTO and attachment/detachment descriptors."""
    from . import excited_state

    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    try:
        return excited_state.summary(_job_paths(job_id), ref, target)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/jobs/{job_id}/excited-analysis/cube")
def job_excited_analysis_cube(job_id: str, kind: str = "nto_hole",
                              ref: int = 0, target: int = 1,
                              rank: int = 0) -> PlainTextResponse:
    """One MRSF physical-root orbital or density as a Gaussian cube."""
    from . import excited_state

    if manager.get(job_id) is None:
        raise HTTPException(status_code=404, detail="job not found")
    try:
        data = excited_state.ExcitedStateData.load(_job_paths(job_id))
        return PlainTextResponse(data.cube(kind, ref, target, rank))
    except (excited_state.AnalysisUnavailable, ValueError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/jobs/{job_id}/spectrum")
def job_spectrum(job_id: str, kind: str = "ir", shape: str = "lorentzian",
                 fwhm: float | None = None, state: int = 1,
                 step: int | None = None) -> dict:
    """One broadened spectrum. Lorentzian is the default line shape."""
    from . import analysis, spectra

    if shape not in spectra.SHAPES:
        raise HTTPException(status_code=400, detail=f"unknown line shape: {shape}")
    summary = _job_summary(job_id)
    if step is not None and kind in {"absorption", "emission", "esa"}:
        selected = next((item for item in analysis.optimization_history(_job_paths(job_id))["steps"]
                         if item["index"] == step), None)
        if selected is None or len(selected["states"]) < 2:
            return {"available": False,
                    "reason": f"no state-resolved spectrum recorded for optimization step {step}"}
        summary = dict(summary)
        summary["states"] = selected["states"]
        summary["transitions"] = selected["transitions"]
        summary["has_states"] = True
        summary["has_oscillators"] = any(row.get("oscillator") is not None for row in selected["states"])
    data = analysis.spectrum(summary, kind, shape=shape, fwhm=fwhm, state=max(1, state))
    if step is not None and kind in {"absorption", "emission", "esa"} and data.get("available"):
        data["title"] = f"{data.get('title', 'Spectrum')} at optimization step {step}"
    return data


MAP_KINDS = ("density", "spin", "alpha", "beta", "esp")


@app.get("/api/jobs/{job_id}/molden/{name}/map")
def molden_map(job_id: str, name: str, kind: str = "density") -> PlainTextResponse:
    """A scalar-field cube: electron density, spin density, or the MEP."""
    from . import molden

    if kind not in MAP_KINDS:
        raise HTTPException(status_code=400, detail=f"unknown map kind: {kind}")
    data = _load_molden(job_id, name)
    if not data.supported:
        raise HTTPException(status_code=422, detail="spherical shells not yet supported")
    if kind == "esp":
        # Prefer charges the run itself exported over ones derived here.
        summary = _job_summary(job_id)
        charges = (summary["charges"].get("resp")
                   or summary["charges"].get("lowdin")
                   or summary["charges"].get("mulliken"))
        return PlainTextResponse(molden.esp_cube(data, charges))
    return PlainTextResponse(
        molden.density_cube(data, "total" if kind == "density" else kind))


@app.get("/api/jobs/{job_id}/molden/{name}/charges")
def molden_charges(job_id: str, name: str) -> dict:
    """Atomic partial charges, from the run when present, else Mulliken."""
    from . import molden

    data = _load_molden(job_id, name)
    if not data.supported:
        raise HTTPException(status_code=422, detail="spherical shells not yet supported")
    summary = _job_summary(job_id)
    for source in ("resp", "lowdin", "mulliken"):
        charges = summary["charges"].get(source)
        if charges and len(charges) == len(data.atoms):
            return {"source": source, "charges": charges}
    return {"source": "mulliken (computed here)",
            "charges": molden.mulliken_charges(data)}


# Hosts the app is allowed to hand to the system browser. The Tauri webview
# blocks window.open, so links are opened by the OS instead — and that is only
# safe for destinations the app itself chose.
EXTERNAL_HOSTS = frozenset({
    "github.com",
    "open-quantum-platform.github.io",
    "www.github.com",
    "docs.openqp.org",
    "app.openqp.org",
    "openqp.org",
    "www.rcsb.org",
    "pubchem.ncbi.nlm.nih.gov",
})


# Progress of an engine download, polled by the UI.
_engine_state: dict = {"status": "idle", "detail": "", "percent": 0}


def _install_engine(assets: list[dict]) -> None:
    try:
        asset = engine.pick_asset(assets)
        if asset is None:
            _engine_state.update(
                status="failed",
                detail=f"no engine archive published for {engine.archive_suffix()}")
            return

        def progress(done: int, total: int) -> None:
            _engine_state.update(
                status="downloading",
                percent=int(done * 100 / total) if total else 0,
                detail=f"{done // 1_000_000} MB of {total // 1_000_000} MB")

        _engine_state.update(status="downloading", percent=0, detail=asset["name"])
        path = engine.install(asset["url"], progress, asset=asset)
        _engine_state.update(status="ready", percent=100, detail=path)
    except Exception as exc:  # noqa: BLE001 — the UI shows whatever went wrong
        _engine_state.update(status="failed", detail=str(exc))


@app.get("/api/engine")
def engine_status() -> dict:
    """Whether a compute engine is available, and where it came from."""
    return engine.status()


@app.post("/api/engine/install")
def engine_install() -> dict:
    """Download the standalone engine archive that matches this machine."""
    if _engine_state["status"] == "downloading":
        return _engine_state
    if not engine.archive_suffix():
        raise HTTPException(status_code=400,
                            detail="no engine archive is published for this platform")
    try:
        assets = [a for a in releases.installed_manifest()["assets"] if a["kind"] == "engine"]
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"engine manifest unavailable: {exc}") from exc
    if not engine.pick_asset(assets):
        raise HTTPException(
            status_code=404,
            detail="this Studio release has no qualified engine archive for this platform")
    _engine_state.update(status="downloading", percent=0, detail="starting…")
    threading.Thread(target=_install_engine, args=(assets,), daemon=True).start()
    return _engine_state


@app.get("/api/engine/status")
def engine_progress() -> dict:
    return _engine_state


class ExternalLink(BaseModel):
    url: str


@app.post("/api/open-external")
def open_external(link: ExternalLink) -> dict:
    """Open a link in the user's normal browser."""
    import webbrowser

    parsed = urllib.parse.urlparse(link.url)
    if parsed.scheme != "https" or parsed.hostname not in EXTERNAL_HOSTS:
        raise HTTPException(status_code=400, detail=f"refusing to open {link.url}")
    # The platform opener is more dependable inside a frozen app than the
    # webbrowser module, which looks for a browser on PATH.
    import subprocess

    if sys.platform == "darwin":
        command = ["open", link.url]
    elif os.name == "nt":
        command = ["cmd", "/c", "start", "", link.url]
    else:
        command = ["xdg-open", link.url]
    try:
        subprocess.Popen(command, start_new_session=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return {"opened": link.url}
    except OSError:
        pass
    if not webbrowser.open(link.url):
        raise HTTPException(status_code=500, detail="no browser could be launched")
    return {"opened": link.url}


@app.post("/api/reveal-results")
def reveal_results() -> dict:
    """Show the results folder in the user's file manager.

    The path is this app's own results directory, never one the page names:
    an endpoint that opened any path a request asked for would let a web page
    open anything on the machine.
    """
    import subprocess

    from . import jobs as jobs_module

    target = str(jobs_module.JOBS_ROOT)
    if sys.platform == "darwin":
        command = ["open", target]
    elif os.name == "nt":
        command = ["explorer", target]
    else:
        command = ["xdg-open", target]
    try:
        subprocess.Popen(command, start_new_session=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as exc:
        raise HTTPException(status_code=500,
                            detail=f"could not open {target}: {exc}") from exc
    return {"opened": target}


class NetworkSettings(BaseModel):
    ca_bundle: str = ""
    insecure: bool = False


@app.get("/api/network")
def get_network() -> dict:
    """What the app is doing about TLS, for the network settings dialog."""
    return network.status()


@app.post("/api/network")
def set_network(settings: NetworkSettings) -> dict:
    if settings.ca_bundle and not Path(settings.ca_bundle).is_file():
        raise HTTPException(status_code=400,
                            detail=f"no such file: {settings.ca_bundle}")
    network.save(settings.ca_bundle, settings.insecure)
    return network.status()


@app.post("/api/network/certificate")
async def upload_certificate(file: UploadFile = _UPLOAD) -> dict:
    """Store a root certificate the user picked in the file dialog.

    A browser file picker hands over the contents, not a path, so the file is
    copied next to the settings and the stored path points there.
    """
    data = (await file.read()).decode("utf-8", errors="replace")
    if "BEGIN CERTIFICATE" not in data:
        raise HTTPException(
            status_code=400,
            detail="that file is not a PEM certificate (no BEGIN CERTIFICATE block)",
        )
    target = network.settings_path().with_name("network-root.pem")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(data)
    network.save(str(target), network.load()["insecure"])
    return network.status()


@app.get("/api/network/test")
def test_network() -> dict:
    """Try the PubChem host and report exactly what happened."""
    try:
        with urllib.request.urlopen(PUBCHEM_URL.format(name="water"), timeout=15,
                                    context=network.context()):
            return {"ok": True, "detail": "PubChem reachable and verified"}
    except urllib.error.HTTPError as exc:
        # The host answered, which is all this test needs to establish.
        return {"ok": True, "detail": f"PubChem reachable (HTTP {exc.code})"}
    except Exception as exc:  # noqa: BLE001 — the message is the whole point
        return {"ok": False, "detail": str(exc)}


RELEASES_API = releases.API + "/latest"
RELEASES_PAGE = releases.PAGE


def _version_tuple(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"\d+", text)[:3]) or (0,)


@app.get("/api/update-check")
def update_check() -> dict:
    """Offer only a published Studio release with a validated binary manifest."""
    try:
        release, manifest = releases.latest()
        latest = manifest["version"]
        available = releases.version_key(latest) > releases.version_key(__version__)
        return {
            "current": __version__, "latest": latest, "available": available,
            "url": release.get("html_url") or RELEASES_PAGE,
            "notes": (release.get("body") or "")[:400],
            "channel": releases.channel(),
            "assets": [a for a in manifest["assets"] if a["kind"] == "installer"],
        }
    except Exception as exc:  # noqa: BLE001 — offline or invalid feed is unavailable
        return {"current": __version__, "available": False,
                "detail": f"could not check the public Studio feed: {exc}",
                "url": RELEASES_PAGE, "channel": releases.channel()}


class UpdateChannel(BaseModel):
    channel: str


@app.get("/api/update/channel")
def update_channel() -> dict:
    return {"channel": releases.channel()}


@app.post("/api/update/channel")
def change_update_channel(req: UpdateChannel) -> dict:
    try:
        releases.set_channel(req.channel)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return update_channel()


# Progress of a running update, polled by the UI.
_update_state: dict = {"status": "idle", "detail": "", "percent": 0}


def _run_update(assets: list[dict]) -> None:
    from . import update

    try:
        asset = update.pick_asset(assets)
        if asset is None:
            _update_state.update(status="failed",
                                 detail="no installer published for this platform")
            return
        target = Path(tempfile.mkdtemp(prefix="oqp-studio-update-download-")) / asset["name"]

        def progress(done: int, total: int) -> None:
            _update_state.update(
                status="downloading",
                percent=int(done * 100 / total) if total else 0,
                detail=f"{done // 1_000_000} MB of {total // 1_000_000} MB",
            )

        _update_state.update(status="downloading", percent=0, detail=asset["name"])
        releases.download(asset, target, progress)

        if sys.platform == "darwin":
            _update_state.update(status="installing", percent=100, detail="unpacking…")
            detail = update.install_macos(
                target, lambda text: _update_state.update(detail=text))
            _update_state.update(status="ready", detail=detail)
        else:
            detail = update.install_elsewhere(target)
            _update_state.update(status="ready", detail=detail)
    except Exception as exc:  # noqa: BLE001 — the UI shows whatever went wrong
        _update_state.update(status="failed", detail=str(exc))


@app.post("/api/update/install")
def update_install() -> dict:
    """Download this platform's installer and stage it for installation."""
    if _update_state["status"] in ("downloading", "installing"):
        return _update_state
    info = update_check()
    if not info.get("available"):
        raise HTTPException(status_code=400, detail="already up to date")
    _update_state.update(status="downloading", percent=0, detail="starting…")
    threading.Thread(target=_run_update, args=(info.get("assets") or [],),
                     daemon=True).start()
    return _update_state


@app.get("/api/update/status")
def update_status() -> dict:
    return _update_state


@app.post("/api/structure/open")
async def open_structure(file: UploadFile = _UPLOAD) -> dict:
    """Read a structure out of whatever file the user already has.

    Handles OpenQP's own inputs and outputs (.oqp, .inp, .log, .json,
    .molden, packed .namd.trj), the common interchange formats (.xyz, .pdb,
    .mol/.sdf, .mol2), and ChemDraw XML. Formats that carry several
    geometries come back as multiple frames.
    """
    from . import structure_io

    data = await file.read()
    name = file.filename or "structure"
    temp_path = None
    try:
        if name.lower().endswith((".trj", ".namd.trj")):
            # The packed record is read through OpenQP's memory-mapped reader,
            # which needs a real path.
            with tempfile.NamedTemporaryFile(suffix=".namd.trj", delete=False) as handle:
                handle.write(data)
                temp_path = handle.name
        structure = structure_io.parse(name, data, path=temp_path)
    except structure_io.UnsupportedFormat as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:  # any parse failure is about the user's file
        raise HTTPException(status_code=422, detail=f"could not read {name}: {exc}") from exc
    finally:
        if temp_path:
            Path(temp_path).unlink(missing_ok=True)

    return {
        "name": name,
        "format": structure.format,
        "frames": [
            {"label": frame.label, "atoms": [[a[0], a[1], a[2], a[3]] for a in frame.atoms]}
            for frame in structure.frames
        ],
    }


class Structure3DRequest(BaseModel):
    molfile: str | None = None
    smiles: str | None = None


@app.post("/api/structure3d")
def structure_3d(req: Structure3DRequest) -> dict:
    """2D sketch (molfile) or SMILES -> 3D coordinates via RDKit ETKDG+MMFF."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise HTTPException(
            status_code=501,
            detail="RDKit is not installed on the backend (pip install rdkit)",
        ) from exc

    if req.molfile:
        mol = Chem.MolFromMolBlock(req.molfile)
    elif req.smiles:
        mol = Chem.MolFromSmiles(req.smiles)
    else:
        raise HTTPException(status_code=422, detail="molfile or smiles required")
    if mol is None:
        raise HTTPException(status_code=422, detail="could not parse the structure")

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        raise HTTPException(status_code=422, detail="3D embedding failed for this structure")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:  # noqa: BLE001, S110 — unoptimized coordinates are still usable
        pass
    conf = mol.GetConformer()
    atoms = [
        [atom.GetSymbol(), *(round(v, 6) for v in conf.GetAtomPosition(i))]
        for i, atom in enumerate(mol.GetAtoms())
    ]
    return {"atoms": atoms}


# In-memory store of geometry previews served back to the 3D viewer.
_previews: OrderedDict[str, str] = OrderedDict()


class PreviewRequest(BaseModel):
    xyz: str


@app.post("/api/preview")
def create_preview(req: PreviewRequest) -> dict:
    preview_id = uuid.uuid4().hex[:12]
    _previews[preview_id] = req.xyz
    while len(_previews) > 50:
        _previews.popitem(last=False)
    return {"url": f"/api/preview/{preview_id}.xyz"}


@app.get("/api/preview/{name}")
def get_preview(name: str) -> PlainTextResponse:
    preview = _previews.get(name.removesuffix(".xyz"))
    if preview is None:
        raise HTTPException(status_code=404, detail="preview not found")
    return PlainTextResponse(preview)


def _frontend_dist() -> Path | None:
    """Locate the built frontend so one server (and one origin) serves both
    the UI and the API — the layout the desktop shell relies on. In dev, the
    Vite server proxies /api instead and this mount is absent."""
    env = os.environ.get("OQP_STUDIO_FRONTEND")
    candidates = [Path(env)] if env else []
    bundle_dir = getattr(sys, "_MEIPASS", None)  # PyInstaller-frozen backend
    if bundle_dir:
        candidates.append(Path(bundle_dir) / "frontend_dist")
    candidates.append(Path(__file__).resolve().parent / "web")  # installed wheel
    candidates.append(Path(__file__).resolve().parents[2] / "frontend" / "dist")
    return next((p for p in candidates if (p / "index.html").is_file()), None)


class _NoStoreStaticFiles(StaticFiles):
    """Serves the UI with caching disabled.

    The desktop shell always loads the same URL (http://127.0.0.1:<port>/)
    under the same bundle identifier, so WebKit would otherwise reuse the
    index.html it cached from a previously installed version and show a stale
    interface after an upgrade.
    """

    async def get_response(self, path: str, scope):
        response = await super().get_response(path, scope)
        response.headers["Cache-Control"] = "no-store, must-revalidate"
        return response


_dist = _frontend_dist()
if _dist is not None:
    app.mount("/", _NoStoreStaticFiles(directory=_dist, html=True), name="frontend")
