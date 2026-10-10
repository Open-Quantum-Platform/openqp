from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import threading
import uuid
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from shutil import rmtree

from pydantic import BaseModel, Field

from .input_files import SUPPORTED_INPUT_SUFFIXES, calculation_log, find_input_file
from .runners import get_runner
from .structure_io import parse_xyz


def _preferred_root() -> Path:
    from . import workspace

    return workspace.preferred()


# Worked out from the environment and the saved setting alone -- importing
# this module must not touch the disk. The directory is created, and any
# fallback chosen, the first time a job actually needs it.
JOBS_ROOT = _preferred_root()
_METADATA_FILE = ".oqp-studio.json"


class JobStatus(str, Enum):
    queued = "queued"
    running = "running"
    cancelling = "cancelling"
    cancelled = "cancelled"
    done = "done"
    not_converged = "not_converged"
    failed = "failed"


class JobRequest(BaseModel):
    input_text: str
    runner: str = "local"
    name: str = "job"
    input_name: str | None = None
    pdb_text: str | None = None
    pdb_name: str | None = None
    threads: int = Field(default=1, ge=1)
    recipe_id: str | None = None


class JobInfo(BaseModel):
    id: str
    name: str
    status: JobStatus
    runner: str
    threads: int = 1
    created_at: str
    exit_code: int | None = None
    error: str | None = None
    group_id: str | None = None
    scan_value: float | None = None
    scan_unit: str | None = None
    scan_state: int | None = None
    recipe_id: str | None = None
    postprocess_status: str | None = None
    postprocess_error: str | None = None


class JobManager:
    """Minimal in-memory job queue: one worker thread per job (Phase 0).

    Phase 1 replaces this with a persistent queue with concurrency limits,
    kill/requeue, and WebSocket progress events.
    """

    def __init__(self) -> None:
        self._jobs: dict[str, JobInfo] = {}
        self._lock = threading.Lock()
        self._processes: dict[str, subprocess.Popen] = {}
        self._cancel_requested: set[str] = set()
        self._ready = False

    def _ensure(self) -> Path:
        """Create the results directory, once, on the first job that needs it.

        Deliberately not done in __init__: this module is imported while the
        server is starting, and on macOS creating a directory in Documents
        waits for the user to grant permission -- which the shell, waiting
        thirty seconds for a port, reads as the backend failing to start.
        """
        global JOBS_ROOT

        if self._ready:
            return JOBS_ROOT
        from . import workspace

        JOBS_ROOT = workspace.ensure(JOBS_ROOT)
        self._ready = True
        self._recover()
        return JOBS_ROOT

    def _recover(self) -> None:
        """Adopt job directories left by earlier runs.

        Outputs outlive the process, so a restart should not empty the job
        list and leave the analysis tab with nothing to open.
        """
        if not JOBS_ROOT.is_dir():
            return
        for job_dir in sorted(JOBS_ROOT.iterdir()):
            if not job_dir.is_dir() or not any(job_dir.iterdir()):
                continue
            created = datetime.fromtimestamp(
                job_dir.stat().st_mtime, timezone.utc).isoformat()
            metadata = self._metadata(job_dir)
            diagnostic = self._optimization_diagnostic(job_dir)
            saved_status = metadata.get("status")
            interrupted = saved_status in {"queued", "running", "cancelling"}
            if interrupted:
                recovered_status = JobStatus.cancelled
            elif saved_status in {status.value for status in JobStatus}:
                recovered_status = JobStatus(saved_status)
            else:
                recovered_status = JobStatus.not_converged if diagnostic else JobStatus.done
            self._jobs[job_dir.name] = JobInfo(
                id=job_dir.name,
                name=str(metadata.get("name") or self._recovered_name(job_dir)),
                status=recovered_status,
                runner=str(metadata.get("runner") or "recovered"),
                threads=int(metadata.get("threads") or 1),
                created_at=str(metadata.get("created_at") or created),
                exit_code=metadata.get("exit_code"),
                error=("Interrupted before the backend restarted" if interrupted
                       else diagnostic or metadata.get("error")),
                group_id=metadata.get("group_id"),
                scan_value=metadata.get("scan_value"),
                scan_unit=metadata.get("scan_unit"),
                scan_state=metadata.get("scan_state"),
                recipe_id=metadata.get("recipe_id"),
                postprocess_status=metadata.get("postprocess_status"),
                postprocess_error=metadata.get("postprocess_error"),
            )

    @staticmethod
    def _metadata(job_dir: Path) -> dict:
        try:
            data = json.loads((job_dir / _METADATA_FILE).read_text())
        except (OSError, ValueError):
            return {}
        return data if isinstance(data, dict) else {}

    @staticmethod
    def _recovered_name(job_dir: Path) -> str:
        input_file = find_input_file(job_dir)
        if input_file.is_file() and input_file.stem.lower() != "input":
            return input_file.stem
        logs = sorted(path for path in job_dir.glob("*.log") if path.name != "job.log")
        return logs[0].stem if logs else job_dir.name

    @staticmethod
    def _write_metadata(job_dir: Path, info: JobInfo) -> None:
        try:
            (job_dir / _METADATA_FILE).write_text(json.dumps({
                "name": info.name,
                "runner": info.runner,
                "threads": info.threads,
                "created_at": info.created_at,
                "group_id": info.group_id,
                "scan_value": info.scan_value,
                "scan_unit": info.scan_unit,
                "scan_state": info.scan_state,
                "recipe_id": info.recipe_id,
                "postprocess_status": info.postprocess_status,
                "postprocess_error": info.postprocess_error,
                "status": info.status.value,
                "exit_code": info.exit_code,
                "error": info.error,
            }))
        except FileNotFoundError:
            # A test teardown or an external cleanup can remove the directory
            # while a runner is returning. There is nowhere left to persist.
            return

    @staticmethod
    def _validate_request(req: JobRequest) -> None:
        from . import engine, host

        engine.check_capabilities(req.input_text, req.runner)
        check = host.admission(req.input_text, req.threads)
        if not check["permitted"]:
            raise ValueError(
                f"RAM limit: estimated {check['estimated_memory_bytes'] / 1024**3:.1f} GiB, "
                f"available {check['memory_available_bytes'] / 1024**3:.1f} GiB. "
                "Choose a smaller calculation or free memory before running."
            )
        if req.recipe_id:
            from . import recipes

            record = recipes.get(req.recipe_id)
            if record is None:
                raise ValueError(f"unknown recipe '{req.recipe_id}'")
            if not record.recipe.available:
                raise ValueError(record.recipe.unavailable_reason or "recipe is unavailable")

    def _prepare(self, req: JobRequest, *, group_id: str | None = None,
                 scan_value: float | None = None, scan_unit: str | None = None,
                 scan_state: int | None = None) -> JobInfo:
        self._ensure()
        job_id = uuid.uuid4().hex[:12]
        job_dir = JOBS_ROOT / job_id
        try:
            job_dir.mkdir(parents=True, exist_ok=True)
            input_name = self._input_name(req)
            (job_dir / input_name).write_text(req.input_text)
            if req.pdb_text is not None:
                pdb_name = self._pdb_name(req)
                (job_dir / pdb_name).write_text(req.pdb_text)
            if req.recipe_id:
                from . import recipes

                recipes.snapshot(req.recipe_id, job_dir / ".oqp-studio-recipe.json")
        except Exception:
            rmtree(job_dir, ignore_errors=True)
            raise
        info = JobInfo(
            id=job_id,
            name=req.name,
            status=JobStatus.queued,
            runner=req.runner,
            threads=req.threads,
            created_at=datetime.now(timezone.utc).isoformat(),
            group_id=group_id,
            scan_value=scan_value,
            scan_unit=scan_unit,
            scan_state=scan_state,
            recipe_id=req.recipe_id,
        )
        with self._lock:
            self._jobs[job_id] = info
        try:
            self._write_metadata(job_dir, info)
        except Exception:
            with self._lock:
                self._jobs.pop(job_id, None)
            rmtree(job_dir, ignore_errors=True)
            raise
        return info

    def submit(self, req: JobRequest) -> JobInfo:
        self._validate_request(req)
        info = self._prepare(req)
        threading.Thread(target=self._run, args=(info.id,), daemon=True).start()
        return info

    def submit_batch(self, requests: list[JobRequest], *, group_id: str,
                     values: list[float], unit: str, state: int | None = None) -> list[JobInfo]:
        """Prepare a scan as normal jobs and run its points serially."""
        if not requests or len(requests) != len(values):
            raise ValueError("scan requests and coordinate values must have equal non-zero length")
        for request in requests:
            self._validate_request(request)
        infos: list[JobInfo] = []
        try:
            for request, value in zip(requests, values):
                infos.append(self._prepare(
                    request, group_id=group_id, scan_value=value, scan_unit=unit,
                    scan_state=state,
                ))
        except Exception:
            with self._lock:
                for info in infos:
                    self._jobs.pop(info.id, None)
            for info in infos:
                rmtree(JOBS_ROOT / info.id, ignore_errors=True)
            raise
        threading.Thread(
            target=self._run_batch, args=([info.id for info in infos],), daemon=True
        ).start()
        return infos

    def _run_batch(self, job_ids: list[str]) -> None:
        for job_id in job_ids:
            self._run(job_id)
            if self._jobs[job_id].status == JobStatus.cancelled:
                break

    @staticmethod
    def _input_name(req: JobRequest) -> str:
        """Choose a default input name or validate the name supplied by the UI."""
        if not req.input_name or not req.input_name.strip():
            # Sectioned legacy input starts with a [section]; otherwise .oqp route style.
            return "input.inp" if req.input_text.lstrip().startswith("[") else "input.oqp"
        supplied = req.input_name.strip()
        name = Path(supplied).name
        if name != supplied or name.startswith(".") or name in {"", ".", ".."}:
            raise ValueError("input file name must not contain a path")
        if Path(name).suffix.lower() not in SUPPORTED_INPUT_SUFFIXES:
            raise ValueError("input file name must end in .oqp or .inp")
        return name

    @staticmethod
    def _pdb_name(req: JobRequest) -> str:
        """Validate the optional PDB asset referenced by a QM/MM input."""
        supplied = (req.pdb_name or "structure.pdb").strip()
        name = Path(supplied).name
        if name != supplied or name.startswith(".") or Path(name).suffix.lower() != ".pdb":
            raise ValueError("PDB file name must be a .pdb file name without a path")
        return name

    def adopt(self, name: str, files: list[tuple[str, bytes]]) -> JobInfo:
        """Register result files computed elsewhere as a job of their own.

        Analysis is built entirely on top of a job directory, so a run made
        on a cluster, by the command-line engine or by an earlier version of
        this app becomes analysable simply by giving its output files a job
        directory of their own -- summaries, spectra, orbitals, normal modes
        and property maps then work on it unchanged.
        """
        self._ensure()
        job_id = uuid.uuid4().hex[:12]
        job_dir = JOBS_ROOT / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        written = 0
        for raw_name, data in files:
            # Browsers send a path when a whole folder is picked; only the
            # file name may reach the job directory.
            safe = Path(raw_name.replace("\\", "/")).name
            if not safe or safe.startswith("."):
                continue
            (job_dir / safe).write_bytes(data)
            written += 1
        if not written:
            job_dir.rmdir()
            raise ValueError("no usable files were given")
        info = JobInfo(
            id=job_id,
            name=name or "imported",
            status=JobStatus.done,
            runner="imported",
            created_at=datetime.now(timezone.utc).isoformat(),
            exit_code=0,
        )
        with self._lock:
            self._jobs[job_id] = info
        self._write_metadata(job_dir, info)
        return info

    def _run(self, job_id: str) -> None:
        info = self._jobs[job_id]
        with self._lock:
            if job_id in self._cancel_requested:
                info.status = JobStatus.cancelled
                info.error = "Cancelled by user"
                self._write_metadata(JOBS_ROOT / job_id, info)
                return
        try:
            info.status = JobStatus.running
            self._write_metadata(JOBS_ROOT / job_id, info)
            exit_code = get_runner(info.runner).run(
                JOBS_ROOT / job_id, info.threads,
                on_start=lambda process: self._register_process(job_id, process),
            )
            info.exit_code = exit_code
            if job_id in self._cancel_requested:
                info.status = JobStatus.cancelled
                info.error = "Cancelled by user"
                return
            diagnostic = self._optimization_diagnostic(JOBS_ROOT / job_id)
            if exit_code == 0 and diagnostic:
                info.status = JobStatus.not_converged
                info.error = diagnostic
            elif exit_code != 0:
                info.status = JobStatus.failed
            elif info.recipe_id:
                from . import recipes

                info.postprocess_status = "running"
                self._write_metadata(JOBS_ROOT / job_id, info)
                try:
                    info.postprocess_status, info.postprocess_error = recipes.run_snapshot(
                        JOBS_ROOT / job_id / ".oqp-studio-recipe.json",
                        JOBS_ROOT / job_id / "postprocess.log",
                    )
                except (OSError, ValueError) as exc:
                    info.postprocess_status = "failed"
                    info.postprocess_error = str(exc)
                info.status = JobStatus.done
            else:
                info.status = JobStatus.done
        except Exception as exc:  # noqa: BLE001 — any failure must land in the job record
            info.error = str(exc)
            info.status = JobStatus.failed
        finally:
            try:
                self._write_metadata(JOBS_ROOT / job_id, info)
            except OSError:
                # The in-memory terminal state must remain usable even when a
                # full or read-only results volume cannot persist it.
                pass
            with self._lock:
                self._processes.pop(job_id, None)

    def _register_process(self, job_id: str, process: subprocess.Popen) -> None:
        with self._lock:
            self._processes[job_id] = process
            cancel_requested = job_id in self._cancel_requested
        if cancel_requested:
            self._terminate_process(process)

    @staticmethod
    def _terminate_process(process: subprocess.Popen) -> None:
        if process.poll() is not None:
            return
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], check=False,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            os.killpg(process.pid, signal.SIGTERM)

    def cancel(self, job_id: str) -> None:
        with self._lock:
            info = self._jobs.get(job_id)
            if info is None:
                raise KeyError(job_id)
            if info.status not in {JobStatus.queued, JobStatus.running, JobStatus.cancelling}:
                raise ValueError("only a queued or running calculation can be cancelled")
            affected = [candidate for candidate in self._jobs.values()
                        if candidate.id == job_id or (
                            info.group_id is not None and candidate.group_id == info.group_id
                        )]
            processes = []
            for candidate in affected:
                if candidate.status not in {
                    JobStatus.queued, JobStatus.running, JobStatus.cancelling
                }:
                    continue
                self._cancel_requested.add(candidate.id)
                process = self._processes.get(candidate.id)
                if process is None and candidate.status == JobStatus.queued:
                    candidate.status = JobStatus.cancelled
                    candidate.error = "Cancelled with scan group" if info.group_id else "Cancelled by user"
                else:
                    candidate.status = JobStatus.cancelling
                    if process is not None:
                        processes.append(process)
                self._write_metadata(JOBS_ROOT / candidate.id, candidate)
        for process in processes:
            self._terminate_process(process)

    @staticmethod
    def _optimization_diagnostic(job_dir: Path) -> str | None:
        """Return the final native-optimizer nonconvergence record, if any."""
        input_file = find_input_file(job_dir)
        logs = [calculation_log(input_file)]
        logs.extend(sorted(path for path in job_dir.glob("*.log") if path not in logs))
        pattern = re.compile(
            r"(?:Native optimization.*did not converge|BaekA.*did not converge|"
            r"Geometry Optimization Has Not Converged|MECP SQP did not converge)",
            re.IGNORECASE,
        )
        matches: list[str] = []
        for log_file in logs:
            if not log_file.is_file():
                continue
            for line in log_file.read_text(errors="replace").splitlines():
                if pattern.search(line):
                    matches.append(" ".join(line.split()))
        return matches[-1] if matches else None

    def restart_input(self, job_id: str) -> dict:
        """Prepare canonical restart input from a nonconverged job's opt.xyz."""
        source = self._jobs.get(job_id)
        self._ensure()
        source = source or self._jobs.get(job_id)
        if source is None:
            raise KeyError(job_id)
        if source.status != JobStatus.not_converged:
            raise ValueError("only a nonconverged geometry job can be restarted")
        source_dir = JOBS_ROOT / job_id
        input_file = find_input_file(source_dir)
        if input_file.suffix.lower() != ".oqp":
            raise ValueError("restart requires OpenQP canonical .oqp input")
        geometry_file = source_dir / "opt.xyz"
        if not geometry_file.is_file():
            raise ValueError("the retained optimization geometry (opt.xyz) is missing")
        frames = parse_xyz(geometry_file.read_text(errors="replace"))
        if not frames or not frames[-1].atoms:
            raise ValueError("could not read the retained optimization geometry")
        geometry = "\n".join(
            f"{element:<2} {x:11.6f} {y:11.6f} {z:11.6f}"
            for element, x, y, z in frames[-1].atoms
        )
        input_text = input_file.read_text(errors="replace")
        restarted, count = re.subn(
            r'geom(?:etry)?\s*=\s*""".*?"""',
            f'geom="""\n{geometry}\n"""',
            input_text,
            count=1,
            flags=re.DOTALL | re.IGNORECASE,
        )
        if count != 1:
            raise ValueError("restart requires an inline geom=\"\"\"...\"\"\" block")
        return {
            "input_text": restarted,
            "input_name": input_file.name,
            "name": f"{source.name} restart",
            "runner": source.runner if source.runner != "recovered" else None,
            "threads": source.threads,
        }

    def rebase(self) -> None:
        """Forget the old results directory and read the new one."""
        with self._lock:
            self._jobs.clear()
        self._recover()

    def busy(self) -> bool:
        return any(job.status in (JobStatus.queued, JobStatus.running)
                   for job in self._jobs.values())

    def get(self, job_id: str) -> JobInfo | None:
        return self._jobs.get(job_id)

    def list(self) -> list[JobInfo]:
        self._ensure()
        return sorted(self._jobs.values(), key=lambda j: j.created_at, reverse=True)

    def delete(self, job_id: str) -> None:
        """Remove a completed project's record and its result directory."""
        self._ensure()
        with self._lock:
            info = self._jobs.get(job_id)
            if info is None:
                raise KeyError(job_id)
            if info.status in (JobStatus.queued, JobStatus.running, JobStatus.cancelling):
                raise ValueError("a running calculation cannot be deleted")
            root = JOBS_ROOT.resolve()
            job_dir = (root / job_id).resolve()
            if job_dir.parent != root:
                raise ValueError("invalid project directory")
            rmtree(job_dir)
            del self._jobs[job_id]

    def files(self, job_id: str) -> list[dict]:
        self._ensure()
        job_dir = JOBS_ROOT / job_id
        if not job_dir.is_dir():
            return []
        return sorted(
            (
                {"name": p.name, "size": p.stat().st_size}
                for p in job_dir.iterdir()
                if p.is_file() and not p.name.startswith(".")
            ),
            key=lambda f: f["name"],
        )

    def file_path(self, job_id: str, name: str) -> Path | None:
        """Resolve a job output file, refusing path escapes."""
        job_dir = (JOBS_ROOT / job_id).resolve()
        candidate = (job_dir / name).resolve()
        if candidate.parent != job_dir or not candidate.is_file():
            return None
        return candidate

    def log_tail(self, job_id: str, max_bytes: int = 65536) -> str:
        job_dir = JOBS_ROOT / job_id
        # job.log is the runner process's stdout/stderr capture. OpenQP writes
        # its actual calculation record separately, so show that record when
        # it has appeared rather than the launcher transcript.
        input_file = find_input_file(job_dir)
        candidates = [calculation_log(input_file)]
        candidates.extend(
            sorted(
                path for suffix in ("*.log", "*.out")
                for path in job_dir.glob(suffix)
                if path not in {job_dir / "job.log", calculation_log(input_file)}
            )
        )
        candidates.append(job_dir / "job.log")
        log_file = next((path for path in candidates if path.is_file()), None)
        if log_file is None:
            return ""
        data = log_file.read_bytes()
        return data[-max_bytes:].decode(errors="replace").replace("\r\n", "\n").replace("\r", "\n")


manager = JobManager()


def set_root(directory: Path) -> Path:
    """Write results somewhere else from now on.

    A job already running keeps writing where it was told -- its runner has
    the old path -- which is why the API refuses to move while one is in
    flight rather than leaving its output stranded.
    """
    global JOBS_ROOT

    JOBS_ROOT = directory
    manager._ready = True
    manager.rebase()
    return JOBS_ROOT
