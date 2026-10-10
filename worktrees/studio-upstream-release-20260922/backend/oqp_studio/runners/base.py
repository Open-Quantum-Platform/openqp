from __future__ import annotations

import abc
import os
import subprocess
from collections.abc import Callable
from pathlib import Path

from ..input_files import find_input_file


class RunnerUnavailable(RuntimeError):
    """Raised when a runner cannot execute on this machine."""


class Runner(abc.ABC):
    """Runs OpenQP on a prepared job directory and streams output to job.log."""

    name: str = "abstract"

    @abc.abstractmethod
    def is_available(self) -> bool:
        """Whether this runner can execute OpenQP on this machine."""

    @abc.abstractmethod
    def build_command(self, input_file: Path) -> list[str]:
        """Command line that runs OpenQP on `input_file`."""

    def run(self, job_dir: Path, threads: int = 1,
            on_start: Callable[[subprocess.Popen[bytes]], None] | None = None) -> int:
        """Blocking execution; the job manager calls this from a worker thread."""
        if not self.is_available():
            raise RunnerUnavailable(f"runner '{self.name}' is not available on this machine")
        # Absolute: the child runs with cwd=job_dir, so a path relative to the
        # server's own directory would not resolve for it — and an engine that
        # writes its results next to the input would put them somewhere else.
        job_dir = job_dir.resolve()
        input_file = find_input_file(job_dir)
        log_file = job_dir / "job.log"
        with log_file.open("ab") as log:
            process_options: dict = {}
            if os.name == "nt":
                process_options["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
            else:
                process_options["start_new_session"] = True
            proc = subprocess.Popen(
                self.build_command(input_file),
                cwd=job_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env={**os.environ, "OMP_NUM_THREADS": str(threads)},
                **process_options,
            )
            if on_start is not None:
                on_start(proc)
            return proc.wait()
