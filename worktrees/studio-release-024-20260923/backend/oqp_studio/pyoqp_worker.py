"""Worker process: run one OpenQP job directory through the pyoqp API.

Usage: python -m oqp_studio.pyoqp_worker <job_dir>

Reads the job's .oqp or .inp input, runs it with oqp.pyoqp.Runner, and writes
the calculation record next to that input (for example, input.log). The parent
also captures launcher stdout/stderr in job.log. A JSON summary of
runner.results() is written to <job_dir>/results.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from .input_files import calculation_log, find_input_file


def _jsonable(value):
    """Best-effort conversion of numpy arrays and other objects to JSON."""
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        return tolist()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def main() -> int:
    job_dir = Path(sys.argv[1]).resolve()
    input_file = find_input_file(job_dir)

    from oqp.pyoqp import Runner  # deferred: only the worker needs OpenQP

    runner = Runner(
        project=input_file.stem,
        input_file=str(input_file),
        log=str(calculation_log(input_file)),
    )
    runner.run()

    summary = _jsonable(runner.results())
    (job_dir / "results.json").write_text(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
