from __future__ import annotations

import runpy
import sys
from pathlib import Path


def run(code_path: Path, job_dir: Path) -> None:
    sys.path.insert(0, str(job_dir))
    runpy.run_path(str(code_path), init_globals={"JOB_DIR": job_dir}, run_name="__main__")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: postprocess_worker CODE_FILE JOB_DIR")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())


if __name__ == "__main__":
    main()
