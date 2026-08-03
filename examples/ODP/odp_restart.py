#!/usr/bin/env python3
"""Run one ODP step, then continue through step two from its manifest."""

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
GEOMETRY = HERE.parent / "geometries" / "CH2O-2bc62dda4b8a.xyz"
TEMPLATE = HERE / "odp_restart.oqp.template"


def run_example(workdir):
    """Run a fresh checkpoint and its generated, one-step-longer manifest."""
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    source = workdir / "odp-restart.oqp"
    geometry = workdir / GEOMETRY.name
    source.write_text(TEMPLATE.read_text(encoding="utf-8"), encoding="utf-8")
    shutil.copy2(GEOMETRY, geometry)

    command = [
        sys.executable, "-m", "oqp.pyoqp", "--nompi", "--omp", "1", "--silent"]
    environment = os.environ.copy()
    source_package = HERE.parents[1] / "pyoqp"
    if source_package.is_dir():
        current = environment.get("PYTHONPATH")
        environment["PYTHONPATH"] = str(source_package) + (
            os.pathsep + current if current else "")
    subprocess.run(
        [*command, str(source)], cwd=workdir, env=environment, check=True)

    manifest = workdir / "odp-restart.restart.oqp"
    if not manifest.is_file():
        raise RuntimeError("fresh ODP run did not generate a restart manifest")
    text = manifest.read_text(encoding="utf-8")
    resumed = re.sub(r"\bnstep=1\b", "nstep=2", text, count=1)
    if resumed == text or "restart=true" not in resumed:
        raise RuntimeError("generated manifest lacks the expected restart contract")
    manifest.write_text(resumed, encoding="utf-8")
    subprocess.run(
        [*command, str(manifest)], cwd=workdir, env=environment, check=True)

    checkpoint = workdir / "odp-restart.namd.restart.npz"
    with np.load(checkpoint, allow_pickle=False) as saved:
        step = int(saved["step"][0])
    if step != 2:
        raise RuntimeError(f"ODP restart stopped at step {step}, expected step 2")
    print(f"ODP restart completed through step {step}: {checkpoint}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workdir", help="keep inputs, manifest, checkpoint, and trajectory here")
    args = parser.parse_args()
    if args.workdir:
        run_example(args.workdir)
    else:
        with tempfile.TemporaryDirectory(prefix="openqp-odp-restart-") as workdir:
            run_example(workdir)


if __name__ == "__main__":
    main()
