#!/usr/bin/env python3
"""Create a matched Baeck-An cohort from validated baseline `.oqp` inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-input-root", type=Path, required=True)
    parser.add_argument("--baseline-jobs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()

    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    rows = []
    calculations = []
    lines = args.baseline_jobs.read_text(encoding="ascii").splitlines()
    if len(lines) != 50:
        parser.error(f"expected 50 baseline identities, found {len(lines)}")

    args.output.mkdir(parents=True)
    for expected_index, line in enumerate(lines):
        fields = line.split("\t")
        if len(fields) != 4 or int(fields[0]) != expected_index:
            parser.error(f"invalid baseline job row {expected_index}: {line}")
        source_relative = Path(fields[3])
        identity = source_relative.name
        source_dir = args.baseline_input_root / source_relative
        destination = args.output / "production_baeck_an" / "baeck_an" / identity
        shutil.copytree(source_dir, destination)

        input_path = destination / "uracil.oqp"
        original = input_path.read_text(encoding="ascii")
        marker = "tdc=npi,rescale=isotropic,nacme_check=off"
        replacement = "tdc=baeck_an,rescale=isotropic,nacme_check=off"
        if original.count(marker) != 1:
            parser.error(f"{source_dir}: expected exactly one baseline TDC marker")
        changed = original.replace(marker, replacement)
        if changed.replace(replacement, marker) != original:
            parser.error(f"{source_dir}: transformation changed more than the TDC provider")
        input_path.write_text(changed, encoding="ascii")

        relative = destination.relative_to(args.output)
        rows.append(
            f"{expected_index}\tproduction_baeck_an:baeck_an:{identity}"
            f"\tbaeck_an\t{relative}\n"
        )
        calculations.append({
            "array_index": expected_index,
            "identity": identity,
            "source_relative": str(source_relative),
            "relative_root": str(relative),
            "input_sha256": sha256(input_path),
            "geometry_sha256": sha256(destination / "geometry.xyz"),
            "velocity_sha256": sha256(destination / "velocity.au"),
            "initial_condition_sha256": sha256(
                destination / "initial-condition.json"),
            "environment_sha256": sha256(destination / "run.env"),
        })

    jobs = args.output / "production_baeck_an-jobs.tsv"
    jobs.write_text("".join(rows), encoding="ascii")
    manifest = {
        "schema": "openqp-uracil-baeck-an-v1",
        "input_format": "canonical .oqp",
        "source_commit": args.source_commit,
        "method": "baeck_an",
        "definition": {
            "tdc": "lagged Baeck-An magnitude with overlap-transported sign",
            "warmup": "overlap NPI on the first/reseeded interval",
            "rescale": "isotropic",
            "nstep": 1000,
            "dt_fs": 0.5,
        },
        "baseline_jobs_sha256": sha256(args.baseline_jobs),
        "jobs_tsv": jobs.name,
        "jobs_tsv_sha256": sha256(jobs),
        "calculations": calculations,
    }
    manifest_path = args.output / "input-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii")
    print(manifest_path)
    print(sha256(manifest_path))


if __name__ == "__main__":
    main()
