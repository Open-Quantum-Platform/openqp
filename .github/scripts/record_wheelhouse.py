#!/usr/bin/env python3
"""Record the exact wheels used for one ephemeral Docker candidate build."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wheel_metadata(path: Path) -> tuple[str, str, list[str]]:
    with zipfile.ZipFile(path) as wheel:
        metadata_names = [
            name for name in wheel.namelist() if name.endswith(".dist-info/METADATA")
        ]
        wheel_names = [
            name for name in wheel.namelist() if name.endswith(".dist-info/WHEEL")
        ]
        if len(metadata_names) != 1 or len(wheel_names) != 1:
            raise RuntimeError(f"invalid wheel metadata layout: {path.name}")
        metadata = wheel.read(metadata_names[0]).decode("utf-8")
        wheel_text = wheel.read(wheel_names[0]).decode("utf-8")
    name_match = re.search(r"^Name:\s*(.+)$", metadata, re.MULTILINE)
    version_match = re.search(r"^Version:\s*(.+)$", metadata, re.MULTILINE)
    tags = re.findall(r"^Tag:\s*(.+)$", wheel_text, re.MULTILINE)
    if not name_match or not version_match or not tags:
        raise RuntimeError(f"wheel lacks Name, Version, or Tag: {path.name}")
    return name_match.group(1).strip(), version_match.group(1).strip(), sorted(tags)


def normalized_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def record(wheelhouse: Path, expected_openqp_version: str) -> dict[str, object]:
    wheels = sorted(wheelhouse.glob("*.whl"))
    if not wheels:
        raise RuntimeError(f"wheelhouse contains no wheels: {wheelhouse}")
    unexpected = sorted(
        path.name for path in wheelhouse.iterdir() if path.is_file() and path.suffix != ".whl"
    )
    # A previous manifest is the only allowed non-wheel file when re-running.
    unexpected = [name for name in unexpected if name != "wheelhouse-manifest.json"]
    if unexpected:
        raise RuntimeError(f"non-wheel dependency artifacts found: {unexpected}")

    distributions: dict[tuple[str, str], str] = {}
    entries = []
    openqp_count = 0
    for path in wheels:
        name, version, tags = wheel_metadata(path)
        key = (normalized_name(name), version)
        previous = distributions.get(key)
        if previous is not None:
            raise RuntimeError(
                f"duplicate wheel distribution/version: {previous}, {path.name}"
            )
        distributions[key] = path.name
        if key[0] == "openqp":
            openqp_count += 1
            if version != expected_openqp_version:
                raise RuntimeError(
                    f"OpenQP wheel version {version} != {expected_openqp_version}"
                )
        entries.append(
            {
                "filename": path.name,
                "name": name,
                "version": version,
                "size": path.stat().st_size,
                "sha256": sha256(path),
                "tags": tags,
            }
        )
    if openqp_count != 1:
        raise RuntimeError(f"expected exactly one OpenQP wheel, found {openqp_count}")

    return {
        "schema_version": 1,
        "wheels": entries,
        "publication_gate": {
            "ready": False,
            "reason": (
                "Dependency wheels were resolved during this candidate build. "
                "Public publication requires a reviewed checked-in hash lock or "
                "reuse of an already verified exact release-run wheelhouse."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheelhouse", type=Path)
    parser.add_argument("--expected-openqp-version", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = record(args.wheelhouse, args.expected_openqp_version)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Recorded {len(manifest['wheels'])} exact Docker wheel inputs")


if __name__ == "__main__":
    main()
