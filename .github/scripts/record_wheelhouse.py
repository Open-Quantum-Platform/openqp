#!/usr/bin/env python3
"""Record the exact wheels used for one ephemeral Docker candidate build."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path, PurePosixPath


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wheel_metadata(path: Path) -> tuple[str, str, list[str]]:
    with zipfile.ZipFile(path) as wheel:
        # Wheel metadata lives in one top-level ``*.dist-info`` directory.
        # Some distributions (notably setuptools) also bundle vendored wheels'
        # metadata below package directories; those entries do not describe the
        # containing wheel and must not make its layout appear ambiguous.
        def root_metadata(name: str, filename: str) -> bool:
            parts = PurePosixPath(name).parts
            return (
                len(parts) == 2
                and parts[0].endswith(".dist-info")
                and parts[1] == filename
            )

        metadata_names = [
            name for name in wheel.namelist() if root_metadata(name, "METADATA")
        ]
        wheel_names = [
            name for name in wheel.namelist() if root_metadata(name, "WHEEL")
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


def wheel_entries(
    wheelhouse: Path, *, allow_existing_manifest: bool = False
) -> list[dict[str, object]]:
    wheels = sorted(wheelhouse.glob("*.whl"))
    if not wheels:
        raise RuntimeError(f"wheelhouse contains no wheels: {wheelhouse}")
    unexpected = sorted(
        path.name for path in wheelhouse.iterdir() if path.is_file() and path.suffix != ".whl"
    )
    if allow_existing_manifest:
        # A previous manifest is the only allowed non-wheel file when re-running.
        unexpected = [name for name in unexpected if name != "wheelhouse-manifest.json"]
    if unexpected:
        raise RuntimeError(f"non-wheel dependency artifacts found: {unexpected}")

    distributions: dict[tuple[str, str], str] = {}
    entries = []
    for path in wheels:
        name, version, tags = wheel_metadata(path)
        key = (normalized_name(name), version)
        previous = distributions.get(key)
        if previous is not None:
            raise RuntimeError(
                f"duplicate wheel distribution/version: {previous}, {path.name}"
            )
        distributions[key] = path.name
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
    return entries


def record(
    wheelhouse: Path,
    expected_openqp_version: str,
    build_wheelhouse: Path,
) -> dict[str, object]:
    entries = wheel_entries(wheelhouse, allow_existing_manifest=True)
    openqp_entries = [
        entry for entry in entries if normalized_name(str(entry["name"])) == "openqp"
    ]
    if len(openqp_entries) == 1 and openqp_entries[0]["version"] != expected_openqp_version:
        raise RuntimeError(
            f"OpenQP wheel version {openqp_entries[0]['version']} != "
            f"{expected_openqp_version}"
        )
    openqp_count = len(openqp_entries)
    if openqp_count != 1:
        raise RuntimeError(f"expected exactly one OpenQP wheel, found {openqp_count}")

    build_entries = wheel_entries(build_wheelhouse)
    if any(normalized_name(str(entry["name"])) == "openqp" for entry in build_entries):
        raise RuntimeError("build-system wheelhouse unexpectedly contains OpenQP")

    return {
        "schema_version": 2,
        "wheels": entries,
        "build_wheels": build_entries,
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
    parser.add_argument("--build-wheelhouse", type=Path, required=True)
    parser.add_argument("--expected-openqp-version", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = record(
        args.wheelhouse,
        args.expected_openqp_version,
        args.build_wheelhouse,
    )
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"Recorded {len(manifest['wheels'])} runtime and "
        f"{len(manifest['build_wheels'])} build wheel inputs"
    )


if __name__ == "__main__":
    main()
