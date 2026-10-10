"""Collect qualified native artifacts into an immutable public release manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

PLATFORMS = {"macos-apple-silicon", "macos-intel", "windows-x64", "linux-x86_64"}
REPOSITORY = "Open-Quantum-Platform/oqp-studio-releases"


def assemble(records: list[Path], output: Path) -> dict:
    if output.exists():
        raise ValueError("refusing to replace an existing release payload")
    data = [(path, json.loads(path.read_text())) for path in records]
    if {d["platform"] for _, d in data} != PLATFORMS or len(data) != len(PLATFORMS):
        raise ValueError("all four native architectures must be qualified exactly once")
    identity = {key: data[0][1][key] for key in ("version", "channel", "studio_commit", "engine")}
    assets = []
    sources = {}
    packaging_commits = {}
    for record, candidate in data:
        packaging_commits[candidate["platform"]] = candidate.get("packaging_commit", candidate["studio_commit"])
        if any(candidate[key] != value for key, value in identity.items()):
            raise ValueError("platforms used different Studio/engine source or channel")
        if candidate["platform"].startswith("macos"):
            checks = candidate.get("native_checks", {})
            if any(checks.get(v, {}).get("signature") != "ad-hoc-verified"
                   for v in ("slim", "with-engine")):
                raise ValueError("both Mac app variants must pass signature validation")
        qualification = candidate["engine_qualification"]
        if not qualification.get("passed") or not candidate.get("sidecar_health"):
            raise ValueError("native qualification did not pass")
        if qualification["engine"]["commit"] != identity["engine"]["commit"]:
            raise ValueError("qualification used a different engine")
        if {a["variant"] for a in candidate["assets"]} != {"slim", "with-engine", "engine"}:
            raise ValueError("each architecture needs both installer variants and the engine")
        for asset in candidate["assets"]:
            name = asset["name"]
            if Path(name).name != name or name in sources:
                raise ValueError("unsafe or duplicate asset")
            # Never allow a Python source distribution/wheel or source archive.
            if not name.startswith(("OQP-Studio-", "openqp-")) or not name.endswith(
                (".dmg", ".app.tar.gz", ".deb", ".AppImage", ".msi", "-setup.exe", ".zip", ".tar.gz")
            ):
                raise ValueError("not an allowed binary distribution asset")
            path = record.parent / name
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if path.stat().st_size != asset["size"] or digest != asset["sha256"]:
                raise ValueError("artifact changed after native qualification")
            sources[name] = path
            assets.append({**asset, "url": f"https://github.com/{REPOSITORY}/releases/download/v{identity['version']}/{name}"})
    output.mkdir(parents=True)
    for name, source in sources.items():
        shutil.copy2(source, output / name)
    manifest = {"schema_version": 1, **identity, "packaging_commits": packaging_commits,
                "assets": sorted(assets, key=lambda a: a["name"])}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = [f"{a['sha256']}  {a['name']}" for a in manifest["assets"]]
    lines.append(f"{hashlib.sha256((output / 'manifest.json').read_bytes()).hexdigest()}  manifest.json")
    (output / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("records", nargs="+", type=Path)
    args = parser.parse_args()
    assemble(args.records, args.output)
