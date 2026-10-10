"""Build native Studio installers, retaining exact artifact hashes and provenance.

Invoke inside an admitted host build lane with an already-qualified frozen engine.
This command never uploads files or creates a public release.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

from bootstrap import environment
from macos_bundle import stage as stage_macos
from macos_bundle import verify as verify_macos
from provenance import source_commits
from qualify_engine import qualify

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256")
    return digest.hexdigest()


def check_sidecar(binary, env, identity):
    request = {"id": 1, "request": {"method": "GET", "path": "/api/health", "headers": {}, "body": ""}}
    result = subprocess.run([str(binary), "--stdio"], input=json.dumps(request) + "\n",
                            text=True, capture_output=True, timeout=120, env=env, check=True)
    replies = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{"id"')]
    reply = next(r for r in replies if r["id"] == 1)["result"]
    assert reply["status"] == 200
    health = json.loads(base64.b64decode(reply["body"]))
    assert health["version"] == identity["version"]
    assert health["build"]["studio_commit"] == identity["studio_commit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, required=True)
    args = parser.parse_args()
    engine = args.engine.resolve()
    identity = json.loads((ROOT / "release-candidate.json").read_text())
    identity["studio_commit"], packaging_commit = source_commits(
        ROOT, os.environ.get("STUDIO_APPLICATION_COMMIT"))
    if json.loads((engine / "engine-identity.json").read_text())["commit"] != identity["engine"]["commit"]:
        raise SystemExit("frozen engine belongs to a different candidate")
    system = platform.system()
    arch = "arm64" if platform.machine().lower() in ("arm64", "aarch64") else "x86_64"
    token = {("Darwin", "arm64"): "macos-apple-silicon", ("Darwin", "x86_64"): "macos-intel",
             ("Linux", "x86_64"): "linux-x86_64", ("Windows", "x86_64"): "windows-x64"}[(system, arch)]
    out = ROOT / "release" / token / identity["studio_commit"]
    out.mkdir(parents=True, exist_ok=False)
    qualification = qualify(engine, ROOT / ".cache" / f"qualification-{token}-{identity['studio_commit'][:12]}")
    env = environment()
    env.update(OQP_STUDIO_CONFIG=str(ROOT / ".cache/package-smoke/network.json"), CI="true")
    if system == "Darwin":
        env["MACOSX_DEPLOYMENT_TARGET"] = "15.0"
    npm = shutil.which("npm", path=env["PATH"])
    npx = shutil.which("npx", path=env["PATH"])
    subprocess.run([npm, "ci"], cwd=ROOT / "frontend", env=env, check=True)
    subprocess.run([npm, "audit", "--audit-level=moderate"], cwd=ROOT / "frontend", env=env, check=True)
    subprocess.run([npm, "run", "build"], cwd=ROOT / "frontend", env=env, check=True)
    if list((ROOT / "frontend/dist").rglob("*.map")):
        raise SystemExit("public payload must not contain source maps")
    (ROOT / "backend/oqp_studio/build-identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    subprocess.run([sys.executable, "-m", "pip", "install", "-e", str(ROOT / "backend"),
                    "rdkit", "pyinstaller==6.19.0"], env=env, check=True)
    subprocess.run([sys.executable, "build_binary.py"], cwd=ROOT / "backend", env=env, check=True)
    cargo = shutil.which("rustc", path=env["PATH"])
    triple = next(line.split(": ")[1] for line in subprocess.check_output([cargo, "-Vv"], env=env, text=True).splitlines() if line.startswith("host:"))
    tauri = ROOT / "shell/src-tauri"
    binaries = tauri / "binaries"
    binaries.mkdir(exist_ok=True)
    ext = ".exe" if system == "Windows" else ""
    source = ROOT / "backend/dist" / ("oqp-studio-backend" + ext)
    if system == "Darwin":
        source = stage_macos(ROOT)
    check_sidecar(source, env, identity)
    shutil.copy2(source, binaries / f"oqp-studio-backend-{triple}{ext}")
    baseline = json.loads((tauri / "tauri.conf.json").read_text())
    assets = []
    native_checks = {}
    for variant in ("slim", "with-engine"):
        config = json.loads(json.dumps(baseline))
        config["bundle"]["resources"] = []
        if system != "Darwin":
            config["bundle"].pop("macOS", None)
        else:
            config["bundle"]["macOS"]["minimumSystemVersion"] = "15.0"
            config["bundle"]["macOS"]["signingIdentity"] = "-"
            config["bundle"]["macOS"]["entitlements"] = str(ROOT / "tools/release/macos-entitlements.plist")
            config["bundle"]["macOS"]["files"] = {
                "Frameworks": "binaries/backend-frameworks",
                "Resources": "binaries/backend-resources",
            }
        if variant == "with-engine":
            if (tauri / "engine").exists():
                shutil.rmtree(tauri / "engine")
            # Preserve qualified Unix framework/runtime links. Windows freezes
            # contain ordinary files and must not require symlink privileges.
            shutil.copytree(engine, tauri / "engine", symlinks=system != "Windows")
            config["bundle"]["resources"] = ["engine/"]
        # NSIS supports both Windows variants without MSI's service-account ICE
        # validation dependency. Do not elevate the CI runner to build an MSI.
        bundles = {"Darwin": "app,dmg", "Linux": "deb,appimage", "Windows": "nsis"}[system]
        if variant == "with-engine":
            bundles = {"Darwin": "app,dmg", "Linux": "deb", "Windows": "nsis"}[system]
        override = ROOT / ".cache" / f"tauri-{variant}.json"
        override.write_text(json.dumps(config, indent=2))
        # Remove only this candidate's prior packaging output, never the Rust cache.
        bundle_dir = tauri / "target/release/bundle"
        if bundle_dir.exists():
            shutil.rmtree(bundle_dir)
        subprocess.run([npx, "--yes", "@tauri-apps/cli@2.11.5", "build", "--config", str(override),
                        "--bundles", bundles], cwd=tauri, env=env, check=True)
        stem = f"OQP-Studio-{identity['version']}-{token}" + ("-with-engine" if variant == "with-engine" else "")
        outputs = []
        if system == "Darwin":
            app = bundle_dir / "macos/OQP Studio.app"
            native_checks[variant] = verify_macos(app)
            check_sidecar(app / "Contents/MacOS/oqp-studio-backend", env, identity)
            archive = out / (stem + ".app.tar.gz")
            with tarfile.open(archive, "w:gz") as stream:
                stream.add(app, arcname=app.name)
            outputs.append(archive)
        for extension in (".dmg", ".AppImage", ".deb", ".msi", "-setup.exe"):
            for path in bundle_dir.glob("*/*" + extension):
                target = out / (stem + extension)
                shutil.copy2(path, target)
                outputs.append(target)
        if not outputs:
            raise SystemExit("native bundler produced no installer")
        for path in outputs:
            assets.append({"name": path.name, "size": path.stat().st_size, "sha256": sha256(path),
                           "platform": system.lower(), "architecture": arch,
                           "variant": variant, "kind": "installer"})
    # Unix tar preserves executable bits and symlinks; macOS/Windows use zip to
    # retain the published archive names expected by engine selection.
    engine_token = {"Darwin": "macos", "Linux": "linux", "Windows": "windows"}[system]
    name = f"openqp-{identity['engine']['version']}-{engine_token}-{arch}"
    if system == "Linux":
        archive = out / (name + ".tar.gz")
        with tarfile.open(archive, "w:gz") as stream:
            for path in engine.iterdir():
                stream.add(path, arcname=path.name)
    elif system == "Darwin":
        archive = out / (name + ".zip")
        subprocess.run(["/usr/bin/zip", "-qry", str(archive), "."], cwd=engine, check=True)
    else:
        archive = Path(shutil.make_archive(str(out / name), "zip", engine))
    assets.append({"name": archive.name, "size": archive.stat().st_size, "sha256": sha256(archive),
                   "platform": system.lower(), "architecture": arch, "variant": "engine", "kind": "engine"})
    record = {**identity, "packaging_commit": packaging_commit,
              "assets": assets, "engine_qualification": qualification,
              "sidecar_health": True, "platform": token, "native_checks": native_checks}
    (out / "qualification.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
