"""Public, version-bound Studio manifests and verified binary downloads."""
from __future__ import annotations

import hashlib
import json
import re
import urllib.request
from pathlib import Path

from . import __version__, network

REPOSITORY = "Open-Quantum-Platform/oqp-studio-releases"
PAGE = "https://open-quantum-platform.github.io/openqp-docs/studio/download/"
API = f"https://api.github.com/repos/{REPOSITORY}/releases"
ASSETS = f"https://github.com/{REPOSITORY}/releases/download/"


def version_key(value: str) -> tuple:
    match = re.fullmatch(r"v?(\d+)\.(\d+)\.(\d+)(?:-([A-Za-z0-9.]+))?", value)
    if not match:
        raise ValueError("invalid Studio version")
    core = tuple(int(match[i]) for i in (1, 2, 3))
    prerelease = match[4]
    parts = tuple((0, int(p)) if p.isdigit() else (1, p)
                  for p in prerelease.split(".")) if prerelease else ()
    return core, prerelease is None, parts


def read_json(url: str) -> dict | list:
    request = urllib.request.Request(url, headers={"User-Agent": "oqp-studio"})
    with urllib.request.urlopen(request, timeout=15, context=network.context()) as response:
        return json.loads(response.read(4_000_001))


def validate(manifest: dict, tag: str) -> dict:
    if manifest.get("schema_version") != 1 or tag != "v" + manifest.get("version", ""):
        raise ValueError("invalid Studio release identity")
    if not re.fullmatch(r"v\d+\.\d+\.\d+(?:-[A-Za-z0-9.]+)?", tag):
        raise ValueError("invalid Studio tag")
    if manifest.get("channel") not in ("stable", "preview"):
        raise ValueError("invalid release channel")
    for sha in (manifest.get("studio_commit"), manifest.get("engine", {}).get("commit")):
        if not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{40}", sha):
            raise ValueError("missing exact source commit")
    names = set()
    for asset in manifest.get("assets", []):
        name = asset.get("name", "")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]+", name) or name in names:
            raise ValueError("invalid or duplicate asset name")
        names.add(name)
        if asset.get("url") != ASSETS + tag + "/" + name:
            raise ValueError("asset must belong to this public Studio release")
        if not re.fullmatch(r"[0-9a-f]{64}", asset.get("sha256", "")):
            raise ValueError("asset checksum missing")
        if not isinstance(asset.get("size"), int) or asset["size"] <= 0:
            raise ValueError("asset size missing")
        if asset.get("kind") not in ("installer", "engine"):
            raise ValueError("invalid asset kind")
        if asset.get("variant") not in ("slim", "with-engine", "engine"):
            raise ValueError("invalid installer variant")
    if not names:
        raise ValueError("release contains no qualified assets")
    return manifest


def manifest_for(release: dict) -> dict:
    tag = release["tag_name"]
    url = ASSETS + tag + "/manifest.json"
    if release.get("draft") or not any(
        a.get("browser_download_url") == url for a in release.get("assets", [])
    ):
        raise ValueError("release has no public manifest")
    manifest = validate(read_json(url), tag)
    if bool(release.get("prerelease")) != (manifest["channel"] == "preview"):
        raise ValueError("release channel disagrees with manifest")
    return manifest


def installed_manifest() -> dict:
    """The current app's engine, never the engine from a newer Studio release."""
    return manifest_for(read_json(API + "/tags/v" + __version__))


def channel() -> str:
    try:
        value = json.loads(network.settings_path().with_name("updates.json").read_text())
        return "preview" if value.get("channel") == "preview" else "stable"
    except (OSError, ValueError):
        return "stable"


def set_channel(value: str) -> None:
    if value not in ("stable", "preview"):
        raise ValueError("choose stable or preview")
    path = network.settings_path().with_name("updates.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"channel": value}))


def latest() -> tuple[dict, dict]:
    if channel() == "stable":
        release = read_json(API + "/latest")
    else:
        # GitHub lists releases newest first; Preview is explicitly opt-in.
        release = max((r for r in read_json(API + "?per_page=100") if not r.get("draft")),
                      key=lambda r: version_key(r["tag_name"]))
    return release, manifest_for(release)


def download(asset: dict, target: Path, progress=None) -> Path:
    """Only return a file after its declared size and SHA-256 match."""
    target.parent.mkdir(parents=True, exist_ok=True)
    digest, done = hashlib.sha256(), 0
    try:
        with urllib.request.urlopen(asset["url"], timeout=120,
                                    context=network.context()) as response, target.open("wb") as out:
            while chunk := response.read(1 << 20):
                done += len(chunk)
                if done > asset["size"]:
                    raise ValueError("download exceeds manifest size")
                digest.update(chunk)
                out.write(chunk)
                if progress:
                    progress(done, asset["size"])
        if done != asset["size"] or digest.hexdigest() != asset["sha256"]:
            raise ValueError("download checksum or size mismatch")
        return target
    except Exception:
        target.unlink(missing_ok=True)
        raise


def build_identity() -> dict:
    path = Path(__file__).with_name("build-identity.json")
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {"version": __version__, "channel": "development"}
