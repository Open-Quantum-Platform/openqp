"""Package a candidate from its private, prequalified engine artifact."""
import hashlib
import json
import os
import platform
import subprocess
import tarfile
import urllib.request
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    if os.environ.get("STUDIO_RELEASE_CANDIDATE") != "1":
        raise SystemExit("native packaging is explicitly selected, never an MR side effect")
    actual = {"Darwin": "macos-apple-silicon" if platform.machine() == "arm64" else "macos-intel",
              "Windows": "windows-x64", "Linux": "linux-x86_64"}[platform.system()]
    if os.environ.get("STUDIO_PACKAGE_PLATFORM") != actual:
        raise SystemExit("native packaging runner does not match the requested platform")
    config = json.loads((ROOT / "release-candidate.json").read_text())
    version = config["engine"]["commit"]
    base = (os.environ["CI_API_V4_URL"] + "/projects/" + os.environ["CI_PROJECT_ID"]
            + f"/packages/generic/studio-engine/{version}/{actual}")
    headers = {"JOB-TOKEN": os.environ["CI_JOB_TOKEN"]}
    cache = ROOT / ".cache"
    cache.mkdir(exist_ok=True)
    archive = cache / "qualified-engine.tar.gz"
    request = urllib.request.Request(base + ".sha256", headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        expected = response.read(1000).decode().split()[0]
    request = urllib.request.Request(base + ".tar.gz", headers=headers)
    digest = hashlib.sha256()
    with urllib.request.urlopen(request, timeout=600) as response, archive.open("wb") as out:
        while chunk := response.read(1 << 20):
            out.write(chunk)
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise SystemExit("private engine artifact checksum mismatch")
    with tarfile.open(archive) as stream:
        stream.extractall(cache / "qualified-engine", filter="data")
    build_env = cache / "package-venv"
    venv.EnvBuilder(with_pip=True).create(build_env)
    python = build_env / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    subprocess.run([str(python), str(ROOT / "tools/release/package.py"), "--engine",
                    str(cache / "qualified-engine/openqp")], check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
