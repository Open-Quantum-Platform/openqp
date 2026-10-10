"""Install candidate-local Node and Rust; never change the host's toolchains."""
from __future__ import annotations

import hashlib
import os
import platform
import subprocess
import tarfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def environment() -> dict:
    cache = ROOT / ".cache/release-tools"
    cache.mkdir(parents=True, exist_ok=True)
    system = {"Darwin": "darwin", "Linux": "linux", "Windows": "win"}[platform.system()]
    arch = "arm64" if platform.machine().lower() in ("arm64", "aarch64") else "x64"
    version = "v24.19.0"
    name = f"node-{version}-{system}-{arch}"
    ext = ".zip" if system == "win" else ".tar.gz"
    archive = cache / (name + ext)
    node_bin = cache / name if system == "win" else cache / name / "bin"
    if not (node_bin / ("node.exe" if system == "win" else "node")).exists():
        base = f"https://nodejs.org/dist/{version}/"
        with urllib.request.urlopen(base + "SHASUMS256.txt", timeout=60) as response:
            checksums = dict(line.split()[::-1] for line in response.read().decode().splitlines())
        urllib.request.urlretrieve(base + archive.name, archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != checksums[archive.name]:
            raise RuntimeError("Node archive checksum mismatch")
        if ext == ".zip":
            with zipfile.ZipFile(archive) as stream:
                stream.extractall(cache)
        else:
            with tarfile.open(archive) as stream:
                stream.extractall(cache, filter="data")
    triple = {
        ("darwin", "arm64"): "aarch64-apple-darwin",
        ("darwin", "x64"): "x86_64-apple-darwin",
        ("linux", "x64"): "x86_64-unknown-linux-gnu",
        ("win", "x64"): "x86_64-pc-windows-msvc",
    }[(system, arch)]
    cargo_home = cache / "cargo"
    env = {**os.environ, "CARGO_HOME": str(cargo_home), "RUSTUP_HOME": str(cache / "rustup"),
           "PATH": str(node_bin) + os.pathsep + str(cargo_home / "bin") + os.pathsep + os.environ["PATH"],
           "CARGO_BUILD_JOBS": os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4"),
           "npm_config_cache": str(cache / "npm")}
    if not (cargo_home / "bin" / ("cargo.exe" if system == "win" else "cargo")).exists():
        binary = "rustup-init.exe" if system == "win" else "rustup-init"
        url = f"https://static.rust-lang.org/rustup/dist/{triple}/{binary}"
        installer = cache / binary
        urllib.request.urlretrieve(url, installer)
        with urllib.request.urlopen(url + ".sha256", timeout=60) as response:
            checksum = response.read().decode().split()[0]
        if hashlib.sha256(installer.read_bytes()).hexdigest() != checksum:
            raise RuntimeError("Rust bootstrap checksum mismatch")
        installer.chmod(0o700)
        subprocess.run([str(installer), "-y", "--no-modify-path", "--profile", "minimal",
                        "--default-toolchain", "1.94.0"], env=env, check=True)
    return env
