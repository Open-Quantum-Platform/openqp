"""Preserve PyInstaller's signed macOS layout when embedding a console sidecar."""
import re
import shutil
import subprocess
import sys
from pathlib import Path


def stage(root: Path) -> Path:
    backend = root / "backend"
    spec = root / ".cache/macos-backend.spec"
    spec.write_text((backend / "oqp-studio-backend.spec").read_text() +
                    '\napp = BUNDLE(coll, name="StudioBackend.app", '
                    'bundle_identifier="org.openqp.studio.backend")\n')
    subprocess.run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--distpath", "dist",
                    "--workpath", "build", str(spec)], cwd=backend, check=True)
    contents = backend / "dist/StudioBackend.app/Contents"
    binaries = root / "shell/src-tauri/binaries"
    for name in ("Frameworks", "Resources"):
        destination = binaries / ("backend-" + name.lower())
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(contents / name, destination, symlinks=True)
    return contents / "MacOS/oqp-studio-backend"


def verify(app: Path) -> dict:
    count = 0
    maximum = (0, 0)
    for path in app.rglob("*"):
        if not path.is_file() or path.is_symlink():
            continue
        with path.open("rb") as stream:
            magic = stream.read(4)
        if magic not in (b"\xcf\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"):
            continue
        metadata = subprocess.check_output(["otool", "-l", str(path)], text=True)
        versions = re.findall(r"\bminos (\d+)\.(\d+)", metadata)
        if not versions:
            versions = re.findall(r"cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version (\d+)\.(\d+)", metadata)
        if not versions:
            raise RuntimeError(f"Missing minimum macOS version: {path}")
        maximum = max(maximum, *(tuple(map(int, value)) for value in versions))
        if maximum > (15, 0):
            raise RuntimeError(f"Runtime requires macOS newer than 15: {path}")
        count += 1
    subprocess.run(["codesign", "--verify", "--deep", "--strict", str(app)], check=True)
    return {"mach_o_count": count, "minimum_macos": list(maximum), "signature": "ad-hoc-verified"}
