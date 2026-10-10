"""Freeze a verified, isolated OpenQP installation without rebuilding it."""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads((ROOT / "release-candidate.json").read_text())
    sha = subprocess.check_output(["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True).strip()
    if sha != config["engine"]["commit"]:
        raise SystemExit("engine checkout is not the selected candidate")
    import oqp
    from oqp import runtime

    pkg = Path(oqp.__file__).resolve().parent
    root, suffix = runtime.resolve_oqp_root()
    library = Path(runtime.library_path(root, suffix)).resolve()
    if not pkg.is_relative_to(Path(sys.prefix)) or not library.is_relative_to(Path(sys.prefix)):
        raise SystemExit("engine must come from this isolated environment")
    if args.output.exists():
        raise SystemExit("refusing to replace an existing frozen engine")
    work = args.output.parent / "freeze-work"
    work.mkdir(parents=True, exist_ok=True)
    entry = work / "openqp_entry.py"
    entry.write_text("import sys\nfrom oqp.pyoqp import main\nif __name__ == '__main__':\n    sys.exit(main())\n")
    sep = ";" if os.name == "nt" else ":"
    binaries = []
    for path in (pkg / "lib").iterdir():
        if path.is_file() and any(s in path.name for s in (".dylib", ".so", ".dll")):
            binaries += ["--add-binary", f"{path}{sep}oqp/lib"]
    if os.name == "nt":
        # PyInstaller cannot infer DLLs dynamically loaded through cffi.
        folders = [Path(sys.prefix) / "Library/bin", Path(sys.prefix) / "Lib/site-packages/mkl"]
        if os.environ.get("MKLROOT"):
            folders += [Path(os.environ["MKLROOT"]) / "bin",
                        Path(os.environ["MKLROOT"]) / "redist/intel64"]
        for folder in folders:
            for path in folder.glob("*.dll"):
                binaries += ["--add-binary", f"{path}{sep}."]
    elif sys.platform.startswith("linux"):
        for pattern in ("libmkl*.so*", "libiomp5.so*"):
            for path in (Path(sys.prefix) / "lib").glob(pattern):
                binaries += ["--add-binary", f"{path}{sep}."]
    command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
               "--name", "openqp", "--distpath", str(args.output.parent),
               "--workpath", str(work / "build"), "--specpath", str(work),
               "--collect-all", "oqp", "--collect-all", "basis_set_exchange",
               "--collect-submodules", "scipy", "--collect-submodules", "numpy",
               "--recursive-copy-metadata", "basis_set_exchange", "--copy-metadata", "OpenQP",
               *binaries, str(entry)]
    subprocess.run(command, check=True)
    built = args.output.parent / "openqp"
    if built != args.output:
        built.rename(args.output)
    shutil.copytree(pkg / "share", args.output / "share", dirs_exist_ok=True)
    shutil.copytree(args.source / "examples/NMR", args.output / "examples/NMR")
    (args.output / "engine-identity.json").write_text(json.dumps(config["engine"], indent=2) + "\n")
    (args.output / "README.txt").write_text(
        f"OpenQP version : {config['engine']['version']}\nOpenQP commit : {sha}\n"
        f"Architecture : {platform.machine()}\n"
        "Qualified gateway snapshot for OQP Studio; version is distinct from the engine release tag.\n"
    )


if __name__ == "__main__":
    main()
