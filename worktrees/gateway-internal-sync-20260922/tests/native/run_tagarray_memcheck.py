#!/usr/bin/env python3
"""Run the installed OpenQP's TagArray C ABI under Valgrind (Linux CI gate).

Requires Valgrind >= 3.22, which detects aligned new/delete mismatches.
No copied dependency implementation, allocator replacement, or suppression is used.
An explicit --library supports negative controls against an older installed build.
"""
import argparse
import hashlib
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--output", type=Path, default=Path("tagarray-memory"))
    args = parser.parse_args()
    version = subprocess.check_output(["valgrind", "--version"], text=True).strip()
    numbers = re.search(r"(\d+)\.(\d+)", version)
    if not numbers or tuple(map(int, numbers.groups())) < (3, 22):
        raise RuntimeError(f"Valgrind >= 3.22 required for alignment checks: {version}")
    library = args.library
    if library is None:
        import oqp
        from oqp.runtime import library_path

        library = library_path(oqp.oqp_root, oqp.suffix)
    if library is None or not library.is_file():
        raise FileNotFoundError(f"Installed liboqp not found: {library}")
    library = library.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("tagarray_lifetime.c")
    executable = output / "tagarray_lifetime"
    compiler = shlex.split(os.environ.get("CC", "cc"))
    subprocess.run(compiler + ["-std=c11", "-Wall", "-Wextra", "-Werror", "-g",
                              str(source), "-ldl", "-o", str(executable)], check=True)
    compiler_version = subprocess.check_output(compiler + ["--version"], text=True).splitlines()[0]
    provenance = (f"{version}\ncompiler={compiler_version}\nlibrary={library}\n"
                  f"sha256={hashlib.sha256(library.read_bytes()).hexdigest()}\n")
    (output / "provenance.txt").write_text(provenance)
    print(provenance, flush=True)
    result = subprocess.run([
        "valgrind", "--tool=memcheck", "--track-origins=yes", "--keep-debuginfo=yes",
        "--leak-check=full", "--show-leak-kinds=definite,indirect",
        "--errors-for-leak-kinds=definite,indirect", "--error-exitcode=97",
        f"--log-file={output / 'valgrind.log'}", str(executable), str(library),
    ], timeout=120, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output / "probe.log").write_text(result.stdout)
    print(result.stdout, end="")
    print((output / "valgrind.log").read_text())
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
