"""Stage the scipy-openblas64 wheel as a plain ILP64 OpenBLAS install tree.

OpenQP's BLAS lookup asks pkg-config for the ``openblas64`` module (ILP64 only;
see validateOpenBLASPkgConfigInt64 in cmake/oqp_functions.cmake).  The
scipy-openblas64 wheel carries exactly that library -- the ILP64 OpenBLAS numpy
and scipy ship -- but publishes it under its own name and inside site-packages,
which differs per interpreter.  Copy it to a fixed prefix and write an
``openblas64.pc`` describing it, so every cibuildwheel interpreter and the
delvewheel repair step can use one stable path.

Usage:  python stage_openblas64_windows.py <prefix>
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    prefix = Path(argv[1])

    import scipy_openblas64

    src = Path(scipy_openblas64.__file__).parent
    version = scipy_openblas64.get_openblas_config().split()[1]

    for sub in ("bin", "lib", "include", "lib/pkgconfig"):
        (prefix / sub).mkdir(parents=True, exist_ok=True)

    copied = {"bin": [], "lib": [], "include": []}
    for dll in (src / "lib").glob("*.dll"):
        shutil.copy2(dll, prefix / "bin" / dll.name)
        copied["bin"].append(dll.name)
    # Some wheel layouts put the DLL next to the import library instead.
    for dll in src.glob("*.dll"):
        shutil.copy2(dll, prefix / "bin" / dll.name)
        copied["bin"].append(dll.name)
    for lib in (src / "lib").glob("*.lib"):
        shutil.copy2(lib, prefix / "lib" / lib.name)
        copied["lib"].append(lib.name)
    for header in (src / "include").glob("*"):
        if header.is_file():
            shutil.copy2(header, prefix / "include" / header.name)
            copied["include"].append(header.name)

    if not copied["bin"] or not copied["lib"]:
        raise SystemExit(
            f"scipy_openblas64 layout not understood under {src}: "
            f"found {copied}"
        )

    # The import library name is what -l must reference.
    stem = Path(copied["lib"][0]).stem

    pc = prefix / "lib" / "pkgconfig" / "openblas64.pc"
    pc.write_text(
        "\n".join(
            [
                f"prefix={prefix.as_posix()}",
                "libdir=${prefix}/lib",
                "includedir=${prefix}/include",
                "",
                "Name: openblas64",
                "Description: ILP64 OpenBLAS staged from the scipy-openblas64 wheel",
                f"Version: {version}",
                "Cflags: -I${includedir}",
                f"Libs: -L${{libdir}} -l{stem}",
                "",
            ]
        ),
        encoding="ascii",
    )

    print(f"staged scipy-openblas64 {version} into {prefix}")
    for kind, names in copied.items():
        print(f"  {kind}: {', '.join(sorted(names))}")
    print(f"  pkg-config: {pc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
