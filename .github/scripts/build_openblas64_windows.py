"""Build an ILP64 OpenBLAS for Windows with unrenamed symbols.

The Linux wheels build OpenBLAS from source with ``INTERFACE64=1`` and no
symbol suffix, so it exports plain ``dgemm_``/``dsytrf_``.  The obvious Windows
shortcut -- the ``scipy-openblas64`` wheel -- cannot replace that: it renames
every symbol (``scipy_`` prefix, ``64_`` suffix; its library is literally
``libscipy_openblas64_.lib``), which numpy and scipy reach through their own
macros but a plain ``DGEMV`` caller such as multicharge never can.

So build it here too, with the Intel toolchain rather than MinGW: OpenQP's
Windows wheels are compiled with ifx/icx, and a gfortran-built OpenBLAS would
drag the MinGW runtime into the wheel alongside Intel's.

Writes, under <prefix>:
  bin/, lib/, include/          the install tree
  lib/pkgconfig/openblas64.pc   the module cmake/oqp_functions.cmake looks for
  openblas-cache.cmake          BLAS_LIBRARIES/LAPACK_LIBRARIES for `cmake -C`,
                                which is how the DFT-D4 subprojects are told
                                where the BLAS is (external/CMakeLists.txt,
                                DFTD4_BLAS_ARGS)

Usage:  python build_openblas64_windows.py <prefix> [version]
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_VERSION = "0.3.30"


def run(cmd: list[str], **kw) -> None:
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def write_metadata(prefix: Path, version: str) -> None:
    """pkg-config module and CMake initial cache for the staged install."""
    libs = sorted((prefix / "lib").glob("*.lib"))
    if not libs:
        raise SystemExit(f"no import library under {prefix / 'lib'}")
    # openblas.lib, not one of the static/aux ones, when several are present.
    lib = next((p for p in libs if p.stem.lower().startswith("openblas")), libs[0])

    pc_dir = prefix / "lib" / "pkgconfig"
    pc_dir.mkdir(parents=True, exist_ok=True)
    (pc_dir / "openblas64.pc").write_text(
        "\n".join(
            [
                f"prefix={prefix.as_posix()}",
                "libdir=${prefix}/lib",
                "includedir=${prefix}/include",
                "",
                "Name: openblas64",
                "Description: ILP64 OpenBLAS built for OpenQP (unrenamed symbols)",
                f"Version: {version}",
                "Cflags: -I${includedir}",
                f"Libs: -L${{libdir}} -l{lib.stem}",
                "",
            ]
        ),
        encoding="ascii",
    )

    (prefix / "openblas-cache.cmake").write_text(
        "\n".join(
            [
                f'set(BLAS_LIBRARIES "{lib.as_posix()}" CACHE FILEPATH "")',
                f'set(LAPACK_LIBRARIES "{lib.as_posix()}" CACHE FILEPATH "")',
                "",
            ]
        ),
        encoding="ascii",
    )
    # Which spelling did it actually export?  ifx on Windows emits UPPERCASE
    # with no underscore by default, OpenBLAS normally uses the Unix `dgemm_`,
    # and the "openblas64" distributions add a `64_` suffix -- link failures
    # here are always one of the three, so say which is present.
    blob = lib.read_bytes()
    spellings = [b"dgemm_64_", b"dgemm_", b"DGEMM", b"dgemm"]
    found = [sp.decode() for sp in spellings if sp in blob]
    print(f"dgemm spellings present in {lib.name}: {found or '(none)'}")
    if not found:
        raise SystemExit(
            f"{lib} exports no recognisable dgemm symbol; the OpenBLAS build "
            "did not produce usable BLAS"
        )

    print(f"import library: {lib}")
    print(f"pkg-config:     {pc_dir / 'openblas64.pc'}")
    print(f"initial cache:  {prefix / 'openblas-cache.cmake'}")


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__)
        return 2
    prefix = Path(argv[1])
    version = argv[2] if len(argv) == 3 else DEFAULT_VERSION

    if list(prefix.glob("lib/*.lib")):
        print(f"OpenBLAS already present in {prefix}; reusing it")
        write_metadata(prefix, version)
        return 0

    src = Path(os.environ.get("RUNNER_TEMP", ".")) / "OpenBLAS"
    build = Path(os.environ.get("RUNNER_TEMP", ".")) / "OpenBLAS-build"
    if not src.exists():
        run(["git", "clone", "--depth", "1", "--branch", f"v{version}",
             "https://github.com/OpenMathLib/OpenBLAS.git", str(src)])
    shutil.rmtree(build, ignore_errors=True)

    cmplr = os.environ.get("CMPLR_ROOT", "").replace("\\", "/")
    if not cmplr:
        raise SystemExit("CMPLR_ROOT is not set; run the oneAPI setup first")

    run([
        "cmake", "-S", str(src), "-B", str(build), "-G", "Ninja",
        f"-DCMAKE_C_COMPILER={cmplr}/bin/icx.exe",
        f"-DCMAKE_Fortran_COMPILER={cmplr}/bin/ifx.exe",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={prefix.as_posix()}",
        # Static.  OpenBLAS's CMake shared build on Windows does not produce a
        # usable export table -- the import library resolves neither `DGEMM`
        # nor `dgemm_` -- and linking the archive into liboqp.dll sidesteps the
        # question entirely, with nothing extra for delvewheel to bundle.
        "-DBUILD_SHARED_LIBS=OFF",
        # The whole point: 8-byte BLAS integers, and NO symbol renaming --
        # OpenQP and the DFT-D4 stack call plain dgemm_/dsytrf_.
        "-DINTERFACE64=ON",
        "-DSYMBOLPREFIX=",
        "-DSYMBOLSUFFIX=",
        # Runtime CPU dispatch, as the Linux wheels do, so one wheel serves
        # every x86-64 machine -- but only over the kernels that build here.
        # OpenBLAS's AVX-512 micro-kernels are compiled without their target
        # feature flags under clang-cl ("always_inline function
        # '_mm512_setzero_ps' requires target feature 'avx512f'"), so the
        # dispatch list stops at AVX2.  Every machine from Sandy Bridge on is
        # still covered; only the AVX-512 fast paths are given up.
        "-DDYNAMIC_ARCH=ON",
        "-DDYNAMIC_LIST=NEHALEM;SANDYBRIDGE;HASWELL;ZEN",
        # OpenBLAS's own thread pool rather than OpenMP.  Its OpenMP check
        # rejects this toolchain outright -- icx advertises -Qiopenmp and ifx
        # -openmp, and it reads two spellings of Intel's single runtime as
        # "multiple OpenMP runtime libraries" -- and a pthread OpenBLAS is what
        # numpy and scipy ship on Windows anyway.  OpenQP drives its thread
        # count through openblas_set_num_threads (source/mathlib/
        # blas_thread_ctl.c), which works either way.
        "-DUSE_OPENMP=OFF",
        "-DUSE_THREAD=ON",
        "-DBUILD_TESTING=OFF",
        "-DBUILD_WITHOUT_LAPACK=OFF",
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        # OpenBLAS's x86-64 micro-kernels call _mm_prefetch(&x[0], ...) with a
        # float*.  On Unix the intrinsic takes const void*, but Windows'
        # winnt.h declares `CHAR CONST *`, and current clang (icx in clang-cl
        # mode) makes the mismatch an error rather than a warning.  The
        # conversion is harmless -- prefetch ignores the pointee type.
        "-DCMAKE_C_FLAGS=/clang:-Wno-incompatible-pointer-types",
        # OpenBLAS's C kernels export the Unix convention (`dgemm_`); build its
        # Fortran LAPACK the same way so the archive is internally consistent
        # and matches what OpenQP emits under /names:lowercase.
        "-DCMAKE_Fortran_FLAGS=/names:lowercase /assume:underscore",
    ])
    run(["cmake", "--build", str(build), "--parallel"])
    run(["cmake", "--install", str(build)])

    write_metadata(prefix, version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
