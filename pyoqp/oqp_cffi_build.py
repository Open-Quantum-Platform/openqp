#!/usr/bin/env python3

from cffi import FFI
import argparse
import re
from glob import glob


def find_file_in_dirs(dirs, name):
    for d in dirs:
        try_file = f"{d}/{name}"
        if glob(try_file):
            return try_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-I", "--include", help="include dir", action="append", nargs="*"
    )
    parser.add_argument(
        "-L", "--libdir", help="required libraries dir", action="append", nargs="*"
    )
    parser.add_argument("--tmpdir", type=str, help="include dir", default=".")
    parser.add_argument(
        "--without-dftd4",
        action="store_true",
        help="liboqp was built with ENABLE_DFTD4=OFF: drop the DFT-D4 entry "
             "points from the binding, since the symbols do not exist",
    )
    parser.add_argument("-l", help="shared libraries", action="append", nargs="*")
    parser.add_argument(
        "other_libraries", help="other libraries", action="append", nargs="*"
    )
    args = parser.parse_args()

    print(args.include)
    paths = [p[0] for p in args.include]
    oqph = find_file_in_dirs(paths, "oqp.h")
    # with open('/'.join([args.include,'oqp.h'])) as oqp_header:
    with open(oqph) as oqp_header:
        defs = oqp_header.read().replace("#include", "//#include")

    if args.without_dftd4:
        # oqp.h declares every entry point unconditionally, but an
        # ENABLE_DFTD4=OFF liboqp does not define these two.  API mode compiles
        # a real reference for each declaration, so leaving them in fails at
        # link time -- or, worse, at import with an undefined symbol.  The ABI
        # path tolerates their absence separately (see oqp/__init__.py).
        defs = re.sub(
            r"^void oqp_dftd4_disp(_v2)?\([^;]*?\);\s*$",
            "",
            defs,
            flags=re.MULTILINE | re.DOTALL,
        )
        for symbol in ("oqp_dftd4_disp(", "oqp_dftd4_disp_v2("):
            assert symbol not in defs, f"failed to drop {symbol} from the cdefs"

    libs_shared = [lib[0] for lib in args.l] if args.l else []
    libs_other = (
        [lib[0] for lib in args.other_libraries] if args.other_libraries != [[]] else []
    )
    ffibuilder = FFI()

    ffibuilder.cdef(defs)

    ffibuilder.set_source(
        "_oqp",
        """
    #include "oqp.h"
""",
        include_dirs=[incdir[0] for incdir in args.include],
        # libraries = ['common', 'oqp.io', 'oqp.mathlib', 'oqp.precision'],
        libraries=libs_shared,
        library_dirs=[libdir[0] for libdir in args.libdir],
        extra_objects=libs_other,
    )

    ffibuilder.compile(tmpdir=args.tmpdir, verbose=True)
