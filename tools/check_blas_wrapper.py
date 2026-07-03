#!/usr/bin/env python3
"""
Enforce PR RULE 1 (portable, source-level): every BLAS/LAPACK routine referenced
from OpenQP source must go through the wrapper layer -- i.e. it must be registered
in OQP_ACCELERATE_ILP64_SYMS (cmake/oqp_functions.cmake), the single list that (a)
the oqp_*_i64 wrappers are built from and (b) the macOS Accelerate link-time alias
interposition uses.

WHY THIS IS THE RIGHT INVARIANT (see commit ad43fae1, the dgels bug):
  OpenQP compiles the whole tree with -fdefault-integer-8 and BLA_SIZEOF_INTEGER=8,
  so on Linux ILP64 backends (MKL ILP64 / openblas64) a bare `call dgemm(...)` binds
  the genuine 8-byte symbol and is fine. The place a raw call breaks is macOS
  Accelerate, whose *classic* Fortran symbols are LP64 (4-byte): every referenced
  name must be interposed onto its $NEWLAPACK$ILP64 variant via the alias list, and
  the wrapper layer (lapack_wrap/blas_wrap + the oqp_linalg renames) is generated in
  lockstep with that same list. The dgels bug was exactly a routine that was CALLED
  but MISSING from OQP_ACCELERATE_ILP64_SYMS: it had no wrapper and no alias, so the
  macOS build's check_accelerate_aliases gate failed and OQP_BLAS_INT=4 builds passed
  4-byte integers to an ILP64 routine.

  check_accelerate_aliases.cmake already catches this -- but only on macOS and only
  at link time. This script is the PORTABLE gate: it runs on any PR on any OS by
  static inspection and names the offending routine + call site before the build.

Usage: check_blas_wrapper.py [REPO_ROOT]   (default: this script's own repo)
Exit: 0 clean, 1 violations, 2 setup error.
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "source")
CMAKE = os.path.join(REPO, "cmake", "oqp_functions.cmake")

# Files where raw BLAS/LAPACK references are correct BY DESIGN: the wrapper bodies
# deliberately call the raw symbol after converting integer args to blas_int.
WRAPPER_FILES = {
    os.path.normpath(os.path.join(SRC, "mathlib", "lapack_wrap.F90")),
    os.path.normpath(os.path.join(SRC, "mathlib", "blas_wrap.F90")),
}

# Canonical BLAS (L1/L2/L3) + core LAPACK routine names (no trailing underscore).
# This is what tells a genuine BLAS/LAPACK reference apart from an ordinary
# internal `call some_subroutine(`. Extend if a new routine family is adopted.
_BLAS1 = """
 srotg srotmg srot srotm sswap sscal scopy saxpy sdot sdsdot snrm2 sasum isamax
 drotg drotmg drot drotm dswap dscal dcopy daxpy ddot dsdot dnrm2 dasum idamax
 crotg csrot cswap cscal csscal ccopy caxpy cdotu cdotc scnrm2 scasum icamax
 zrotg zdrot zswap zscal zdscal zcopy zaxpy zdotu zdotc dznrm2 dzasum izamax
"""
_BLAS2 = """
 sgemv sgbmv ssymv ssbmv sspmv strmv stbmv stpmv strsv stbsv stpsv sger ssyr
 sspr ssyr2 sspr2
 dgemv dgbmv dsymv dsbmv dspmv dtrmv dtbmv dtpmv dtrsv dtbsv dtpsv dger dsyr
 dspr dsyr2 dspr2
 cgemv cgbmv chemv chbmv chpmv ctrmv ctbmv ctpmv ctrsv ctbsv ctpsv cgeru cgerc
 cher chpr cher2 chpr2
 zgemv zgbmv zhemv zhbmv zhpmv ztrmv ztbmv ztpmv ztrsv ztbsv ztpsv zgeru zgerc
 zher zhpr zher2 zhpr2
"""
_BLAS3 = """
 sgemm ssymm ssyrk ssyr2k strmm strsm
 dgemm dsymm dsyrk dsyr2k dtrmm dtrsm
 cgemm csymm chemm csyrk cher2k csyr2k cherk ctrmm ctrsm
 zgemm zsymm zhemm zsyrk zher2k zsyr2k zherk ztrmm ztrsm
"""
_LAPACK = """
 sgesv sgetrf sgetrs sgetri sgels sgesvd ssyev ssyevd sspev sspevx ssytrf
 ssytri ssytrs ssysv sgeqrf sorgqr sormqr spotrf spotri strtri stpttr strttp
 sgglse slamch
 dgesv dgetrf dgetrs dgetri dgels dgesvd dsyev dsyevd dspev dspevx dsytrf
 dsytri dsytrs dsysv dgeqrf dorgqr dormqr dpotrf dpotri dtrtri dtpttr dtrttp
 dgglse dlamch
 cgesv cgetrf cgetrs cgetri cgels cgesvd cheev cheevd chpev chpevx chetrf
 chetri chetrs chesv cgeqrf cungqr cunmqr cpotrf cpotri ctrtri ctpttr ctrttp
 cgglse
 zgesv zgetrf zgetrs zgetri zgels zgesvd zheev zheevd zhpev zhpevx zhetrf
 zhetri zhetrs zhesv zgeqrf zungqr zunmqr zpotrf zpotri ztrtri ztpttr ztrttp
 zgglse
 xerbla lsame
"""
CANON = set((_BLAS1 + _BLAS2 + _BLAS3 + _LAPACK).split())

# CANON is the *known* set, but it cannot list every LAPACK routine (there are
# hundreds). So a second, list-free heuristic catches a genuinely new routine:
# a `call NAME(` whose NAME has the LAPACK shape [precision][matrix-type][op],
# is NOT defined anywhere in the scanned OpenQP source (=> it is an external
# symbol resolved from the BLAS/LAPACK library), and is NOT registered. Anchoring
# on the real LAPACK matrix-type digraphs (ge, sy, he, po, tr, ...) -- not just
# "[sdcz] + letters" -- is what keeps ordinary externals like `call solver(` or
# `call second(` from being mistaken for LAPACK, while still flagging a new
# `call dgesdd(` / `call dpotrs(` / `call dsyevr(` that CANON does not list.
# (BLAS L1/L2/L3, which do not follow this scheme, are fully covered by CANON.)
# Restricting this to `call NAME(` (never bare function-style refs, which produce
# fragment false positives on names like `dftd4_calc`) and excluding
# OpenQP-defined procedures completes the false-positive guard.
BLAS_CALL_PAT = re.compile(
    r"^[sdcz](?:bd|di|gb|ge|gg|gt|hb|he|hg|hp|hs|la|op|or|pb|po|pp|pt|"
    r"sb|sp|st|sy|tb|tg|tp|tr|tz|un|up)[a-z0-9]{2,4}$")
# Belt-and-suspenders: [sdcz]-initial subroutines that survive the digraph anchor
# yet are not BLAS/LAPACK externals (none known today; kept as an explicit escape).
INTRINSIC_SUB = set()

DEF_RX = re.compile(r"\b(?:subroutine|function)\s+([a-z][a-z0-9_]*)", re.I)


def scan_files():
    """Fortran compiled into liboqp: all of source/, plus the tests/fortran/*.F90
    bind(C) self-test harnesses that source/CMakeLists.txt appends to the library
    target (they ship inside liboqp, so a raw call there evades the macOS alias
    gate exactly like one in source/). Wrapper bodies are excluded by the caller."""
    files = []
    for root, _, names in os.walk(SRC):
        for fn in names:
            if fn.endswith((".F90", ".f90")):
                files.append(os.path.normpath(os.path.join(root, fn)))
    cml = os.path.join(SRC, "CMakeLists.txt")
    try:
        with open(cml, encoding="utf-8") as fh:
            cml_txt = fh.read()
    except OSError:
        cml_txt = ""
    for rel in re.findall(r"tests/fortran/([A-Za-z0-9_]+\.[Ff]90)", cml_txt):
        p = os.path.normpath(os.path.join(REPO, "tests", "fortran", rel))
        if os.path.exists(p):
            files.append(p)
    return files


def defined_procedures(files):
    """Names of every subroutine/function DEFINED in the scanned source. A
    `call NAME(` to a name in this set is an internal routine, never a BLAS/LAPACK
    external, so the shape heuristic must not flag it."""
    defined = set()
    for path in files:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                for m in DEF_RX.finditer(strip_comment(line)):
                    defined.add(m.group(1).lower())
    return defined


def registered_syms():
    """Names in OQP_ACCELERATE_ILP64_SYMS (the wrapper/alias source of truth)."""
    with open(CMAKE, encoding="utf-8") as fh:
        txt = fh.read()
    m = re.search(r"set\(OQP_ACCELERATE_ILP64_SYMS(.*?)CACHE INTERNAL", txt, re.S)
    if not m:
        return None
    return set(w.lower() for w in re.findall(r"\b[a-z][a-z0-9]+\b", m.group(1)))


def strip_comment(line):
    """Drop a trailing free-form Fortran '!' comment, respecting string literals."""
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def main():
    # Optional argv[1] = repository root to scan. CI passes the PR-head tree here
    # so the TRUSTED base-branch copy of this script can lint untrusted PR code
    # (a PR must not be able to weaken the gate by editing this file). Default:
    # the tree this script lives in.
    global REPO, SRC, CMAKE, WRAPPER_FILES
    if len(sys.argv) > 1:
        REPO = os.path.abspath(sys.argv[1])
        SRC = os.path.join(REPO, "source")
        CMAKE = os.path.join(REPO, "cmake", "oqp_functions.cmake")
        WRAPPER_FILES = {
            os.path.normpath(os.path.join(SRC, "mathlib", "lapack_wrap.F90")),
            os.path.normpath(os.path.join(SRC, "mathlib", "blas_wrap.F90")),
        }

    reg = registered_syms()
    if not reg:
        print(f"ERROR: could not parse OQP_ACCELERATE_ILP64_SYMS from {CMAKE}",
              file=sys.stderr)
        return 2

    call_rx = re.compile(r"\bcall\s+([a-z][a-z0-9]+)\s*\(", re.I)
    func_rx = re.compile(r"(?<![\w%])([a-z][a-z0-9]+)\s*\(", re.I)
    defkw = re.compile(r"\b(subroutine|function|end|module|use|interface)\b", re.I)

    files = [p for p in scan_files() if p not in WRAPPER_FILES]
    defined = defined_procedures(files)

    violations = []
    for path in files:
        with open(path, encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
        for i, line in enumerate(lines, 1):
            code = strip_comment(line)
            # subroutine calls: flag known-BLAS and, via the shape heuristic,
            # a new BLAS/LAPACK external that CANON does not list. Use finditer,
            # not search: Fortran ';' statement separators allow more than one
            # `call` per line (a style the tree uses, e.g. hf_hessian.F90), and a
            # non-first call must be inspected too.
            for m in call_rx.finditer(code):
                name = m.group(1).lower()
                if name not in reg and (
                    name in CANON
                    or (BLAS_CALL_PAT.match(name)
                        and name not in defined
                        and name not in INTRINSIC_SUB)
                ):
                    violations.append(
                        (os.path.relpath(path, REPO), i, name, line.strip()))
            # BLAS/LAPACK *function* refs (ddot, dnrm2, dlamch, ...): limited to
            # CANON so name fragments never produce false positives.
            if not defkw.search(code):
                for fm in func_rx.finditer(code):
                    fn = fm.group(1).lower()
                    if fn in CANON and fn not in reg:
                        violations.append(
                            (os.path.relpath(path, REPO), i, fn, line.strip()))

    if violations:
        seen, uniq = set(), []
        for v in violations:
            k = (v[0], v[1], v[2])
            if k not in seen:
                seen.add(k)
                uniq.append(v)
        print("PR RULE 1 VIOLATION: BLAS/LAPACK routine(s) referenced from the "
              "Fortran compiled into liboqp but NOT registered in "
              "OQP_ACCELERATE_ILP64_SYMS.\n"
              "Such a routine has no oqp_*_i64 wrapper and no macOS Accelerate "
              "ILP64 alias, so the macOS build's check_accelerate_aliases gate will "
              "fail (OpenQP is ILP64-only, so this is the failure that bites).\n"
              "Fix: add the base name to OQP_ACCELERATE_ILP64_SYMS in "
              "cmake/oqp_functions.cmake, add its oqp_<name>_i64 wrapper in "
              "source/mathlib/{blas,lapack}_wrap.F90, and rename it in "
              "source/mathlib/oqp_linalg.F90 (see how dgels was added in ad43fae1).\n"
              "If a flagged name is an internal routine (not BLAS/LAPACK), it will "
              "have a `subroutine`/`function` definition in the scanned source and "
              "so should not be reported -- please open an issue if it is.\n")
        for rel, ln, name, txt in uniq:
            print(f"  {rel}:{ln}: {name}  ->  {txt}")
        print(f"\n{len(uniq)} unregistered BLAS/LAPACK reference(s).")
        return 1

    print(f"OK: every BLAS/LAPACK routine referenced in source/ and the "
          f"tests/fortran/*.F90 harnesses linked into liboqp is registered in "
          f"OQP_ACCELERATE_ILP64_SYMS ({len(reg)} symbols).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
