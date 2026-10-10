"""GPU-free cross-validation of the MRSF/UMRSF METC contraction port.

The METC contraction is encoded **twice, independently**:

  * the CUDA kernels ``mrsf_metc_kernel`` / ``umrsf_metc_kernel`` in
    ``source/gpu_metc_cuda.cu`` (the device path), and
  * the Fortran host loops in ``int2_mrsf_data_t_update`` /
    ``int2_umrsf_data_t_update`` in ``source/tdhf_mrsf_lib.F90`` (the CPU
    reference / fallback path).

These must compute the *same* contraction, or a CPU/GPU parity run on real
hardware will silently disagree. A transcription error -- a swapped index
permutation, a wrong matrix-column (``m``) range, or a flipped Coulomb/exchange
sign -- is precisely the kind of bug that this test catches **without a GPU**.

Strategy: parse both source files, normalize every accumulation into a term
``(target_row, target_col, source_row, source_col, kind)`` keyed by Davidson
pass and matrix index ``m`` (0-based), and assert the Fortran and CUDA term
multisets are identical for every (pass, m). ``kind`` is ``C`` (Coulomb, ``+``)
or ``X`` (exchange, ``-``); the parser also asserts the literal sign matches the
kind, so a sign bug fails too.

This is a *computational-accuracy* check of the port itself, runnable in the
cloud (no CUDA toolchain, no built binary). It complements -- and must pass
before -- the on-hardware CPU/GPU parity validation (tools/gpu).
"""

import collections
import os
import re

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO_ROOT, "src", "metc")

IDX = "ijkl"  # the only loop indices that may appear as f3/d3 row/col


def _read(name):
    path = os.path.join(REPO_ROOT, "tests", "fixtures", name) if name.endswith(".F90") else os.path.join(SRC, name)
    with open(path) as fh:
        return fh.read()


# --------------------------------------------------------------------------
# CUDA kernel parsing.
# --------------------------------------------------------------------------
_CU_ADD4 = re.compile(
    r"add4\(f3,\s*f,\s*m,\s*([ijkl]),\s*([ijkl]),\s*nf,\s*nmatrix,\s*nbf,"
    r"\s*(-?)\s*(cval|xval)\s*\*\s*get4\(d3,\s*f,\s*m,\s*([ijkl]),\s*([ijkl])")
_CU_IF = re.compile(r"if\s*\((.*)\)\s*\{")
_CU_ELSEIF = re.compile(r"\}\s*else\s+if\s*\((.*)\)\s*\{")


def _parse_cond(cond):
    """Return (pass_set or None, m_set or None) for a C condition string."""
    passes = set(int(x) for x in re.findall(r"cur_pass\s*==\s*(\d+)", cond)) or None
    m_set = None
    lt = re.search(r"m\s*<\s*(\d+)", cond)
    eqs = re.findall(r"m\s*==\s*(\d+)", cond)
    if lt:
        m_set = set(range(int(lt.group(1))))
    elif eqs:
        m_set = set(int(x) for x in eqs)
    return passes, m_set


def _kind_and_sign(neg, coef):
    """Map (leading-minus, coef name) -> kind, asserting sign matches kind."""
    sign = "-" if neg == "-" else "+"
    if coef == "cval":
        assert sign == "+", "Coulomb term must be additive"
        return "C"
    assert coef == "xval"
    assert sign == "-", "exchange term must be subtractive"
    return "X"


def _slice(text, start_marker, end_marker):
    a = text.index(start_marker)
    b = text.index(end_marker, a + len(start_marker))
    return text[a:b]


def parse_cuda_kernel(text, kernel_name, end_marker):
    """Parse one CUDA kernel into {pass: {m: Counter(term)}}."""
    body = _slice(text, "void " + kernel_name, end_marker)
    result = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    stack = []  # list of (pass_set or None, m_set or None)

    def eff(which):
        cur = None
        for frame in stack:
            v = frame[which]
            if v is None:
                continue
            cur = set(v) if cur is None else (cur & v)
        return cur

    for raw in body.splitlines():
        line = raw.strip()
        if not line:
            continue
        m_elseif = _CU_ELSEIF.search(line)
        if m_elseif:
            if stack:
                stack.pop()
            stack.append(_parse_cond(m_elseif.group(1)))
            continue
        if line.startswith("}"):
            if stack:
                stack.pop()
            continue
        m_if = _CU_IF.match(line)
        if m_if:
            stack.append(_parse_cond(m_if.group(1)))
            continue
        if line.endswith("{"):  # function signature / other opener
            stack.append((None, None))
            continue
        for tr, tc, neg, coef, sr, sc in _CU_ADD4.findall(line):
            kind = _kind_and_sign(neg, coef)
            term = (tr, tc, sr, sc, kind)
            passes = eff(0) or {1, 2}
            ms = eff(1)
            assert ms is not None, "every add4 must sit under an m-condition: " + line
            for p in passes:
                for mm in ms:
                    result[p][mm][term] += 1
    return result


# --------------------------------------------------------------------------
# Fortran update parsing.
# --------------------------------------------------------------------------
_F_TERM = re.compile(
    r"f3\(\s*[^,]+,\s*([^,]+),\s*([ijkl]),\s*([ijkl])\)\s*=\s*"
    r"f3\([^)]*\)\s*([+\-])\s*(cval|xval)\s*\*\s*"
    r"d3\(\s*[^,]+,\s*([^,]+),\s*([ijkl]),\s*([ijkl])\)")


def _f_mrange(tok):
    """Expand a Fortran matrix-column range token into a 0-based m set."""
    tok = tok.strip()
    if tok.startswith(":"):          # ":4"  -> columns 1..4 -> {0,1,2,3}
        return set(range(0, int(tok[1:])))
    if ":" in tok:                    # "1:8" -> {0..7}; "9:10" -> {8,9}
        a, b = tok.split(":")
        return set(range(int(a) - 1, int(b)))
    return {int(tok) - 1}             # "11" -> {10}


def parse_fortran_update(text, sub_name):
    """Parse one Fortran update subroutine into {pass: {m: Counter(term)}}."""
    body = _slice(text, "subroutine " + sub_name, "end subroutine")
    result = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    passes = {1, 2}  # default outside an explicit cur_pass block
    for raw in body.splitlines():
        line = raw.strip()
        if line.startswith("!"):
            continue
        if "cur_pass==1" in line.replace(" ", ""):
            passes = {1}
        elif "cur_pass==2" in line.replace(" ", ""):
            passes = {2}
        elif re.match(r"end\s*if", line, re.I):
            passes = {1, 2}
        m = _F_TERM.search(line)
        if not m:
            continue
        mtok, tr, tc, sign, coef, src_mtok, sr, sc = m.groups()
        neg = "-" if sign == "-" else ""
        kind = _kind_and_sign(neg, coef)
        # LHS and RHS matrix ranges must reference the same columns.
        assert _f_mrange(mtok) == _f_mrange(src_mtok), \
            "f3/d3 column ranges differ on line: " + line
        term = (tr, tc, sr, sc, kind)
        for p in passes:
            for mm in _f_mrange(mtok):
                result[p][mm][term] += 1
    return result


# --------------------------------------------------------------------------
# Comparison helpers + tests.
# --------------------------------------------------------------------------
def _compare(fortran, cuda, label):
    passes = set(fortran) | set(cuda)
    assert passes, "%s: no terms parsed at all" % label
    for p in sorted(passes):
        fm, cm = fortran.get(p, {}), cuda.get(p, {})
        ms = set(fm) | set(cm)
        for mm in sorted(ms):
            fc, cc = fm.get(mm, collections.Counter()), cm.get(mm, collections.Counter())
            assert fc == cc, (
                "%s mismatch at pass=%d m=%d:\n  only in Fortran: %s\n  only in CUDA:    %s"
                % (label, p, mm, fc - cc, cc - fc))


def test_mrsf_cuda_matches_fortran():
    cu = parse_cuda_kernel(_read("gpu_metc_cuda.cu"),
                           "mrsf_metc_kernel", "void umrsf_metc_kernel")
    fo = parse_fortran_update(_read("tdhf_mrsf_lib.F90"), "int2_mrsf_data_t_update")
    _compare(fo, cu, "MRSF")


def test_umrsf_cuda_matches_fortran():
    cu = parse_cuda_kernel(_read("gpu_metc_cuda.cu"),
                           "umrsf_metc_kernel", "void umrsf_metc_combined_coulomb_kernel")
    fo = parse_fortran_update(_read("tdhf_mrsf_lib.F90"), "int2_umrsf_data_t_update")
    _compare(fo, cu, "UMRSF")


def test_mrsf_term_counts_are_sane():
    """Guard against a silently-empty parse (e.g. a regex drift)."""
    cu = parse_cuda_kernel(_read("gpu_metc_cuda.cu"),
                           "mrsf_metc_kernel", "void umrsf_metc_kernel")
    # pass 1: m in 0..3 carry 8 Coulomb + 8 exchange; m in 4..6 carry 8 exchange.
    assert sum(cu[1][m].total() for m in cu[1]) == 4 * 16 + 3 * 8
    # pass 2: only m == 6, 8 exchange terms.
    assert cu[2][6].total() == 8
    assert all(t[4] == "X" for t in cu[2][6])


def test_umrsf_term_counts_are_sane():
    cu = parse_cuda_kernel(_read("gpu_metc_cuda.cu"),
                           "umrsf_metc_kernel", "void umrsf_metc_combined_coulomb_kernel")
    # pass 1: m 0..7 -> 16 each; m 8,9 -> 8 each; m 10 -> 8.
    assert sum(cu[1][m].total() for m in cu[1]) == 8 * 16 + 2 * 8 + 8
    # pass 2: only m == 10, 8 exchange terms.
    assert cu[2][10].total() == 8


def test_parser_detects_injected_divergence():
    """The parser must actually be able to see a difference -- mutate one CUDA
    term and confirm the comparison fails (so a real divergence can't slip by)."""
    cu = parse_cuda_kernel(_read("gpu_metc_cuda.cu"),
                           "mrsf_metc_kernel", "void umrsf_metc_kernel")
    fo = parse_fortran_update(_read("tdhf_mrsf_lib.F90"), "int2_mrsf_data_t_update")
    # Corrupt a single CUDA term at pass 1, m 0.
    bad_term = next(iter(cu[1][0]))
    cu[1][0][bad_term] += 1
    raised = False
    try:
        _compare(fo, cu, "MRSF-injected")
    except AssertionError:
        raised = True
    assert raised, "comparison failed to detect an injected term divergence"


def load_tests(loader, tests, pattern):
    import unittest
    return unittest.TestSuite(unittest.FunctionTestCase(f) for name, f in globals().items()
                              if name.startswith('test_') and callable(f))
