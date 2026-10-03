#!/usr/bin/env python3
"""
Enforce PR RULE 8 (portable, source-level): every derivative of the electron-nucleus
attraction must see the ECP-screened nuclear charges, and every one-electron gradient
or coupling assembly that contracts a density with dV_en/dx must also contract it with
the ECP derivative.

WHY: with an effective core potential the nucleus carries Z - N_core electrons' worth of
charge for the valence basis (`infos%basis%ecp_zn_num` holds N_core per atom), and the
potential acquires the semilocal ECP term whose derivative is `grad_1e_ecp` (packed
density) or `ecp_deriv_ints` (AO derivative matrices). hf_gradient, tdhf_gradient and
the CASSCF/PT2 gradients got both when ECPs arrived (#46); sf_1e_grad (SF/MRSF
gradient) kept the full Z for eight months, and the analytic MRSF NAC 1e term had the
full Z and no ECP derivative. Neither was visible to any test: the only ECP MRSF
gradient example removed no core electrons (LANL2DZ on C/H), the heavy-atom ECP examples
were energy-only, and everything else was all-electron. HBr/LANL2DZ exposed 6.3 Ha/Bohr
(gradient) and a factor 140 (NAC).

WHAT IS CHECKED (per call site, per enclosing subroutine):
  1. The charge argument of grad_en_hellman_feynman / grad_en_pulay / der_nucattr_matrix
     either contains `ecp_zn_num` itself, or is a name whose defining `=>`/`=` line inside
     the same subroutine contains `ecp_zn_num`.
  2. A subroutine that calls grad_en_pulay or grad_en_hellman_feynman also references
     grad_1e_ecp or ecp_deriv_ints (the ECP derivative of the same density).

This is a static gate: it proves the question was asked, not that the answer is right.
Run an ECP example with a heavy atom whose core is actually removed (HBr/LANL2DZ,
NaCl/SBKJC) against finite differences or the numerical NAC before trusting a new term.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHARGE_CALLS = ("grad_en_hellman_feynman", "grad_en_pulay", "der_nucattr_matrix")
ECP_DERIV = ("grad_1e_ecp", "ecp_deriv_ints")

SUB_RE = re.compile(r"^\s*(?:pure\s+|elemental\s+|recursive\s+)*(?:subroutine|function)\s+(\w+)", re.I)
END_RE = re.compile(r"^\s*end\s+(?:subroutine|function)\b", re.I)


def join_continuations(lines):
    """Yield (first_line_number, logical_line) with Fortran '&' continuations joined,
    comments stripped."""
    buf, start = "", None
    for no, raw in enumerate(lines, 1):
        code = raw.split("!", 1)[0].rstrip()
        if start is None:
            start = no
        stripped = code.lstrip()
        if stripped.startswith("&"):
            stripped = stripped[1:]
        if code.rstrip().endswith("&"):
            buf += stripped.rstrip()[:-1]
            continue
        buf += stripped
        yield start, buf
        buf, start = "", None


def split_args(text):
    depth, cur, out = 0, "", []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            if depth == 0:
                break
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur.strip()); cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def check_file(path):
    problems = []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    logical = list(join_continuations(lines))
    # split into subroutine scopes
    scopes, cur, cur_name, cur_start = [], [], None, 0
    for no, text in logical:
        m = SUB_RE.match(text)
        if m and cur_name is None:
            cur_name, cur_start, cur = m.group(1), no, []
        cur.append((no, text))
        if END_RE.match(text) and cur_name is not None:
            scopes.append((cur_name, cur_start, cur)); cur_name, cur = None, []
    for name, start, body in scopes:
        texts = [t for _, t in body]
        charge_calls = []
        for no, text in body:
            for call in CHARGE_CALLS:
                m = re.search(r"\bcall\s+%s\s*\(" % call, text, re.I)
                if not m:
                    continue
                args = split_args(text[m.end():])
                if len(args) < 3:
                    problems.append(f"{path}:{no}: {call}: cannot parse arguments"); continue
                zarg = args[2]
                charge_calls.append((no, call, zarg))
        if not charge_calls:
            continue
        for no, call, zarg in charge_calls:
            if "ecp_zn_num" in zarg:
                continue
            ident = re.match(r"[A-Za-z_]\w*", zarg.strip())
            ok = False
            if ident:
                defn = re.compile(r"(?:^|[,(\s])%s\s*(=>|=)\s*(.*)" % re.escape(ident.group(0)), re.I)
                for t in texts:
                    dm = defn.search(t)
                    if dm and "ecp_zn_num" in dm.group(2):
                        ok = True; break
            if not ok:
                problems.append(f"{path}:{no}: {name}: {call} receives '{zarg}' which is not "
                                f"ECP-screened (expected '... - ecp_zn_num' in the argument or in "
                                f"its defining line inside {name})")
        needs_deriv = any(c in ("grad_en_pulay", "grad_en_hellman_feynman") for _, c, _ in charge_calls)
        if needs_deriv and not any(re.search(r"\b(%s)\b" % "|".join(ECP_DERIV), t, re.I) for t in texts):
            problems.append(f"{path}:{start}: {name}: contracts a density with dV_en/dx but never "
                            f"with the ECP derivative (grad_1e_ecp / ecp_deriv_ints)")
    return problems


def main():
    src = ROOT / "source"
    problems = []
    for path in sorted(src.rglob("*.F90")) + sorted(src.rglob("*.f90")):
        if path.name == "grd1.F90":      # the integral routines themselves
            continue
        problems.extend(check_file(path))
    if problems:
        print("ECP charge/derivative gate (rule 8): FAIL")
        for p in problems:
            print("  " + p)
        return 1
    print("ECP charge/derivative gate (rule 8): PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
