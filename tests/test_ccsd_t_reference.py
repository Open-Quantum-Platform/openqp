"""End-to-end CCSD(T) validation against PySCF reference energies.

Runs the whole native path -- OpenQP's own integral engine, the packed AO
store, the cc_ao2mo transformation and the cc_lib solver -- and compares the
correlation energies with values computed by PySCF on the same geometry and
basis (see tests/data/ccsd_t_pyscf_validation.json).

Skipped unless a built OpenQP is importable, so it is inert in source-only
checkouts and becomes a real regression gate wherever the library exists.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "ccsd_t_pyscf_validation.json"
DATA_OS = Path(__file__).parent / "data" / "ccsd_t_open_shell_validation.json"


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    import oqp
    return hasattr(oqp, "ccsd_t_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library with ccsd_t_energy not importable",
)


def _load():
    doc = json.loads(DATA.read_text())
    return doc["tolerance"], doc["cases"]


def _load_os():
    doc = json.loads(DATA_OS.read_text())
    return doc["tolerance"], doc["cases"]


def _write_input(path, case, triples=True):
    atoms = "\n".join(
        "   %d   %.10f   %.10f   %.10f" % (a[0], a[1], a[2], a[3])
        for a in case["atoms"]
    )
    method = "ccsd(t)" if triples else "ccsd"
    mult = case.get("spin", 0) + 1
    scf_type = case.get("scf_type", "rhf")
    path.write_text(
        "[input]\nsystem=\n%s\ncharge=0\nruntype=energy\nbasis=%s\nmethod=%s\n\n"
        "[guess]\ntype=huckel\n\n"
        "[scf]\nmultiplicity=%d\ntype=%s\nconv=1.0e-11\n\n"
        "[cc]\nnfzc=%d\nconv=1e-9\n"
        % (atoms, case["basis"], method, mult, scf_type, case["nfzc"])
    )


def _run(tmp_path, case, triples=True):
    inp = tmp_path / "case.inp"
    _write_input(inp, case, triples)
    log = inp.with_suffix(".log")
    if log.exists():
        log.unlink()
    subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    assert log.exists(), "no log produced -- the run did not reach the CC module"
    text = log.read_text()

    def grab(pattern):
        m = re.search(pattern, text)
        return float(m.group(1)) if m else None

    return {
        "nbf": grab(r"nbf = (\d+)"),
        "scf": grab(r"Final \w+ energy is\s+(-?[\d.]+)"),
        "ccsd": grab(r"E\(CCSD, correlation\)\s*=\s*(-?[\d.]+)"),
        "triples": grab(r"E\(\(T\), correction\)\s*=\s*(-?[\d.]+)"),
    }


@pytest.mark.parametrize("case", _load()[1], ids=lambda c: c["name"])
def test_ccsd_t_matches_pyscf(tmp_path, case):
    tol, _ = _load()
    got = _run(tmp_path, case)

    assert got["nbf"] == case["nbf"], (
        "basis size differs from the reference (%s vs %s). For Pople sets this "
        "usually means a Cartesian/spherical mismatch, in which case the two "
        "codes are not running the same basis and the energies below are not "
        "comparable." % (got["nbf"], case["nbf"])
    )

    # A correlation energy is only meaningful against its reference determinant,
    # so the SCF has to line up before the CC numbers are worth comparing.
    assert abs(got["scf"] - case["scf_reference"]) < tol["scf"], (
        "SCF reference differs; CC comparison would be meaningless"
    )

    assert abs(got["ccsd"] - case["e_ccsd"]) < tol["ccsd"]
    assert abs(got["triples"] - case["e_triples"]) < tol["triples"]


def test_ccsd_without_triples_matches_ccsd_part(tmp_path):
    """method=ccsd must reproduce the CCSD part and not report a (T) term."""
    tol, cases = _load()
    case = next(c for c in cases if c["name"] == "H2O/cc-pVDZ")
    got = _run(tmp_path, case, triples=False)
    assert abs(got["ccsd"] - case["e_ccsd"]) < tol["ccsd"]
    assert got["triples"] is None


@pytest.mark.parametrize("case", _load_os()[1], ids=lambda c: c["name"])
def test_open_shell_ccsd_t_matches_pyscf(tmp_path, case):
    """UHF references go through the spin-orbital solver in cc_uhf_lib."""
    tol, _ = _load_os()
    got = _run(tmp_path, case)

    assert abs(got["scf"] - case["scf_reference"]) < tol["scf"], (
        "SCF reference differs; CC comparison would be meaningless"
    )
    assert abs(got["ccsd"] - case["e_ccsd"]) < tol["ccsd"]
    assert abs(got["triples"] - case["e_triples"]) < case.get(
        "triples_tolerance", tol["triples"]
    )
