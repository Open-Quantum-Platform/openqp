"""ROHF MP2 must include the second-order singles term.

Semicanonicalisation diagonalises the occupied-occupied and virtual-virtual
Fock blocks but leaves the occupied-virtual block alone, and for ROHF that
block does not vanish: the ROHF stationarity condition fixes one combination
of the two spins, not each spin separately. The resulting

    E_S = sum_(ia,sigma) |f^sigma_ia|^2 / (e^sigma_i - e^sigma_a)

is a genuine part of E(2) on such a reference. Omitting it left the ROHF MP2
energy short by 1-12 mHa depending on the system, silently.

The values below were checked against an independent evaluation of the same
expression built from PySCF's semicanonical Fock matrices, agreeing to ~1e-10.

Skipped unless a built OpenQP is importable.
"""

import re
import subprocess
import sys

import pytest


def _have_native_oqp():
    try:
        import oqp
    except Exception:
        return False
    return hasattr(oqp, "mp2_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(), reason="native OpenQP library not importable"
)

TEMPLATE = """[input]
system=
{atoms}
charge=0
method=mp2
basis={basis}
runtype=energy
functional=

[guess]
type=huckel

[scf]
type={scf}
multiplicity={mult}
conv=1.0e-11
"""

OH = " 8 0.0 0.0 0.0\n 1 0.0 0.0 0.97"
NH2 = " 7 0.0 0.0 0.0\n 1 0.0 0.8 0.57\n 1 0.0 -0.8 0.57"
O2 = " 8 0.0 0.0 0.0\n 8 0.0 0.0 1.21"
H2O = (" 8 0.0 0.0 -0.041061554\n"
       " 1 -0.533194329 0.533194329 -0.614469223\n"
       " 1 0.533194329 -0.533194329 -0.614469223")


def _run(tmp_path, atoms, basis, scf, mult, name="mp2"):
    inp = tmp_path / (name + ".inp")
    inp.write_text(TEMPLATE.format(atoms=atoms, basis=basis, scf=scf, mult=mult))
    log = inp.with_suffix(".log")
    if log.exists():
        log.unlink()
    subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)],
                   capture_output=True, cwd=str(tmp_path), timeout=1800)
    assert log.exists(), "no log produced"
    text = log.read_text(errors="replace")

    def grab(pattern):
        m = re.search(pattern, text)
        return float(m.group(1)) if m else None

    return {
        "singles": grab(r"E\(MP2, singles\)\s*=\s*(-?[\d.]+)"),
        "corr": grab(r"E\(MP2, correlation\)\s*=\s*(-?[\d.]+)"),
    }


@pytest.mark.parametrize("name,atoms,basis,mult,e_singles", [
    ("OH/cc-pVDZ", OH, "cc-pvdz", 2, -0.0027250078),
    ("NH2/6-31G", NH2, "6-31g", 2, -0.0010842111),
    ("O2/cc-pVDZ", O2, "cc-pvdz", 3, -0.0121300022),
])
def test_rohf_singles_matches_the_independent_value(tmp_path, name, atoms,
                                                    basis, mult, e_singles):
    got = _run(tmp_path, atoms, basis, "rohf", mult)
    assert got["singles"] is not None, "ROHF MP2 must report a singles term"
    assert abs(got["singles"] - e_singles) < 1.0e-7, (name, got)


def test_closed_shell_gives_the_same_answer_through_rhf_and_rohf(tmp_path):
    """The strongest internal check: for a closed shell ROHF is RHF, so the
    singles term must vanish and the correlation energy must be identical."""
    rhf = _run(tmp_path, H2O, "cc-pvdz", "rhf", 1, name="rhf")
    rohf = _run(tmp_path, H2O, "cc-pvdz", "rohf", 1, name="rohf")

    assert rhf["corr"] is not None and rohf["corr"] is not None
    assert abs(rhf["corr"] - rohf["corr"]) < 1.0e-9, (rhf, rohf)
    # Zero singles are not printed at all, so neither run should show the line.
    assert rhf["singles"] is None
    assert rohf["singles"] is None


def test_uhf_is_untouched(tmp_path):
    """A canonical UHF reference has f_ia = 0, so nothing may change there."""
    got = _run(tmp_path, OH, "cc-pvdz", "uhf", 2)
    assert got["singles"] is None, "canonical UHF must have no singles term"
    assert abs(got["corr"] - (-0.1510087705)) < 1.0e-7, got
