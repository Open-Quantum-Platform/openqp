"""A level shift must not reach published orbital energies on the pFON and SOSCF paths.

The DIIS loop publishes the orbitals of the previous Fock.  After a level
shift that Fock carried it, so the converged exit refreshes the orbital
energies from the final unshifted Fock.  The pFON exit used to leave the loop
before that refresh: RHF and ROHF published every virtual +vshift (H2O/6-31G
PBE, vshift=0.3: LUMO+n 1.210 instead of 0.910 Eh).  The decks are the
examples/SCF h2o *_pbe_pfon / *_pbe_soscf systems (H2O/6-31G PBE); with HF
instead of PBE the ROHF triplet lands on a different SCF solution when shifted
(-75.6256 vs -75.7192 Eh), which would test the state, not the published
energies.  ROHF under SOSCF also
published shifted energies before the zero-shift switch was fixed.

Each case runs the same deck with and without vshift and compares all alpha
(and beta) orbital energies; the shift is 0.3 Eh (0.15 Eh on ROHF open
shells), far above the 1e-5 Eh tolerance.
"""

import os
import re
import subprocess
import sys

import pytest

DECK = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
method=hf
functional=pbe

[guess]
type=huckel

[scf]
type={scf}
multiplicity={mult}
conv=1.0e-8
maxit=200
{extra}
"""

CASES = {
    "rhf-pfon": ("rhf", 1, "pfon=True"),
    "rohf-pfon": ("rohf", 3, "pfon=True"),
    "uhf-pfon": ("uhf", 3, "pfon=True"),
    "rhf-soscf": ("rhf", 1, "converger_type=soscf"),
    "rohf-soscf": ("rohf", 3, "converger_type=soscf"),
}


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(), reason="native OpenQP library not importable")


def _orbital_energies(log, spin):
    blocks = re.findall(spin + r" Orbitals -+\s*\n((?:.*\n)*?)\s*\n\s*=", log)
    assert blocks, log[-3000:]
    values = []
    for line in blocks[-1].splitlines():
        toks = line.split()
        if toks and all(re.match(r"-?\d+\.\d{6,}$", t) for t in toks):
            values += [float(t) for t in toks]
    return values


def _run(tmp_path, name, deck):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)],
                          capture_output=True, text=True, cwd=str(tmp_path),
                          env=env, timeout=900)
    log = inp.with_suffix(".log").read_text() if inp.with_suffix(".log").exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    return log


@pytest.mark.parametrize("case", sorted(CASES))
def test_level_shift_not_in_published_energies(tmp_path, case):
    scf, mult, extra = CASES[case]
    plain = _run(tmp_path, "plain", DECK.format(scf=scf, mult=mult, extra=extra))
    shifted = _run(tmp_path, "shifted",
                   DECK.format(scf=scf, mult=mult, extra=extra + "\nvshift=0.3"))
    energy = re.compile(r"Final \w+ energy is\s+(-?\d+\.\d+)")
    assert float(energy.findall(shifted)[-1]) == pytest.approx(
        float(energy.findall(plain)[-1]), abs=1e-7)
    for spin in (("Alpha", "Beta") if scf == "uhf" else ("Alpha",)):
        eps_plain = _orbital_energies(plain, spin)
        eps_shift = _orbital_energies(shifted, spin)
        assert len(eps_shift) == len(eps_plain)
        assert eps_shift == pytest.approx(eps_plain, abs=1e-5)
