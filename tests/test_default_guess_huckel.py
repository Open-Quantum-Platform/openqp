"""The default initial guess must stay Hückel.

Large, hard cases (e.g. the 1BNA benchmark: a charge -22 anion, 758 atoms,
RHF/STO-3G) converge smoothly from the Hückel guess in ~29 iterations but
oscillate from SAP, needing a DIIS reset to recover.  Hückel is also the
guess every pre-2025 OpenQP release defaulted to, so keeping it preserves
legacy reproducibility.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"


def test_default_guess_is_huckel():
    source = OQPDATA.read_text()
    assert "'type': {'type': string, 'default': 'huckel'}" in source
    assert "'type': {'type': string, 'default': 'sap'}" not in source
