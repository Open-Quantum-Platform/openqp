"""IR intensities are reported in km/mol and Raman activities in bohr^4/amu.

The native kernel (source/modules/vibrational_intensities.F90) multiplies
|d mu / d Q|^2 by a conversion factor.  The dipole derivatives it receives are
in atomic units (e*bohr/bohr = e) and the normal modes are Cartesian
displacements per unit sqrt(amu)*bohr, so |d mu / d Q|^2 is in e^2/amu and the
factor is N_A / (12 eps0 c^2) * e^2/amu = 974.880 km/mol.  Until 2026-10-20
both the kernel and pyoqp/oqp/library/frequency.py used 42.255, the factor for
(Debye/Angstrom)^2/amu, so every "IR(km/mol)" value in the log, hess.json and
freq.molden was 23.07 times too small.

Raman activities 45 abar'^2 + 7 gamma'^2 are formed from polarizability
derivatives in bohr^3 per sqrt(amu)*bohr, so they are in bohr^4/amu (multiply
by 0.078415 for Angstrom^4/amu); the label was 'a.u.' until 2026-10-20.

The end-to-end test compares water RHF/6-31G with an independent PySCF 2.14
calculation at the same (PySCF-optimized) geometry and with the same isotope
masses: analytic PySCF Hessian; dipole derivatives by central differences
(1e-3 bohr) of the SCF dipole; polarizability by central finite fields
(5e-4 au); IR from scipy.constants with no literal conversion factor.  Script:
pyscf_ir_raman_ref.py in the 2026-10-20 ir-intensity-units session folder.
The old factor fails the IR comparison by a factor of 23.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest


# PySCF-optimized RHF/6-31G water (max |gradient| 5e-9 Eh/bohr), Angstrom.
H2O_XYZ = """3
water RHF/6-31G
O  0.0000000000  0.0000000000  0.0823973117
H  0.0000000000  0.7851665083 -0.4517486558
H  0.0000000000 -0.7851665083 -0.4517486558
"""

# PySCF 2.14, masses O 15.99491, H 1.007825 (the OpenQP masses).
PYSCF_FREQS = np.array([1737.009, 3988.502, 4145.433])          # cm-1
PYSCF_IR = np.array([123.05524, 2.95518, 54.30713])             # km/mol
PYSCF_RAMAN = np.array([135.53278, 1149.47891, 511.23172])      # bohr^4/amu


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import main  # noqa: F401
        return True
    except Exception:
        return False


def test_python_conversion_factor_is_km_mol_per_e2_amu():
    from oqp.library.frequency import IR_INTENSITY_CONVERSION_KM_MOL

    # Independent SI route through eps0 (the definition uses alpha and h).
    n_a, e, amu = 6.02214076e23, 1.602176634e-19, 1.66053906660e-27
    eps0, c = 8.8541878128e-12, 299792458.0
    expected = n_a * e ** 2 / amu / (12.0 * eps0 * c ** 2) / 1000.0
    assert IR_INTENSITY_CONVERSION_KM_MOL == pytest.approx(expected, rel=1e-8)
    assert IR_INTENSITY_CONVERSION_KM_MOL == pytest.approx(974.880, abs=1e-3)


@pytest.mark.skipif(not _runtime_available(), reason="compiled OpenQP runtime not available")
def test_water_rhf_ir_raman_match_pyscf(tmp_path):
    from oqp.library.frequency import IR_INTENSITY_CONVERSION_KM_MOL

    (tmp_path / "h2o.xyz").write_text(H2O_XYZ, encoding="utf-8")
    deck = tmp_path / "h2o_ir.oqp"
    deck.write_text('rhf/6-31g hess(S0,type=analytical) ir raman scf(conv=1e-10)\n'
                    'geom="h2o.xyz"\n', encoding="utf-8")
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "2")
    proc = subprocess.run(
        [sys.executable, "-c", "import sys; from oqp.pyoqp import main; sys.exit(main())",
         deck.name],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, (proc.stdout[-3000:], proc.stderr[-3000:])

    data = json.loads((tmp_path / "h2o_ir.hess.json").read_text(encoding="utf-8"))
    assert data["vibrational_intensity_metadata"]["status"] == "computed"
    assert data["vibrational_intensity_metadata"]["backend"] == "native_openqp_analytic"
    assert data["vibrational_intensity_metadata"]["ir_units"] == "km/mol"
    assert data["vibrational_intensity_metadata"]["raman_units"] == "bohr^4/amu"
    assert "IR(km/mol)      Raman(bohr^4/amu)" in (tmp_path / "h2o_ir.log").read_text(
        encoding="utf-8", errors="replace")
    freqs = np.asarray(data["freqs"], dtype=float)
    order = np.argsort(freqs)
    freqs = freqs[order]
    ir = np.asarray(data["infrared_intensities"], dtype=float)[order]
    raman = np.asarray(data["raman_activities"], dtype=float)[order]
    mode_dipoles = np.asarray(data["infrared_mode_dipole_derivatives"], dtype=float)[order]

    # Measured |diff| 0.011 cm-1 (same Hessian theory, different integral codes).
    np.testing.assert_allclose(freqs, PYSCF_FREQS, atol=0.1)
    # The native kernel and the Python helper use the same factor.
    np.testing.assert_allclose(
        ir, IR_INTENSITY_CONVERSION_KM_MOL * np.sum(mode_dipoles ** 2, axis=1), rtol=1e-10)
    # External reference; the old 42.255 factor is off by 23x.
    np.testing.assert_allclose(ir, PYSCF_IR, rtol=0.01)
    np.testing.assert_allclose(raman, PYSCF_RAMAN, rtol=0.01)
