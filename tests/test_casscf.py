"""End-to-end test for the clean CASSCF module (``method=casscf`` -> oqp.library.casscf).

CASSCF optimizes the active-space CI and the orbitals together (generalized-Fock
gradient + Newton).  The clean ``casscf.py`` replaced the retired CASSCF
development scaffolding; this validates it against PySCF CASSCF.

    H2O/STO-3G CAS(4,4) ncore=3 : OpenQP -75.0085688841 vs PySCF -75.0085688625 (2.2e-8)
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401  probe the dual-path oqp.utils shadow
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


def _run_casscf(tmp_path, project, *, system, basis, cas):
    from oqp.pyoqp import Runner

    runner = Runner(
        project=project,
        input_file=None,
        log=str(tmp_path / f"{project}.log"),
        input_dict={
            "input": {
                "system": system,
                "charge": "0",
                "basis": basis,
                "method": "casscf",
                "runtype": "energy",
            },
            "guess": {"type": "hcore"},
            "scf": {
                "type": "rhf",
                "multiplicity": "1",
                "maxit": "80",
                "forced_attempt": "3",
                "save_molden": "False",
            },
            "cas": cas,
            "ci": {
                "nroot": "1",
                "solver": "dense",
                "eig_tol": "1.0e-10",
                "integral_backend": "native",
                "target_spin": "any",
                "root_tracking": "energy",
            },
            "casscf": {
                "max_macro_iterations": "50",
                "optimizer": "newton",
                "root": "0",
                "gradient_norm_tol": "1.0e-7",
                "energy_decrease_tol": "1.0e-10",
                "step_norm_tol": "1.0e-8",
                "max_rotation_norm": "2.0e-1",
                "level_shift": "1.0e-3",
                "canonicalize": "true",
            },
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )
    runner.run(test_mod=True)
    return runner


def _casscf_energy(runner):
    return float(np.asarray(runner.mol.data["OQP::CASSCF_ENERGIES"], dtype=float)[0])


_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"


def test_casscf_h2o_cas44_matches_pyscf(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASSCF tests")

    runner = _run_casscf(
        tmp_path,
        "casscf_h2o_cas44",
        system=_H2O,
        basis="sto-3g",
        cas={"active_electrons": "4", "active_orbitals": "4", "frozen_core": "3", "max_det": "5000"},
    )
    e_casscf = _casscf_energy(runner)

    # PySCF CASSCF for the same active space; orbital optimization agrees to ~1e-8.
    assert e_casscf == pytest.approx(-75.0085688625, abs=1.0e-7)
    # CASSCF lowers the energy below the RHF reference (variational orbital opt).
    assert e_casscf < -74.963
    assert runner.mol.energies[0] == pytest.approx(e_casscf, abs=1.0e-12)
