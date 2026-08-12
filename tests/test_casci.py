"""End-to-end tests for the clean CASCI module (``method=casci`` -> oqp.library.casci).

CASCI is the active-space full CI at fixed RHF orbitals (a special case of the
verified determinant FCI).  The clean ``casci.py`` was extracted from the retired
19k-line development scaffolding; these tests validate it against PySCF CASCI /
OpenMolcas CIonly references.

    H4/STO-3G  CAS(2,2) frozen_core=1 : -2.1147334880 (PySCF)  / -2.1147334918 (Molcas)
    H2O/STO-3G CAS(2,2) frozen_core=4 : -74.9642716960 (PySCF/Molcas, benchmark)
    H2/STO-3G  CAS(2,2) frozen_core=0 : -1.1372838344 (FCI in the full STO-3G space)

The frozen_core=0 case is also a regression pin: with no core orbitals the
published ``*_TRIAL_CORE_ORBITALS`` tensor is empty, and tagarray aborts the
process on a zero-byte record, so the run used to die with SIGABRT *after*
printing correct energies.  Empty tensors are now skipped instead of stored.
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


def _run_casci(tmp_path, project, *, system, basis, cas):
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
                "method": "casci",
                "runtype": "energy",
            },
            "guess": {"type": "hcore"},
            "scf": {
                "type": "rhf",
                "multiplicity": "1",
                "maxit": "60",
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
            },
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )
    runner.run(test_mod=True)
    return runner


def _casci_energy(runner):
    return float(np.asarray(runner.mol.data["OQP::CASCI_ENERGIES"], dtype=float)[0])


_H4 = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"
_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"
_H2 = "\nH 0 0 0\nH 0 0 0.740"

# (name, system, basis, [cas], PySCF/Molcas CASCI energy)
_CASCI_CASES = [
    ("h4", _H4, "sto-3g",
     {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "5000"},
     -2.1147334880),
    ("h2o", _H2O, "sto-3g",
     {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "4", "max_det": "5000"},
     -74.9642716960),
    ("h2_no_frozen_core", _H2, "sto-3g",
     {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "0", "max_det": "5000"},
     -1.1372838344),
]


@pytest.mark.parametrize("name,system,basis,cas,ref", _CASCI_CASES, ids=[c[0] for c in _CASCI_CASES])
def test_casci_matches_reference(tmp_path, name, system, basis, cas, ref):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASCI tests")

    runner = _run_casci(tmp_path, f"casci_{name}", system=system, basis=basis, cas=cas)
    e_casci = _casci_energy(runner)

    # CASCI = active-space FCI; matches the independent reference to ~1e-8.
    assert e_casci == pytest.approx(ref, abs=1.0e-7)
    # the total energy stored on the molecule equals the CASCI ground-state energy
    assert runner.mol.energies[0] == pytest.approx(e_casci, abs=1.0e-12)


def test_casci_without_frozen_core_publishes_no_empty_core_tensor(tmp_path):
    """A zero-length tagarray record aborts the process, so it must be skipped.

    Reaching the assertions at all is most of the point: before the fix this
    run died with SIGABRT inside tagarray after the energies were printed.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASCI tests")

    runner = _run_casci(
        tmp_path,
        "casci_empty_core",
        system=_H2,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2",
             "frozen_core": "0", "max_det": "5000"},
    )

    assert _casci_energy(runner) == pytest.approx(-1.1372838344, abs=1.0e-7)

    # No core orbitals exist, so the tensor is absent rather than zero-length.
    with pytest.raises(Exception):
        runner.mol.data["OQP::CASCI_TRIAL_CORE_ORBITALS"]

    # The non-empty companion tensor is still published.
    active = np.asarray(runner.mol.data["OQP::CASCI_TRIAL_ACTIVE_ORBITALS"], dtype=int)
    assert active.size == 2


def test_casci_memory_guard_combines_resident_integrals_and_ci_workspace():
    from oqp.library.casci import CASCI
    from oqp.library.fci import ActiveSpacePlan, FCISettings

    driver = CASCI.__new__(CASCI)
    driver.settings = FCISettings(
        active_electrons=2,
        active_orbitals=30,
        max_det=1000,
        max_memory=110,
        nroot=1,
        solver="davidson",
    )
    plan = ActiveSpacePlan(
        norb=30,
        active=tuple(range(30)),
        core=(),
        nelec=(1, 1),
        metadata={},
    )
    with pytest.raises(ValueError, match="combined integral and CI working set"):
        driver._check_combined_ci_memory(30, plan)
