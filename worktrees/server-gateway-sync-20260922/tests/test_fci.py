import importlib
import os
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

# Allow the pure-Python FCI components (solve_fci and the input checker) to be
# imported and tested without a compiled liboqp build. Tests that genuinely need
# the native backend skip themselves when it is unavailable.
os.environ.setdefault("OQP_BACKEND_OPTIONAL", "1")
ROOT = Path(__file__).resolve().parents[1]
PYOQP = ROOT / "pyoqp"
if str(PYOQP) not in sys.path:
    sys.path.insert(0, str(PYOQP))


def _use_real_oqp_package():
    # Several lightweight unit tests stub oqp.* modules in sys.modules. The FCI
    # tests exercise the real installed package, so they must not inherit those
    # stubs when the whole suite runs in filename order.
    for name in list(sys.modules):
        if name == "oqp" or name.startswith("oqp."):
            sys.modules.pop(name, None)
    importlib.invalidate_caches()


def _backend_available() -> bool:
    """True when the compiled liboqp backend loaded (required for Runner tests)."""
    _use_real_oqp_package()
    import oqp

    return bool(getattr(oqp, "BACKEND_AVAILABLE", False))


def _pyscf_reference(atom, charge=0, spin=0, basis="sto-3g"):
    pyscf = pytest.importorskip("pyscf")
    from pyscf import ao2mo, fci, gto, scf

    mol = gto.Mole()
    mol.atom = atom
    mol.unit = "Angstrom"
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build(cart=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    assert mf.converged

    coeff = mf.mo_coeff
    h1e = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, coeff), coeff.shape[1])
    energy_ref, _ = fci.FCI(mol, coeff).kernel()
    return h1e, eri, mol.nelec, mol.energy_nuc(), energy_ref


def _pyscf_casci_reference(
    atom, active_electrons, active_orbitals, charge=0, spin=0, basis="sto-3g"
):
    pytest.importorskip("pyscf")
    from pyscf import gto, mcscf, scf

    mol = gto.Mole()
    mol.atom = atom
    mol.unit = "Angstrom"
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build(cart=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    assert mf.converged

    return mcscf.CASCI(mf, active_orbitals, active_electrons).kernel()[0]


def _pyscf_fci_roots(atom, nroots, charge=0, spin=0, basis="sto-3g"):
    pyscf = pytest.importorskip("pyscf")
    from pyscf import ao2mo, fci, gto, scf
    from pyscf.fci import spin_op

    mol = gto.Mole()
    mol.atom = atom
    mol.unit = "Angstrom"
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build(cart=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    assert mf.converged

    coeff = mf.mo_coeff
    h1e = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, coeff), coeff.shape[1])
    energies, ci_roots = fci.FCI(mol, coeff).kernel(nroots=nroots)
    s2_ref = [
        spin_op.spin_square(ci, coeff.shape[1], mol.nelec)[0]
        for ci in ci_roots
    ]
    return h1e, eri, mol.nelec, mol.energy_nuc(), np.asarray(energies), np.asarray(s2_ref)


def _assert_spin_spectrum_matches(energies, s2, ref_energies, ref_s2):
    unmatched = list(range(len(ref_energies)))
    for energy, spin in zip(energies, s2):
        candidates = [
            idx
            for idx in unmatched
            if abs(ref_energies[idx] - energy) < 1.0e-8
            and abs(ref_s2[idx] - spin) < 1.0e-8
        ]
        assert candidates, (
            f"no PySCF spin match for energy={energy:.12f}, s2={spin:.12f}; "
            f"remaining={[(ref_energies[i], ref_s2[i]) for i in unmatched]}"
        )
        unmatched.remove(candidates[0])


def test_compute_s2_for_analytic_high_spin_determinant():
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, compute_s2, fci_spin_diagnostics

    dets = _determinants(2, (2, 0))
    ci = np.array([1.0])

    assert compute_s2(ci, dets, 2, (2, 0)) == pytest.approx(2.0, abs=1.0e-12)

    s2, multiplicity = fci_spin_diagnostics(ci[:, None], dets, 2, (2, 0))
    assert s2[0] == pytest.approx(2.0, abs=1.0e-12)
    assert multiplicity[0] == 3


@pytest.mark.parametrize(
    ("ci_vector", "message"),
    [
        (np.asarray([True], dtype=bool), "CI vector"),
        (["1.0"], "CI vector"),
        (np.asarray([1.0 + 0.0j], dtype=complex), "CI vector"),
    ],
)
def test_compute_s2_rejects_nonreal_ci_vectors(ci_vector, message):
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, compute_s2

    dets = _determinants(2, (2, 0))

    with pytest.raises(ValueError, match=message):
        compute_s2(ci_vector, dets, 2, (2, 0))


@pytest.mark.parametrize(
    "ci_vectors",
    [
        np.asarray([[True]], dtype=bool),
        [["1.0"]],
        np.asarray([[1.0 + 0.0j]], dtype=complex),
    ],
)
def test_fci_spin_diagnostics_rejects_nonreal_ci_vectors(ci_vectors):
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, fci_spin_diagnostics

    dets = _determinants(2, (2, 0))

    with pytest.raises(ValueError, match="CI vectors"):
        fci_spin_diagnostics(ci_vectors, dets, 2, (2, 0))


@pytest.mark.parametrize(
    ("target_spin", "canonical", "multiplicity"),
    [
        ("doublet", "doublet", 2),
        ("Quartet", "quartet", 4),
        ("quintet", "quintet", 5),
        ("sextet", "sextet", 6),
        ("7", "7", 7),
    ],
)
def test_target_spin_named_multiplicities_are_canonicalized(
    target_spin,
    canonical,
    multiplicity,
):
    _use_real_oqp_package()
    from oqp.library.fci import _target_spin_multiplicity, _target_spin_setting

    normalized = _target_spin_setting(target_spin, "target_spin")

    assert normalized == canonical
    assert _target_spin_multiplicity(target_spin) == multiplicity


def test_target_spin_filter_accepts_named_higher_multiplicity():
    _use_real_oqp_package()
    from oqp.library.fci import _filter_roots_by_target_spin

    energies, coeffs, s2, multiplicity, root_indices = _filter_roots_by_target_spin(
        np.asarray([-1.0, -0.8, -0.6], dtype=float),
        np.eye(3, dtype=float),
        np.asarray([0.0, 2.0, 6.0], dtype=float),
        np.asarray([1, 3, 5], dtype=np.int64),
        target_spin="quintet",
        requested_nroot=1,
        ci_label="FCI",
        ci_section="[fci]",
    )

    np.testing.assert_allclose(energies, [-0.6])
    np.testing.assert_allclose(coeffs, np.eye(3, dtype=float)[:, [2]])
    np.testing.assert_allclose(s2, [6.0])
    assert multiplicity.tolist() == [5]
    assert root_indices.tolist() == [2]


def test_solve_fci_filters_named_target_spin_roots():
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, fci_spin_diagnostics, solve_fci

    h1e = np.diag([0.0, 1.0])
    eri = np.zeros((2, 2, 2, 2), dtype=float)
    eri[0, 1, 1, 0] = 0.5
    eri[1, 0, 0, 1] = 0.5

    energies, coeffs = solve_fci(
        h1e,
        eri,
        (1, 1),
        nroot=1,
        max_det=100,
        solver="dense",
        target_spin="triplet",
    )
    s2, multiplicity = fci_spin_diagnostics(
        coeffs,
        _determinants(2, (1, 1)),
        2,
        (1, 1),
    )

    np.testing.assert_allclose(energies, [0.5])
    assert coeffs.shape == (4, 1)
    np.testing.assert_allclose(s2, [2.0])
    assert multiplicity.tolist() == [3]


def test_solve_fci_davidson_target_spin_expands_past_subspace_override():
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, fci_spin_diagnostics, solve_fci

    h1e = np.diag([0.0, 1.0])
    eri = np.zeros((2, 2, 2, 2), dtype=float)
    eri[0, 1, 1, 0] = 0.5
    eri[1, 0, 0, 1] = 0.5

    energies, coeffs = solve_fci(
        h1e,
        eri,
        (1, 1),
        nroot=1,
        max_det=100,
        solver="davidson",
        davidson_subspace=1,
        target_spin="triplet",
    )
    s2, multiplicity = fci_spin_diagnostics(
        coeffs,
        _determinants(2, (1, 1)),
        2,
        (1, 1),
    )

    np.testing.assert_allclose(energies, [0.5])
    assert coeffs.shape == (4, 1)
    np.testing.assert_allclose(s2, [2.0])
    assert multiplicity.tolist() == [3]


def test_solve_fci_target_spin_reports_no_matching_roots():
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e = np.diag([0.0, 1.0])
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    with pytest.raises(ValueError, match="target_spin=quintet.*4 solved roots"):
        solve_fci(
            h1e,
            eri,
            (1, 1),
            nroot=1,
            max_det=100,
            solver="dense",
            target_spin="quintet",
        )


def test_solve_fci_target_spin_ignores_ambiguous_roots_above_the_selection():
    # The degenerate open-shell pair at E=1 comes back as spin-mixed
    # determinants labelled "doublet" (impossible for two electrons), but the
    # requested singlet ground state at E=0 is pure and unambiguous, so the
    # mislabelled roots above the selection window must not abort the solve.
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e = np.diag([0.0, 1.0])
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    with pytest.warns(RuntimeWarning, match="wrong parity"):
        energies, coeffs = solve_fci(
            h1e,
            eri,
            (1, 1),
            nroot=1,
            max_det=100,
            solver="dense",
            target_spin="singlet",
        )

    np.testing.assert_allclose(energies, [0.0], atol=1.0e-12)
    assert coeffs.shape == (4, 1)


def test_solve_fci_target_spin_still_refuses_ambiguity_inside_the_selection():
    # Asking for two singlets pushes the selection window past the mislabelled
    # degenerate pair at E=1; the labels then decide the answer and the solve
    # must refuse rather than silently return energies of spin mixtures.
    _use_real_oqp_package()
    from oqp.library.fci import SpinLabelAmbiguityError, solve_fci

    h1e = np.diag([0.0, 1.0])
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    with pytest.raises(SpinLabelAmbiguityError, match="target_spin=singlet"), \
            pytest.warns(RuntimeWarning, match="wrong parity"):
        solve_fci(
            h1e,
            eri,
            (1, 1),
            nroot=2,
            max_det=100,
            solver="dense",
            target_spin="singlet",
        )


def test_ci_vector_log_entries_apply_threshold_and_active_labels():
    _use_real_oqp_package()
    from oqp.library.fci import _ci_vector_log_entries, _determinants

    determinants = _determinants(2, (1, 1))
    coeffs = np.array(
        [
            [0.95, 0.02],
            [0.03, -0.60],
            [0.01, 0.04],
            [-0.20, 0.80],
        ],
        dtype=float,
    )

    payload = _ci_vector_log_entries(
        coeffs,
        determinants,
        2,
        threshold=0.05,
        active_orbital_indices="3,5",
        root_indices=[0, 2],
    )

    assert payload["threshold"] == pytest.approx(0.05, abs=0.0)
    assert payload["root_indices"] == (0, 2)
    assert [entry["index"] for entry in payload["entries"]] == [0, 1, 3]
    assert payload["entries"][0]["alpha"] == "3"
    assert payload["entries"][0]["beta"] == "3"
    assert payload["entries"][1]["alpha"] == "3"
    assert payload["entries"][1]["beta"] == "5"
    assert payload["entries"][2]["alpha"] == "5"
    assert payload["entries"][2]["beta"] == "5"


@pytest.mark.parametrize(
    ("threshold", "root_indices", "message"),
    [
        (True, [0, 1], "CI print threshold"),
        ("0.05", [0, 1], "CI print threshold"),
        (float("nan"), [0, 1], "CI print threshold"),
        (0.05, [True, False], "CI root indices"),
        (0.05, ["0", "1"], "CI root indices"),
        (0.05, [0.0, 1.0], "CI root indices"),
        (0.05, [0], "CI root indices length"),
    ],
)
def test_ci_vector_log_entries_rejects_non_native_diagnostic_inputs(
    threshold,
    root_indices,
    message,
):
    _use_real_oqp_package()
    from oqp.library.fci import _ci_vector_log_entries, _determinants

    determinants = _determinants(2, (1, 1))
    coeffs = np.zeros((4, 2), dtype=float)

    with pytest.raises(ValueError, match=message):
        _ci_vector_log_entries(
            coeffs,
            determinants,
            2,
            threshold=threshold,
            root_indices=root_indices,
        )


def test_dump_log_renders_fci_ci_vector_rows(tmp_path):
    _use_real_oqp_package()
    from oqp.utils.file_utils import dump_log

    class FakeData:
        _data = SimpleNamespace(
            mol_prop=SimpleNamespace(natom=2, charge=0, mult=1),
            control=SimpleNamespace(scftype=1, maxit=50),
            tddft=SimpleNamespace(mult=1),
        )

    log_file = tmp_path / "fci_ci_vectors.log"
    mol = SimpleNamespace(
        log=str(log_file),
        config={
            "input": {"method": "casci", "basis": "sto-3g", "functional": ""},
            "scf": {
                "forced_attempt": 1,
                "conv": 1.0e-8,
                "incremental": False,
                "diis_type": "cdiis",
                "cdiis_switch": 0.3,
                "vdiis_cdiis_switch": 0.0,
                "vdiis_vshift_switch": 0.0,
                "vshift_cdiis_switch": 0.0,
                "vshift": 0.0,
            },
            "tdhf": {
                "type": "none",
                "maxit": 0,
                "maxit_zv": 0,
                "conv": 0.0,
                "nstate": 0,
                "zvconv": 0.0,
                "nvdav": 0,
            },
        },
        data=FakeData(),
    )

    dump_log(
        mol,
        title="PyOQP: Complete Active Space Configuration Interaction",
        section="fci",
        info={
            "method": "casci",
            "ci_label": "CASCI",
            "active_electrons": 2,
            "active_orbitals": 2,
            "frozen_core": 0,
            "determinants": 4,
            "orbital_source": "rhf",
            "orbital_selection": "sequential",
            "active_orbital_indices": "1,2",
            "core_orbital_indices": "",
            "energies": [-1.0],
            "s2": [0.0],
            "multiplicity": [1],
            "state_average": {
                "energy": -0.75,
                "weights": [0.25, 0.75],
                "roots": [0, 1],
                "root_indices": [0, 1],
            },
            "ci_vector_log": {
                "threshold": 0.05,
                "root_indices": (0,),
                "entries": [
                    {"index": 0, "alpha": "1", "beta": "1", "coefficients": [0.95]},
                ],
            },
        },
    )

    log_text = log_file.read_text(encoding="utf-8")
    assert "PyOQP LOG | ENERGIES AND STATES" in log_text
    assert "PyOQP CASCI CI vectors (abs coeff >= 0.05)" in log_text
    assert "PyOQP CASCI fixed-orbital state average" in log_text
    assert "PyOQP state-average energy:" in log_text
    assert re.search(r"PyOQP\s+det\s+alpha occ\s+beta occ\s+root\s+0", log_text)
    assert re.search(r"PyOQP\s+0\s+1\s+1\s+0.950000", log_text)


def test_dump_log_renders_state_specific_casscf_macroiteration_summary(tmp_path):
    _use_real_oqp_package()
    from oqp.utils.file_utils import dump_log

    class FakeData:
        _data = SimpleNamespace(
            mol_prop=SimpleNamespace(natom=2, charge=0, mult=1),
            control=SimpleNamespace(scftype=1, maxit=50),
            tddft=SimpleNamespace(mult=1),
        )

    log_file = tmp_path / "casscf_macroiteration.log"
    mol = SimpleNamespace(
        log=str(log_file),
        config={
            "input": {"method": "casscf", "basis": "sto-3g", "functional": ""},
            "scf": {
                "forced_attempt": 1,
                "conv": 1.0e-8,
                "incremental": False,
                "diis_type": "cdiis",
                "cdiis_switch": 0.3,
                "vdiis_cdiis_switch": 0.0,
                "vdiis_vshift_switch": 0.0,
                "vshift_cdiis_switch": 0.0,
                "vshift": 0.0,
            },
            "tdhf": {
                "type": "none",
                "maxit": 0,
                "maxit_zv": 0,
                "conv": 0.0,
                "nstate": 0,
                "zvconv": 0.0,
                "nvdav": 0,
            },
        },
        data=FakeData(),
    )

    dump_log(
        mol,
        title="PyOQP: CASSCF Macroiteration Convergence",
        section="casscf_macroiteration",
        info={
            "mode": "state-specific",
            "optimizer": "microiteration",
            "root": 1,
            "max_macro_iterations": 3,
            "n_iterations": 2,
            "n_accepted": 1,
            "converged": False,
            "orbitals_mutated": True,
            "stop_reason": "max_iterations",
            "loop_stop_reason": "max_iterations",
            "root_ambiguous": False,
            "root_swapped": False,
            "energy_initial": -7.0,
            "energy_final": -7.125,
            "energy_decrease": 0.125,
            "gradient_norm_initial": 0.25,
            "gradient_norm_final": 0.03125,
            "final_state_energy": -7.125,
        },
    )

    log_text = log_file.read_text(encoding="utf-8")
    assert "PyOQP: CASSCF Macroiteration Convergence" in log_text
    assert "PyOQP CASSCF mode:                  state-specific" in log_text
    assert "PyOQP CASSCF target root:           1" in log_text
    assert "PyOQP CASSCF macro iterations:      2 / 3" in log_text
    assert "PyOQP CASSCF converged:             no" in log_text
    assert "PyOQP CASSCF final macro energy:    -7.1250000000" in log_text
    assert "PyOQP CASSCF final state energy:    -7.1250000000" in log_text


def test_fci_energy_saves_ci_vector_npz_artifact(tmp_path, monkeypatch):
    _use_real_oqp_package()
    import oqp.library.fci as fci_module

    class FakeFCI(fci_module.FCI):
        def _native_mo_integrals(self):
            h1e = np.zeros((2, 2), dtype=float)
            eri = np.zeros((2, 2, 2, 2), dtype=float)
            plan = fci_module.ActiveSpacePlan(
                norb=2, active=(0, 1), core=(), nelec=(1, 1), metadata={})
            return h1e, eri, plan, 0.0, {
                "determinants": 4,
                "active_electrons": 2,
                "active_orbitals": 2,
                "frozen_core": 0,
                "active_orbital_indices": "1,2",
                "core_orbital_indices": "",
                "orbital_source": "rhf",
                "orbital_selection": "sequential",
            }

    mol = SimpleNamespace(
        project_name="fci-artifact",
        log_path=str(tmp_path),
        config={
            "input": {"functional": ""},
            "scf": {"type": "rhf"},
            "fci": {"nroot": 2, "save_ci_vectors": True},
        },
        data={"nelec_A": 1, "nelec_B": 1},
        mol_energy=SimpleNamespace(nenergy=0.0, energy=None),
        energies=None,
    )
    monkeypatch.setattr(fci_module, "dump_log", lambda *args, **kwargs: None)

    FakeFCI(mol).energy(ref_energy=[0.0])

    artifact = tmp_path / "fci-artifact.fci_ci_vectors.npz"
    assert artifact.exists()
    with np.load(artifact) as payload:
        assert payload["ci_vectors"].shape == (4, 2)
        np.testing.assert_allclose(payload["energies"], mol.data["OQP::FCI_ENERGIES"])
        np.testing.assert_allclose(payload["s2"], mol.data["OQP::FCI_S2"])
        assert payload["multiplicity"].tolist() == mol.data["OQP::FCI_MULTIPLICITY"].tolist()
        assert payload["determinant_count"].tolist() == [4]
        assert payload["active_electrons"].tolist() == [1, 1]


def test_fci_energy_saves_rdm_npz_artifact(tmp_path, monkeypatch):
    _use_real_oqp_package()
    import oqp.library.fci as fci_module

    class FakeFCI(fci_module.FCI):
        def _native_mo_integrals(self):
            h1e = np.zeros((2, 2), dtype=float)
            eri = np.zeros((2, 2, 2, 2), dtype=float)
            plan = fci_module.ActiveSpacePlan(
                norb=2, active=(0, 1), core=(), nelec=(1, 1), metadata={})
            return h1e, eri, plan, 0.0, {
                "determinants": 4,
                "active_electrons": 2,
                "active_alpha_electrons": 1,
                "active_beta_electrons": 1,
                "active_orbitals": 2,
                "frozen_core": 0,
                "active_orbital_indices": "1,2",
                "core_orbital_indices": "",
                "orbital_source": "rhf",
                "orbital_selection": "sequential",
            }

    mol = SimpleNamespace(
        project_name="fci-rdm",
        log_path=str(tmp_path),
        config={
            "input": {"functional": ""},
            "scf": {"type": "rhf"},
            "fci": {"nroot": 2, "save_rdm": True},
        },
        data={"nelec_A": 1, "nelec_B": 1},
        mol_energy=SimpleNamespace(nenergy=0.0, energy=None),
        energies=None,
    )
    monkeypatch.setattr(fci_module, "dump_log", lambda *args, **kwargs: None)

    FakeFCI(mol).energy(ref_energy=[0.0])

    np.testing.assert_allclose(np.trace(mol.data["OQP::FCI_RDM1"]), 2.0)
    assert mol.data["OQP::FCI_RDM1_ROOTS"].shape == (2, 2, 2)
    assert mol.data["OQP::FCI_RDM2_ROOTS"].shape == (2, 2, 2, 2, 2)
    assert mol.data["OQP::FCI_NATURAL_OCCUPATIONS_ROOTS"].shape == (2, 2)
    np.testing.assert_allclose(
        np.trace(mol.data["OQP::FCI_RDM1_ROOTS"], axis1=1, axis2=2),
        [2.0, 2.0],
    )
    np.testing.assert_allclose(
        np.sum(mol.data["OQP::FCI_NATURAL_OCCUPATIONS_ROOTS"], axis=1),
        [2.0, 2.0],
    )

    artifact = tmp_path / "fci-rdm.fci_rdm.npz"
    assert artifact.exists()
    with np.load(artifact) as payload:
        assert payload["rdm1"].shape == (2, 2)
        assert payload["rdm2"].shape == (2, 2, 2, 2)
        assert payload["rdm1_roots"].shape == (2, 2, 2)
        assert payload["rdm2_roots"].shape == (2, 2, 2, 2, 2)
        assert payload["natural_occupations"].shape == (2,)
        assert payload["natural_occupation_roots"].shape == (2, 2)
        np.testing.assert_allclose(
            np.sum(payload["natural_occupation_roots"], axis=1),
            [2.0, 2.0],
        )
        np.testing.assert_allclose(payload["energies"], mol.data["OQP::FCI_ENERGIES"])
        assert payload["root_indices"].tolist() == [0, 1]
        assert payload["active_electrons"].tolist() == [1, 1]


def test_ci_energy_stores_fixed_orbital_state_average(monkeypatch):
    _use_real_oqp_package()
    import oqp.library.fci as fci_module

    class FakeEnergy:
        nenergy = 0.0
        energy = None

    class FakeMol:
        config = {"input": {"functional": ""}, "scf": {"type": "rhf"}}
        data = {"nelec_A": 1, "nelec_B": 1}
        mol_energy = FakeEnergy()
        energies = None

    class FakeCI(fci_module.FCI):
        method_label = "casci"
        data_prefix = "CASCI"
        active_section = "[cas]"
        ci_section = "[ci]"

        def __init__(self, mol):
            self.mol = mol
            self.settings = fci_module.FCISettings(
                nroot=3,
                state_average_enabled=True,
                state_average_weights=(0.25, 0.75),
                state_average_nstate=2,
                state_average_target_roots=(0, 2),
                state_average_equal_weights=False,
            )

        def _native_mo_integrals(self):
            h1e = np.zeros((2, 2), dtype=float)
            eri = np.zeros((2, 2, 2, 2), dtype=float)
            plan = fci_module.ActiveSpacePlan(
                norb=2, active=(0, 1), core=(), nelec=(1, 1), metadata={})
            return h1e, eri, plan, 0.0, {"determinants": 4}

    # The CI solve is one call now (native or the Python driver behind it), so
    # the fabricated roots are injected there rather than at solve_fci.  The
    # driver asks for want_roots, so the double returns the fourth element as
    # well: each returned root's index among the roots the solve computed,
    # which is the identity map for an unfiltered window like this one.
    monkeypatch.setattr(
        fci_module,
        "solve_active_ci",
        lambda *args, **kwargs: (
            np.asarray([-1.0, -0.5, 0.25], dtype=float),
            np.eye(4, 3, dtype=float),
            np.asarray([0.0, 0.0, 2.0], dtype=float),
            np.asarray([0, 1, 2], dtype=np.int64),
        ),
    )
    monkeypatch.setattr(fci_module, "_determinants", lambda norb, nelec: [0, 1, 2, 3])
    monkeypatch.setattr(fci_module, "dump_log", lambda *args, **kwargs: None)

    driver = FakeCI(FakeMol())
    driver.energy(ref_energy=[-0.8])
    stored = driver.mol.data
    assert driver.mol.mol_energy.energy == pytest.approx(-0.0625)
    np.testing.assert_allclose(stored["OQP::CASCI_STATE_AVERAGE_ENERGY"], [-0.0625])
    np.testing.assert_allclose(stored["OQP::CASCI_STATE_AVERAGE_WEIGHTS"], [0.25, 0.75])
    assert stored["OQP::CASCI_STATE_AVERAGE_ROOTS"].tolist() == [0, 2]
    assert stored["OQP::CASCI_STATE_AVERAGE_ROOT_INDICES"].tolist() == [0, 2]


@pytest.mark.parametrize(
    ("settings_kwargs", "energies", "root_indices", "message"),
    [
        (
            {"state_average_equal_weights": False, "state_average_weights": (True, False)},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.weights",
        ),
        (
            {"state_average_equal_weights": False, "state_average_weights": ("0.25", "0.75")},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.weights",
        ),
        (
            {"state_average_target_roots": (False,)},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.target_roots",
        ),
        (
            {"state_average_target_roots": ("0",)},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.target_roots",
        ),
        (
            {"state_average_target_roots": (0.5,)},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.target_roots",
        ),
        (
            {"state_average_target_roots": (0, 0)},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.target_roots",
        ),
        (
            {"state_average_nstate": True},
            np.asarray([-1.0, -0.5]),
            None,
            "state_average.nstate",
        ),
        (
            {},
            np.asarray([True, False], dtype=bool),
            None,
            "state-average energies",
        ),
        (
            {},
            ["-1.0", "-0.5"],
            None,
            "state-average energies",
        ),
        (
            {},
            np.asarray([-1.0, -0.5]),
            np.asarray([True, False], dtype=bool),
            "state_average.root_indices",
        ),
        (
            {},
            np.asarray([-1.0, -0.5]),
            ["0", "1"],
            "state_average.root_indices",
        ),
    ],
)
def test_state_average_payload_rejects_non_native_direct_values(
    settings_kwargs,
    energies,
    root_indices,
    message,
):
    _use_real_oqp_package()
    from oqp.library.fci import FCISettings, _state_average_payload

    settings_options = {
        "nroot": 2,
        "state_average_enabled": True,
        "state_average_nstate": 2,
    }
    settings_options.update(settings_kwargs)
    settings = FCISettings(**settings_options)

    with pytest.raises(ValueError, match=message):
        _state_average_payload(settings, energies, root_indices)


def test_h2_closed_shell_ground_state_is_singlet():
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, fci_spin_diagnostics, solve_fci

    h1e, eri, nelec, ecore, _ = _pyscf_reference("H 0 0 0; H 0 0 0.740")
    _, coeffs = solve_fci(h1e, eri, nelec, ecore=ecore, nroot=1, max_det=1000)
    dets = _determinants(h1e.shape[0], nelec)

    s2, multiplicity = fci_spin_diagnostics(coeffs, dets, h1e.shape[0], nelec)

    assert s2[0] == pytest.approx(0.0, abs=1.0e-10)
    assert multiplicity[0] == 1


def test_h2_full_spin_spectrum_matches_live_pyscf_spin_reference():
    _use_real_oqp_package()
    from oqp.library.fci import _determinants, fci_spin_diagnostics, solve_fci

    nroots = 4
    h1e, eri, nelec, ecore, ref_energies, ref_s2 = _pyscf_fci_roots(
        "H 0 0 0; H 0 0 0.740",
        nroots,
    )
    energies, coeffs = solve_fci(
        h1e,
        eri,
        nelec,
        ecore=ecore,
        nroot=nroots,
        max_det=1000,
        solver="dense",
    )
    dets = _determinants(h1e.shape[0], nelec)

    s2, multiplicity = fci_spin_diagnostics(coeffs, dets, h1e.shape[0], nelec)

    _assert_spin_spectrum_matches(energies, s2, ref_energies, ref_s2)
    assert sorted(multiplicity.tolist()) == [1, 1, 1, 3]


@pytest.mark.parametrize(
    ("atom", "charge", "expected"),
    [
        # Reference totals generated from live PySCF 2.13.0 output.
        ("H 0 0 0; H 0 0 0.740", 0, -1.137283834489),
        ("He 0 0 0; H 0 0 1.4632", 1, -2.826674836464),
    ],
)
def test_dense_fci_solver_matches_live_pyscf_reference(atom, charge, expected):
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e, eri, nelec, ecore, energy_ref = _pyscf_reference(atom, charge=charge)
    assert energy_ref == pytest.approx(expected, abs=5.0e-11)

    energies, _ = solve_fci(
        h1e,
        eri,
        nelec,
        ecore=ecore,
        nroot=1,
        max_det=1000,
    )

    assert energies[0] == pytest.approx(energy_ref, abs=1.0e-10)


def test_input_checker_accepts_fci_energy():
    _use_real_oqp_package()
    from oqp.utils.input_checker import check_input_values

    report = check_input_values(
        {
            "input": {
                "system": "\nH 0 0 0\nH 0 0 0.740",
                "basis": "sto-3g",
                "method": "fci",
                "runtype": "energy",
            },
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "fci": {"nroot": 1, "max_det": 1000},
            "properties": {"scf_prop": []},
        },
        raise_error=False,
        emit=False,
    )

    assert report.ok, report.to_text()


@pytest.mark.parametrize(
    ("name", "system", "atom", "charge", "abs_tol"),
    [
        (
            "h2_fci",
            "\nH 0 0 0\nH 0 0 0.740",
            "H 0 0 0; H 0 0 0.740",
            "0",
            1.0e-8,
        ),
        (
            "h2o_fci",
            "\nO 0.000000 0.000000 0.000000"
            "\nH 0.000000 -0.757000 0.587000"
            "\nH 0.000000 0.757000 0.587000",
            "O 0.000000 0.000000 0.000000;"
            " H 0.000000 -0.757000 0.587000;"
            " H 0.000000 0.757000 0.587000",
            "0",
            5.0e-8,
        ),
    ],
)
def test_openqp_fci_matches_live_pyscf_reference(tmp_path, name, system, atom, charge, abs_tol):
    pytest.importorskip("pyscf")
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end FCI tests")
    from oqp.pyoqp import Runner

    _, _, _, _, energy_ref = _pyscf_reference(atom, charge=int(charge))

    runner = Runner(
        project=name,
        input_file=None,
        log=str(tmp_path / f"{name}.log"),
        input_dict={
            "input": {
                "system": system,
                "charge": charge,
                "basis": "sto-3g",
                "method": "fci",
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
            "properties": {"scf_prop": ""},
            "fci": {"nroot": "1", "max_det": "1000"},
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )

    runner.run(test_mod=True)

    assert np.asarray(runner.mol.energies)[0] == pytest.approx(energy_ref, abs=abs_tol)


def test_openqp_fci_mid_size_active_space_matches_live_pyscf_casci(tmp_path):
    pytest.importorskip("pyscf")
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end FCI tests")
    from oqp.pyoqp import Runner

    system = "\nLi 0.000000 0.000000 0.000000\nH 0.000000 0.000000 1.595000"
    atom = "Li 0.000000 0.000000 0.000000; H 0.000000 0.000000 1.595000"
    energy_ref = _pyscf_casci_reference(atom, active_electrons=4, active_orbitals=4)

    log_file = tmp_path / "lih_cas44_fci.log"
    runner = Runner(
        project="lih_cas44_fci",
        input_file=None,
        log=str(log_file),
        input_dict={
            "input": {
                "system": system,
                "charge": "0",
                "basis": "sto-3g",
                "method": "fci",
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
            "properties": {"scf_prop": ""},
            "fci": {
                "nroot": "1",
                "active_electrons": "4",
                "active_orbitals": "4",
                "max_det": "1000",
            },
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )

    runner.run(test_mod=True)

    assert np.asarray(runner.mol.energies)[0] == pytest.approx(energy_ref, abs=5.0e-8)
    assert runner.mol.data["OQP::FCI_DET_COUNT"].tolist() == [36]

    log_text = log_file.read_text(encoding="utf-8")
    assert re.search(r"PyOQP active electrons:\s+4\b", log_text)
    assert re.search(r"PyOQP active orbitals:\s+4\b", log_text)
    assert re.search(r"PyOQP frozen core orbitals:\s+0\b", log_text)
    assert re.search(r"PyOQP determinant count:\s+36\b", log_text)


def test_openqp_fci_end_to_end_frozen_core_active_space_matches_live_pyscf_casci(tmp_path):
    pytest.importorskip("pyscf")
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end FCI tests")
    from oqp.pyoqp import Runner

    system = (
        "\nO 0.000000 0.000000 0.000000"
        "\nH 0.000000 -0.757000 0.587000"
        "\nH 0.000000 0.757000 0.587000"
    )
    atom = (
        "O 0.000000 0.000000 0.000000;"
        " H 0.000000 -0.757000 0.587000;"
        " H 0.000000 0.757000 0.587000"
    )
    energy_ref = _pyscf_casci_reference(atom, active_electrons=2, active_orbitals=2)

    log_file = tmp_path / "h2o_cas22_fci.log"
    runner = Runner(
        project="h2o_cas22_fci",
        input_file=None,
        log=str(log_file),
        input_dict={
            "input": {
                "system": system,
                "charge": "0",
                "basis": "sto-3g",
                "method": "fci",
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
            "properties": {"scf_prop": ""},
            "fci": {
                "nroot": "1",
                "active_electrons": "2",
                "active_orbitals": "2",
                "max_det": "1000",
            },
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )

    runner.run(test_mod=True)

    assert np.asarray(runner.mol.energies)[0] == pytest.approx(energy_ref, abs=5.0e-8)
    assert runner.mol.data["OQP::FCI_DET_COUNT"].tolist() == [4]

    log_text = log_file.read_text(encoding="utf-8")
    assert re.search(r"PyOQP active electrons:\s+2\b", log_text)
    assert re.search(r"PyOQP active orbitals:\s+2\b", log_text)
    assert re.search(r"PyOQP frozen core orbitals:\s+4\b", log_text)
    assert re.search(r"PyOQP determinant count:\s+4\b", log_text)


def test_openqp_fci_stores_spin_diagnostics_and_logs_multiplicity(tmp_path):
    pytest.importorskip("pyscf")
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end FCI tests")
    from oqp.pyoqp import Runner

    log_file = tmp_path / "h2_fci_spin.log"
    runner = Runner(
        project="h2_fci_spin",
        input_file=None,
        log=str(log_file),
        input_dict={
            "input": {
                "system": "\nH 0 0 0\nH 0 0 0.740",
                "charge": "0",
                "basis": "sto-3g",
                "method": "fci",
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
            "properties": {"scf_prop": ""},
            "fci": {"nroot": "4", "max_det": "1000"},
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )

    runner.run(test_mod=True)

    s2 = runner.mol.data["OQP::FCI_S2"]
    multiplicity = runner.mol.data["OQP::FCI_MULTIPLICITY"]

    assert s2.dtype == np.float64
    assert s2.flags.c_contiguous
    assert s2.shape == (4,)
    assert multiplicity.dtype == np.int64
    assert multiplicity.flags.c_contiguous
    assert multiplicity.shape == (4,)
    assert sorted(multiplicity.tolist()) == [1, 1, 1, 3]
    assert any(np.isclose(s2, 2.0, atol=1.0e-8))

    log_text = log_file.read_text(encoding="utf-8")
    assert "<S^2>" in log_text
    assert "multiplicity" in log_text
    assert "PyOQP note: FCI roots span multiple spin multiplicities" in log_text


def test_active_space_with_dropped_virtuals_matches_casci():
    """Frozen core + active orbitals smaller than the remaining MO space (so at
    least one virtual is dropped) must reproduce a PySCF CASCI energy."""
    pytest.importorskip("pyscf")
    from pyscf import ao2mo, gto, mcscf, scf

    _use_real_oqp_package()
    from oqp.library.fci import FCISettings, _active_space, solve_fci

    mol = gto.Mole()
    mol.atom = "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"
    mol.unit = "Angstrom"
    mol.basis = "sto-3g"
    mol.verbose = 0
    mol.build(cart=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    assert mf.converged

    coeff = mf.mo_coeff
    h1e = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, coeff), coeff.shape[1])

    # CAS(2e,2o) on H2O/STO-3G: 4 frozen-core orbitals, 2 active, >=1 dropped virtual.
    energy_cas = mcscf.CASCI(mf, 2, 2).kernel()[0]

    settings = FCISettings(active_electrons=2, active_orbitals=2)
    h_act, eri_act, nelec_act, ecore_act, meta = _active_space(
        h1e, eri, tuple(mol.nelec), float(mol.energy_nuc()), settings
    )

    assert meta["frozen_core"] == 4
    assert meta["active_orbitals"] == 2
    assert meta["norb"] - meta["frozen_core"] - meta["active_orbitals"] >= 1
    assert meta["determinants"] == 4

    energies, _ = solve_fci(h_act, eri_act, nelec_act, ecore=ecore_act, nroot=1, max_det=1000)
    assert energies[0] == pytest.approx(energy_cas, abs=1.0e-9)


@pytest.mark.parametrize(
    ("target", "value", "message"),
    [
        ("h1e", np.zeros((2, 2), dtype=bool), "h1e"),
        ("h1e", [["0.0", "0.0"], ["0.0", "0.0"]], "h1e"),
        ("h1e", np.zeros((2, 2), dtype=complex), "h1e"),
        ("eri", np.zeros((2, 2, 2, 2), dtype=bool), "eri"),
        ("eri", np.full((2, 2, 2, 2), "0.0"), "eri"),
        ("eri", np.zeros((2, 2, 2, 2), dtype=complex), "eri"),
        ("nelec", (True, False), "nelec"),
        ("nelec", ("1", "1"), "nelec"),
        ("nelec", (1.0, 1.0), "nelec"),
        ("ecore", True, "ecore"),
    ],
)
def test_active_space_rejects_non_native_numeric_inputs(target, value, message):
    _use_real_oqp_package()
    from oqp.library.fci import FCISettings, _active_space

    h1e = np.zeros((2, 2), dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)
    nelec = (1, 1)
    ecore = 0.0
    if target == "h1e":
        h1e = value
    elif target == "eri":
        eri = value
    elif target == "nelec":
        nelec = value
    else:
        ecore = value

    with pytest.raises(ValueError, match=message):
        _active_space(h1e, eri, nelec, ecore, FCISettings(active_electrons=2, active_orbitals=2))


def test_solver_rejects_oversized_determinant_space():
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    norb = 8  # (4 alpha, 4 beta) -> C(8,4)^2 = 4900 determinants
    h1e = np.zeros((norb, norb))
    eri = np.zeros((norb, norb, norb, norb))
    with pytest.raises(ValueError, match="max_det"):
        solve_fci(h1e, eri, (4, 4), max_det=100)


def test_solver_rejects_oversized_dense_hamiltonian():
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    norb = 8  # 4900 determinants -> dense H ~ 192 MiB; cap the budget at 1 MiB.
    h1e = np.zeros((norb, norb))
    eri = np.zeros((norb, norb, norb, norb))
    with pytest.raises(ValueError, match="max_memory"):
        solve_fci(h1e, eri, (4, 4), max_det=100000, max_memory=1)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"nroot": True}, "nroot"),
        ({"nroot": 1.5}, "nroot"),
        ({"max_det": True}, "max_det"),
        ({"max_memory": 1.5}, "max_memory"),
        ({"solver": True}, "solver"),
        ({"solver": []}, "solver"),
        ({"davidson_maxiter": True}, "davidson_maxiter"),
        ({"davidson_subspace": 1.5}, "davidson_subspace"),
        ({"davidson_subspace": True}, "davidson_subspace"),
        ({"ecore": True}, "ecore"),
        ({"ecore": float("nan")}, "ecore"),
        ({"eig_tol": True}, "eig_tol"),
        ({"eig_tol": 0.0}, "eig_tol"),
        ({"eig_tol": float("nan")}, "eig_tol"),
        ({"integral_cutoff": True}, "integral_cutoff"),
        ({"integral_cutoff": -1.0e-12}, "integral_cutoff"),
        ({"integral_cutoff": float("nan")}, "integral_cutoff"),
    ],
)
def test_solve_fci_rejects_invalid_option_literals(kwargs, message):
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    norb = 2
    h1e = np.zeros((norb, norb))
    eri = np.zeros((norb, norb, norb, norb))
    options = {"max_det": 100}
    options.update(kwargs)

    with pytest.raises(ValueError, match=message):
        solve_fci(h1e, eri, (1, 1), **options)


@pytest.mark.parametrize(
    "nelec",
    [
        True,
        "2",
        (True, False),
        ("1", "1"),
        (1.0, 1.0),
        (1, -1),
    ],
)
def test_solve_fci_rejects_non_native_electron_counts(nelec):
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e = np.zeros((2, 2), dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    with pytest.raises(ValueError, match="nelec"):
        solve_fci(h1e, eri, nelec, max_det=100)


@pytest.mark.parametrize(
    ("target", "value", "message"),
    [
        ("h1e", np.zeros((2, 2), dtype=bool), "h1e"),
        ("h1e", [["0.0", "0.0"], ["0.0", "0.0"]], "h1e"),
        ("h1e", np.zeros((2, 2), dtype=complex), "h1e"),
        ("eri", np.zeros((2, 2, 2, 2), dtype=bool), "eri"),
        ("eri", np.full((2, 2, 2, 2), "0.0"), "eri"),
        ("eri", np.zeros((2, 2, 2, 2), dtype=complex), "eri"),
    ],
)
def test_solve_fci_rejects_nonreal_integral_arrays(target, value, message):
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e = np.zeros((2, 2), dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)
    if target == "h1e":
        h1e = value
    else:
        eri = value

    with pytest.raises(ValueError, match=message):
        solve_fci(h1e, eri, (1, 1), max_det=100)


@pytest.mark.parametrize(
    ("override", "bad_path"),
    [
        ({"input": {"runtype": "grad"}}, "input.runtype"),
        ({"scf": {"type": "uhf", "multiplicity": 1}}, "scf.type"),
        ({"input": {"functional": "b3lyp"}}, "input.functional"),
        ({"input": {"d4": True}}, "input.d4"),
        ({"fci": {"nroot": 0}}, "fci.nroot"),
        ({"fci": {"max_det": 0}}, "fci.max_det"),
        ({"fci": {"max_memory": 0}}, "fci.max_memory"),
        ({"fci": {"frozen_core": -1}}, "fci.active_space"),
        ({"fci": {"integral_backend": "pyscf"}}, "fci.integral_backend"),
    ],
)
def test_input_checker_rejects_invalid_fci(override, bad_path):
    _use_real_oqp_package()
    from oqp.utils.input_checker import check_input_values

    config = {
        "input": {
            "system": "\nH 0 0 0\nH 0 0 0.740",
            "basis": "sto-3g",
            "method": "fci",
            "runtype": "energy",
        },
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": 1},
        "fci": {"nroot": 1, "max_det": 1000},
        "properties": {"scf_prop": []},
    }
    for section, values in override.items():
        config.setdefault(section, {}).update(values)

    report = check_input_values(config, raise_error=False, emit=False)

    assert not report.ok, "expected the invalid FCI input to be rejected"
    assert any(d.path == bad_path for d in report.diagnostics), report.to_text()


def test_input_checker_warns_fci_scf_properties_are_reference_only():
    _use_real_oqp_package()
    from oqp.utils.input_checker import check_input_values

    report = check_input_values(
        {
            "input": {
                "system": "\nH 0 0 0\nH 0 0 0.740",
                "basis": "sto-3g",
                "method": "fci",
                "runtype": "energy",
            },
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "fci": {"nroot": 1, "max_det": 1000},
            "properties": {"scf_prop": ["el_mom"]},
        },
        raise_error=False,
        emit=False,
    )

    assert report.ok, report.to_text()
    assert any(
        d.path == "properties.scf_prop" and "RHF reference" in d.message
        for d in report.warnings
    ), report.to_text()


@pytest.mark.parametrize(
    ("atom", "charge", "expected"),
    [
        ("H 0 0 0; H 0 0 0.740", 0, -1.137283834489),
        ("He 0 0 0; H 0 0 1.4632", 1, -2.826674836464),
    ],
)
def test_davidson_solver_matches_live_pyscf_reference(atom, charge, expected):
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e, eri, nelec, ecore, energy_ref = _pyscf_reference(atom, charge=charge)
    assert energy_ref == pytest.approx(expected, abs=5.0e-11)

    energies, _ = solve_fci(
        h1e, eri, nelec, ecore=ecore, nroot=1, max_det=1000, solver="davidson"
    )
    assert energies[0] == pytest.approx(energy_ref, abs=1.0e-9)


def test_dense_and_davidson_agree_for_h2o():
    """H2O/STO-3G full FCI (441 determinants) exercises the real Davidson
    iteration; both solvers must agree with each other and with PySCF."""
    pytest.importorskip("pyscf")
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    atom = (
        "O 0.000000 0.000000 0.000000;"
        " H 0.000000 -0.757000 0.587000;"
        " H 0.000000 0.757000 0.587000"
    )
    h1e, eri, nelec, ecore, energy_ref = _pyscf_reference(atom)

    e_dense, _ = solve_fci(h1e, eri, nelec, ecore=ecore, max_det=5000, solver="dense")
    e_davidson, _ = solve_fci(h1e, eri, nelec, ecore=ecore, max_det=5000, solver="davidson")

    assert e_dense[0] == pytest.approx(energy_ref, abs=1.0e-9)
    assert e_davidson[0] == pytest.approx(energy_ref, abs=1.0e-9)
    assert e_davidson[0] == pytest.approx(e_dense[0], abs=1.0e-9)


def test_davidson_multiroot_matches_dense():
    pytest.importorskip("pyscf")
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e, eri, nelec, ecore, _ = _pyscf_reference("H 0 0 0; H 0 0 0.740")
    e_dense, _ = solve_fci(h1e, eri, nelec, ecore=ecore, nroot=3, max_det=1000, solver="dense")
    e_dav, _ = solve_fci(h1e, eri, nelec, ecore=ecore, nroot=3, max_det=1000, solver="davidson")
    for k in range(3):
        assert e_dav[k] == pytest.approx(e_dense[k], abs=1.0e-9)


def test_solve_fci_passes_davidson_subspace_to_backend(monkeypatch):
    _use_real_oqp_package()
    import oqp.library.fci as fci_module

    h1e = np.zeros((2, 2), dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)
    seen = {}

    def fake_davidson(hamiltonian, diag, nroot, *, tol, max_iter, max_subspace):
        seen["shape"] = hamiltonian.shape
        seen["nroot"] = nroot
        seen["tol"] = tol
        seen["max_iter"] = max_iter
        seen["max_subspace"] = max_subspace
        return np.arange(nroot, dtype=float), np.eye(hamiltonian.shape[0], nroot, dtype=float)

    monkeypatch.setattr(fci_module, "_davidson", fake_davidson)

    energies, coeffs = fci_module.solve_fci(
        h1e,
        eri,
        (1, 1),
        nroot=2,
        max_det=10,
        max_memory=1,
        eig_tol=3.0e-8,
        solver="davidson",
        davidson_maxiter=7,
        davidson_subspace=3,
        ci_section="[ci]",
    )

    assert seen == {
        "shape": (4, 4),
        "nroot": 2,
        "tol": 3.0e-8,
        "max_iter": 7,
        "max_subspace": 3,
    }
    np.testing.assert_allclose(energies, [0.0, 1.0])
    assert coeffs.shape == (4, 2)


def test_solve_fci_rejects_davidson_subspace_smaller_than_nroot():
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    h1e = np.zeros((2, 2), dtype=float)
    eri = np.zeros((2, 2, 2, 2), dtype=float)

    with pytest.raises(ValueError, match=r"\[ci\] davidson_subspace.*nroot"):
        solve_fci(
            h1e,
            eri,
            (1, 1),
            nroot=2,
            solver="davidson",
            davidson_subspace=1,
            ci_section="[ci]",
        )


def test_auto_solver_crosses_dense_memory_wall():
    """With a tiny max_memory the dense path is rejected, but solver=auto must
    fall back to Davidson (sparse) and still return the correct energy."""
    pytest.importorskip("pyscf")
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    atom = (
        "O 0.000000 0.000000 0.000000;"
        " H 0.000000 -0.757000 0.587000;"
        " H 0.000000 0.757000 0.587000"
    )
    h1e, eri, nelec, ecore, energy_ref = _pyscf_reference(atom)  # 441 dets, dense H ~1.5 MiB

    # auto: dense would exceed 1 MiB, so it must use Davidson and succeed.
    energies, _ = solve_fci(
        h1e, eri, nelec, ecore=ecore, max_det=5000, max_memory=1, solver="auto"
    )
    assert energies[0] == pytest.approx(energy_ref, abs=1.0e-9)

    # dense explicitly requested with the same tiny budget must be rejected.
    with pytest.raises(ValueError, match="max_memory"):
        solve_fci(h1e, eri, nelec, ecore=ecore, max_det=5000, max_memory=1, solver="dense")


def test_dense_solver_counts_all_live_eigensolver_matrices():
    _use_real_oqp_package()
    from oqp.library.fci import solve_fci

    # CAS(4e,12o) has 4356 determinants.  Two dense matrices fit below
    # 300 MiB, but the Hamiltonian, symmetric copy and writable LAPACK copy do
    # not; reject before allocating the spin tensor or Hamiltonian.
    h1e = np.zeros((12, 12), dtype=float)
    eri = np.zeros((12, 12, 12, 12), dtype=float)
    with pytest.raises(ValueError, match="matrix working set"):
        solve_fci(
            h1e,
            eri,
            (2, 2),
            nroot=600,
            solver="dense",
            max_det=5000,
            max_memory=300,
        )


@pytest.mark.parametrize(
    ("override", "bad_path"),
    [
        ({"fci": {"solver": "lanczos"}}, "fci.solver"),
        ({"fci": {"davidson_maxiter": 0}}, "fci.davidson_maxiter"),
        ({"fci": {"nroot": 3, "davidson_subspace": 2}}, "fci.davidson_subspace"),
        ({"fci": {"ci_print_threshold": -1.0e-3}}, "fci.ci_print_threshold"),
        ({"fci": {"target_spin": "quartet-ish"}}, "fci.target_spin"),
    ],
)
def test_input_checker_rejects_invalid_solver_options(override, bad_path):
    _use_real_oqp_package()
    from oqp.utils.input_checker import check_input_values

    config = {
        "input": {
            "system": "\nH 0 0 0\nH 0 0 0.740",
            "basis": "sto-3g",
            "method": "fci",
            "runtype": "energy",
        },
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": 1},
        "fci": {"nroot": 1, "max_det": 1000},
        "properties": {"scf_prop": []},
    }
    for section, values in override.items():
        config.setdefault(section, {}).update(values)

    report = check_input_values(config, raise_error=False, emit=False)

    assert not report.ok
    assert any(d.path == bad_path for d in report.diagnostics), report.to_text()


@pytest.mark.parametrize("target_spin", ["triplet", "quintet"])
def test_input_checker_accepts_fci_target_spin_filter(target_spin):
    _use_real_oqp_package()
    from oqp.utils.input_checker import check_input_values

    config = {
        "input": {
            "system": "\nH 0 0 0\nH 0 0 0.740",
            "basis": "sto-3g",
            "method": "fci",
            "runtype": "energy",
        },
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": 1},
        "fci": {"nroot": 2, "max_det": 1000, "target_spin": target_spin},
        "properties": {"scf_prop": []},
    }

    report = check_input_values(config, raise_error=False, emit=False)

    assert report.ok, report.to_text()
    assert not any(d.path == "fci.target_spin" for d in report.diagnostics), report.to_text()


@pytest.mark.parametrize("norb", [1, 2, 3, 4, 6, 8])
def test_spin_orbital_integrals_match_element_reference(norb):
    """The strided spin-block build is a permuted copy of the same spatial
    tensor, so it must reproduce the element-by-element reference bit for bit."""
    _use_real_oqp_package()
    from oqp.library.fci import (
        _spin_orbital_integrals,
        _spin_orbital_integrals_reference,
    )

    rng = np.random.default_rng(1234 + norb)
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4)

    h_ref, g_ref = _spin_orbital_integrals_reference(h1e, eri)
    h_new, g_new = _spin_orbital_integrals(h1e, eri)

    assert np.array_equal(h_new, h_ref)
    assert np.array_equal(g_new, g_ref)


def test_davidson_reaches_the_requested_residual_tolerance():
    """A preconditioned correction has norm ~||r||/gap, so the acceptance test
    on it must be relative.  With an absolute 1e-8 cutoff the iteration used to
    discard every correction once the residual reached ~1e-8 and returned
    silently unconverged, three orders of magnitude short of eig_tol."""
    _use_real_oqp_package()
    from oqp.library.fci import _davidson, _DenseMatrixOperator

    rng = np.random.default_rng(5)
    ndet = 400
    matrix = rng.standard_normal((ndet, ndet)) * 0.02
    matrix = 0.5 * (matrix + matrix.T) + np.diag(np.arange(ndet) * 0.05)

    tol = 1.0e-11
    eigvals, eigvecs = _davidson(
        _DenseMatrixOperator(matrix), np.diag(matrix).copy(), 2,
        tol=tol, max_iter=200,
    )
    residual = np.linalg.norm(
        matrix @ eigvecs - eigvecs * eigvals[None, :], axis=0)
    assert np.max(residual) <= tol

    exact = np.linalg.eigvalsh(matrix)[:2]
    np.testing.assert_allclose(eigvals, exact, atol=1.0e-12)


def test_dense_solver_lowest_roots_shortcut_matches_the_full_solve():
    """solve_fci's dense path extracts only the roots it returns; the result
    must equal the full ndet-root LAPACK solve in energy and in CI vector."""
    _use_real_oqp_package()
    import oqp.library.fci as fci

    rng = np.random.default_rng(41)
    norb = 6
    h1e = rng.standard_normal((norb, norb)) * 0.3
    h1e = h1e + h1e.T - np.diag(np.arange(norb) * 1.5)
    eri = rng.standard_normal((norb,) * 4) * 0.05
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)

    kwargs = dict(nroot=3, solver="dense", max_det=100000, max_memory=2048)
    # CAS(6,6) is 400 determinants, comfortably past the shortcut's floor, so
    # the comparison below is not vacuous.
    taken = []
    saved = fci._lowest_dense_roots

    def _spy(*args, **kw):
        result = saved(*args, **kw)
        taken.append(result is not None)
        return result

    fci._lowest_dense_roots = _spy
    try:
        e_new, v_new = fci.solve_fci(h1e, eri, (3, 3), **kwargs)
    finally:
        fci._lowest_dense_roots = saved
    assert taken == [True]

    fci._lowest_dense_roots = lambda *args, **kw: None
    try:
        e_ref, v_ref = fci.solve_fci(h1e, eri, (3, 3), **kwargs)
    finally:
        fci._lowest_dense_roots = saved

    np.testing.assert_allclose(e_new, e_ref, atol=1.0e-11)
    overlap = np.abs(np.einsum("ij,ij->j", v_new, v_ref))
    np.testing.assert_allclose(overlap, np.ones(3), atol=1.0e-9)


def test_lowest_dense_roots_declines_small_and_wide_windows():
    """The shortcut is a cost optimization, not a semantic change: below the
    measured crossover, and when the window is a large fraction of the
    spectrum, it declines so the caller keeps the full solve."""
    _use_real_oqp_package()
    from oqp.library.fci import _lowest_dense_roots

    rng = np.random.default_rng(9)
    small = rng.standard_normal((100, 100))
    small = 0.5 * (small + small.T)
    assert _lowest_dense_roots(small, 1, 1.0e-10, 100) is None

    wide = rng.standard_normal((400, 400))
    wide = 0.5 * (wide + wide.T)
    assert _lowest_dense_roots(wide, 60, 1.0e-10, 100) is None
