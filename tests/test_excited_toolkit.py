"""Tests for the MRSF excited-state analysis & interoperability toolkit.

Pure-Python unit tests (compare harness, descriptors, GTO-evaluator guard) load
the modules directly from their files, so they run without a compiled ``liboqp``
(mirroring tests/test_quantum_fcidump.py). The integration test exercises the
public ``oqp.interop`` API end-to-end and is skipped when the compiled extension
or the fixture inputs are unavailable.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
_INPUT_DIR = os.path.join(_HERE, "fixtures", "excited_toolkit")


def _load_file(modname, relpath):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(_ROOT, relpath))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


compare = _load_file("_xt_compare", "pyoqp/oqp/interop/compare.py")
descriptors = _load_file("_xt_descriptors", "pyoqp/oqp/analysis/descriptors.py")
gto_grid = _load_file("_xt_gto_grid", "pyoqp/oqp/analysis/gto_grid.py")


# --------------------------------------------------------------------------
# compare harness  (Codex finding: fail on missing/mismatched states)
# --------------------------------------------------------------------------
def test_compare_equal_arrays_pass():
    rows, ok = compare.compare_results(
        {"excitation_energies_ev": [4.0, 8.0]},
        {"excitation_energies_ev": [4.0001, 8.0]},
        {"excitation_energies_ev": 1e-3})
    assert ok and rows[0][3] == "PASS"


def test_compare_length_mismatch_is_failure():
    # external reported fewer states than OQP -> must FAIL, not prefix-PASS
    rows, ok = compare.compare_results(
        {"excitation_energies_ev": [4.0, 8.0, 9.0]},
        {"excitation_energies_ev": [4.0]},
        {"excitation_energies_ev": 1e-3})
    assert not ok
    assert rows[0][3] == "FAIL" and "length" in rows[0][4]


def test_compare_scalar_and_missing():
    rows, ok = compare.compare_results(
        {"scf_energy_ha": -1.0, "x": None},
        {"scf_energy_ha": -1.0 + 1e-9},
        {"scf_energy_ha": 1e-6, "x": 1e-6})
    by = {r[0]: r for r in rows}
    assert by["scf_energy_ha"][3] == "PASS"
    assert by["x"][3] == "N/A"


# --------------------------------------------------------------------------
# descriptors
# --------------------------------------------------------------------------
def test_participation_ratio_bounds():
    assert abs(descriptors.participation_ratio([1.0, 0.0, 0.0]) - 1.0) < 1e-12
    assert abs(descriptors.participation_ratio([0.5, 0.5]) - 2.0) < 1e-12
    assert abs(descriptors.participation_ratio([0.25] * 4) - 4.0) < 1e-12


class _PhysicalPairStates:
    """Two-root fixture whose auxiliary-reference amplitudes must not be used."""

    nstates = 2
    nbf = 2
    S = np.eye(2)
    C = np.eye(2)
    energies = np.array([-10.0, -9.8])

    class _Mol:
        @staticmethod
        def get_atoms():
            return np.array([8, 6])

    mol = _Mol()

    @staticmethod
    def tdm_mo(i, j):
        assert (i, j) == (0, 1)
        # Row = particle AO (C p_z), column = hole AO (O p_x).
        return np.array([[0.0, 0.0], [1.0, 0.0]])

    @classmethod
    def tdm_ao(cls, i, j):
        return cls.tdm_mo(i, j)

    @staticmethod
    def amplitude_matrix(_n):
        raise AssertionError("physical-root analysis used the high-spin reference")


class _TwoAO:
    nbf = 2
    ao_atom = np.array([0, 1])
    # O p_x is in the xy plane; C p_z is perpendicular.
    ao_index = [(0, (1, 0, 0), 1.0), (1, (0, 0, 1), 1.0)]
    coords = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])


def test_fragment_ct_uses_physical_root_pair_and_orientation():
    result = descriptors.fragment_ct_matrix(
        _PhysicalPairStates(), _TwoAO(), 1, [[0], [1]], ref=0)
    np.testing.assert_allclose(result["Omega"], [[0.0, 1.0], [0.0, 0.0]])
    assert result["le_fraction"] == 0.0
    assert result["ct_fraction"] == 1.0
    assert result["pair"] == (0, 1)


def test_fragment_ct_rejects_incomplete_or_overlapping_partition():
    states, ao = _PhysicalPairStates(), _TwoAO()
    with pytest.raises(ValueError, match="do not assign"):
        descriptors.fragment_ct_matrix(states, ao, 1, [[0]])
    with pytest.raises(ValueError, match="more than one"):
        descriptors.fragment_ct_matrix(states, ao, 1, [[0], [0, 1]])


# --------------------------------------------------------------------------
# GTO evaluator guard  (Codex finding: spherical shells unsupported)
# --------------------------------------------------------------------------
class _FakeData:
    def __init__(self, basis):
        self._b = basis

    def get_basis(self):
        return self._b


class _FakeMol:
    def __init__(self, basis, sysv):
        self.data = _FakeData(basis)
        self._sys = np.asarray(sysv, dtype=float)

    def get_system(self):
        return self._sys


def _one_d_shell_basis(nbf):
    # a single d shell on one atom: Cartesian=6 AOs, spherical=5 AOs
    return {"nbf": nbf, "nsh": 1, "centers": np.array([0]), "angs": np.array([2]),
            "ncontr": np.array([1]), "alpha": np.array([1.0]), "coef": np.array([1.0])}


def test_gto_grid_rejects_spherical_basis():
    mol = _FakeMol(_one_d_shell_basis(5), [0.0, 0.0, 0.0])     # 5 == spherical d
    with pytest.raises(NotImplementedError):
        gto_grid.AOBasis(mol)


def test_gto_grid_accepts_cartesian_basis():
    mol = _FakeMol(_one_d_shell_basis(6), [0.0, 0.0, 0.0])     # 6 == Cartesian d
    ao = gto_grid.AOBasis(mol)
    assert ao.nbf == 6 and len(ao.ao_index) == 6


# --------------------------------------------------------------------------
# public API surface
# --------------------------------------------------------------------------
def test_interop_public_api_shape():
    oqp = pytest.importorskip("oqp")           # needs the compiled extension
    import oqp.interop as I
    for name in ("MRSFExcitedStates", "nto_excitation", "attachment_detachment",
                 "participation_ratio", "tozer_lambda", "fragment_ct_matrix",
                 "analyze_mrsf_transition", "infer_molecular_plane",
                 "CubeExporter", "to_qcschema", "validate_qcschema",
                 "dump_fcidump", "verify_fcidump_fci", "from_openqp",
                 "parse_output", "parse_pyscf_tddft", "compare_results"):
        assert hasattr(I, name), f"oqp.interop missing {name}"


# --------------------------------------------------------------------------
# end-to-end integration (skipped without a build / inputs)
# --------------------------------------------------------------------------
def _have_input(name):
    return os.path.isfile(os.path.join(_INPUT_DIR, name))


def _input_path(name):
    return os.path.join(_INPUT_DIR, name)


@pytest.mark.skipif(not os.environ.get("OPENQP_ROOT"),
                    reason="needs a built OpenQP (OPENQP_ROOT)")
def test_integration_keystone_qcschema_fcidump(tmp_path):
    pytest.importorskip("oqp")
    if not (_have_input("ch2o_mrsf.inp") and _have_input("h2_rhf_sto3g.inp")):
        pytest.skip("excited-toolkit fixture inputs not present")
    from oqp.pyoqp import Runner
    from oqp.interop import (MRSFExcitedStates, to_qcschema, validate_qcschema,
                             dump_fcidump, verify_fcidump_fci)

    # keystone: reconstructed transition dipole matches OQP's own to <=1e-6
    inp = _input_path("ch2o_mrsf.inp")
    r = Runner(project="ch2o_t", input_file=inp,
               log=str(tmp_path / "ch2o.log"), silent=1, usempi=False)
    r.run()
    st = MRSFExcitedStates(r.mol)
    worst = 0.0
    for i in range(st.nstates):
        for j in range(i + 1, st.nstates):
            worst = max(worst, float(np.max(np.abs(
                st.transition_dipole(i, j) - st.dip_oqp[:, i, j]))))
    assert worst < 1e-6

    # qcschema validates and the excited-state arrays are aligned (S0->Sn)
    payload = to_qcschema(r.mol, states=st)
    res = validate_qcschema(payload)
    assert res.success
    ex = res.extras["oqp"]
    assert len(ex["excitation_energies_ev"]) == st.nstates - 1
    assert len(ex["oscillator_strengths"]) == len(ex["excitation_energies_ev"])
    assert len(ex["total_state_energies_hartree"]) == st.nstates

    # FCIDUMP delegated to oqp.quantum, verified against PySCF FCI
    pytest.importorskip("pyscf")
    inp2 = _input_path("h2_rhf_sto3g.inp")
    r2 = Runner(project="h2_t", input_file=inp2,
                log=str(tmp_path / "h2.log"), silent=1, usempi=False)
    r2.run()
    path = str(tmp_path / "h2.FCIDUMP")
    dump_fcidump(path, r2.mol)
    ver = verify_fcidump_fci(path, r2.mol)
    assert ver["diff"] < 1e-8


@pytest.mark.skipif(not os.environ.get("OPENQP_ROOT"),
                    reason="needs a built OpenQP (OPENQP_ROOT)")
def test_integration_state_analysis_on_real_mrsf_roots(tmp_path):
    """Run the state analysis on genuine MRSF densities, not a hand-made 2x2.

    The unit tests inject a 2x2 matrix with ``tdm_ao == tdm_mo``, so they cannot
    see the tagarray reshape, the root ordering, the MO->AO transform or the
    nonorthogonal overlap.  Formaldehyde exercises all of them.
    """
    pytest.importorskip("oqp")
    if not _have_input("ch2o_mrsf.inp"):
        pytest.skip("excited-toolkit fixture inputs not present")
    from oqp.pyoqp import Runner
    from oqp.interop import (AOBasis, MRSFExcitedStates, analyze_mrsf_transition,
                             infer_molecular_plane)

    r = Runner(project="ch2o_sa", input_file=_input_path("ch2o_mrsf.inp"),
               log=str(tmp_path / "ch2o_sa.log"), silent=1, usempi=False)
    r.run()
    st = MRSFExcitedStates(r.mol)
    ao = AOBasis(r.mol)
    whole = [list(range(len(np.asarray(r.mol.get_atoms()).ravel())))]
    # C, O | H, H -- two fragments, so Omega is 2x2 and its off-diagonal
    # elements are what the hole/particle index convention actually controls.
    # A single fragment would make every transpose assertion below vacuous.
    split = [[0, 1], [2, 3]]

    # The fixture geometry is planar in yz, so the pi axis is x.  This is pure
    # geometry and holds regardless of the electronic structure.
    normal = infer_molecular_plane(ao, np.asarray(r.mol.get_atoms(), dtype=int))
    assert abs(abs(float(normal[0])) - 1.0) < 1e-6

    off_diagonal_seen = False
    for n in range(1, st.nstates):
        rep = analyze_mrsf_transition(st, ao, n, whole)
        ct = rep["fragment_ct"]
        assert ct["well_defined"]
        # Omega is a partition of the squared 1-TDM norm.
        assert abs(ct["Omega"].sum() - ct["total"]) < 1e-8 * max(ct["total"], 1.0)
        # A single fragment cannot support charge transfer.
        assert abs(rep["le_fraction"] - 1.0) < 1e-8
        fractions = rep["orbital_character"]["fractions"]
        assert rep["orbital_character"]["classified"]
        assert abs(sum(fractions.values()) - 1.0) < 1e-8
        assert all(np.isfinite(v) for v in fractions.values())

        # gamma^{n->0} is the transpose of gamma^{0->n}, so hole and particle
        # swap and Omega must transpose with them.  A stray transpose in the
        # MO->AO path or in the fragment loop breaks this -- but only shows up
        # once Omega has off-diagonal weight, hence the two-fragment split.
        pair = analyze_mrsf_transition(st, ao, n, split)["fragment_ct"]
        back = analyze_mrsf_transition(st, ao, 0, split, ref=n)["fragment_ct"]
        np.testing.assert_allclose(back["Omega"], pair["Omega"].T,
                                   rtol=1e-8, atol=1e-12)
        assert abs(back["le_fraction"] - pair["le_fraction"]) < 1e-8
        off_diagonal_seen |= bool(
            np.max(np.abs(pair["Omega"] - np.diag(np.diag(pair["Omega"]))))
            > 1e-6 * pair["total"])
        assert abs(rep["energy_gap"]
                   + analyze_mrsf_transition(st, ao, 0, whole,
                                             ref=n)["energy_gap"]) < 1e-10

    # Guard the guard: if every Omega were diagonal the transpose check above
    # would have proved nothing.
    assert off_diagonal_seen

    # Formaldehyde's low MRSF roots are the classic n->pi* / pi->pi* pair: the
    # acceptor is pi* in both, which the plane projection must reproduce.
    s1 = analyze_mrsf_transition(st, ao, 1, whole)["orbital_character"]
    to_pi_star = sum(v for k, v in s1["fractions"].items() if k.endswith("->pi*"))
    assert to_pi_star > 0.5


def test_fragment_ct_accepts_numpy_fragment_definitions():
    # Fragments are naturally built with numpy helpers, and `not ndarray`
    # raises "truth value of an array ... is ambiguous" for any fragment with
    # more than one atom (Codex review of #290).
    states, ao = _PhysicalPairStates(), _TwoAO()
    result = descriptors.fragment_ct_matrix(
        states, ao, 1, [np.array([0]), np.array([1])], ref=0)
    np.testing.assert_allclose(result["Omega"], [[0.0, 1.0], [0.0, 0.0]])
    # ... including a multi-atom fragment, which is the case that raised.
    with pytest.raises(ValueError, match="more than one"):
        descriptors.fragment_ct_matrix(
            states, ao, 1, [np.array([0, 1]), np.array([1])], ref=0)
