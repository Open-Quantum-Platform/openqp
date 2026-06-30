"""Tests for the MRSF excited-state analysis & interoperability toolkit.

Pure-Python unit tests (compare harness, descriptors, GTO-evaluator guard) load
the modules directly from their files, so they run without a compiled ``liboqp``
(mirroring tests/test_quantum_fcidump.py). The integration test exercises the
public ``oqp.interop`` API end-to-end and is skipped when the compiled extension
or the benchmark inputs are unavailable.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))


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
                 "CubeExporter", "to_qcschema", "validate_qcschema",
                 "dump_fcidump", "verify_fcidump_fci", "from_openqp",
                 "parse_output", "parse_pyscf_tddft", "compare_results"):
        assert hasattr(I, name), f"oqp.interop missing {name}"


# --------------------------------------------------------------------------
# end-to-end integration (skipped without a build / inputs)
# --------------------------------------------------------------------------
def _have_input(name):
    return os.path.isfile(os.path.join(_ROOT, "benchmark", "inputs", name))


@pytest.mark.skipif(not os.environ.get("OPENQP_ROOT"),
                    reason="needs a built OpenQP (OPENQP_ROOT)")
def test_integration_keystone_qcschema_fcidump(tmp_path):
    pytest.importorskip("oqp")
    if not (_have_input("ch2o_mrsf.inp") and _have_input("h2_rhf_sto3g.inp")):
        pytest.skip("benchmark inputs not present")
    from oqp.pyoqp import Runner
    from oqp.interop import (MRSFExcitedStates, to_qcschema, validate_qcschema,
                             dump_fcidump, verify_fcidump_fci)

    # keystone: reconstructed transition dipole matches OQP's own to <=1e-6
    inp = os.path.join(_ROOT, "benchmark", "inputs", "ch2o_mrsf.inp")
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
    inp2 = os.path.join(_ROOT, "benchmark", "inputs", "h2_rhf_sto3g.inp")
    r2 = Runner(project="h2_t", input_file=inp2,
                log=str(tmp_path / "h2.log"), silent=1, usempi=False)
    r2.run()
    path = str(tmp_path / "h2.FCIDUMP")
    dump_fcidump(path, r2.mol)
    ver = verify_fcidump_fci(path, r2.mol)
    assert ver["diff"] < 1e-8
