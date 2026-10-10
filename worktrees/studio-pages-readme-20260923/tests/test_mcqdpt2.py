"""End-to-end and unit tests for the native QDPT family (MRMP2/MCQDPT2/XMCQDPT2).

The QDPT family re-labels and re-routes the validated determinant PT2 machinery
of ``oqp.library.caspt2_dyall`` under the GAMESS conventions:

* ``mrmp2``    == single-state ``caspt2`` with the (default) diagonal-Fock H0;
* ``mcqdpt2``  == multistate QDPT on the SINGLE-SET state-averaged canonical
  orbitals with one common diag(eps) H0 and ket-state-specific E0 (Nakano);
  by definition it never takes the multi-set OpenMolcas-style MS-CASPT2 branch;
* ``xmcqdpt2`` == Granovsky's extended QDPT (XZERO): the SA-Fock model-space
  rotation + the same common H0 -- the identical code path as ``xms-caspt2``.
* ``[pt2] edshft`` is the GAMESS ISA denominator regularization d -> d + s/d,
  mutually exclusive with the real/imaginary level shifts.

These tests pin the label/routing equivalences, the single-set-vs-multi-set
gap, the ISA intruder taming, and the option-validation surface.
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


def _run_pt2(tmp_path, project, *, system, basis, cas, pt2, ci=None, method="caspt2"):
    from oqp.pyoqp import Runner

    ci_block = {
        "nroot": "1",
        "solver": "dense",
        "eig_tol": "1.0e-10",
        "integral_backend": "native",
        "integral_cutoff": "5.0e-11",
    }
    if ci:
        ci_block.update(ci)
    runner = Runner(
        project=project,
        input_file=None,
        log=str(tmp_path / f"{project}.log"),
        input_dict={
            "input": {
                "system": system,
                "charge": "0",
                "basis": basis,
                "method": method,
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
            "cas": cas,
            "ci": ci_block,
            "pt2": pt2,
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )
    runner.run(test_mod=True)
    return runner


_H4_LINEAR = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"
_H2 = "\nH 0 0 0\nH 0 0 0.740"
_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"}


def _corrections(runner):
    return np.asarray(runner.mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"], dtype=float)


# --------------------------------------------------------------------------- label/routing equivalences
def test_mrmp2_equals_single_state_caspt2_fock(tmp_path):
    """MRMP2 is the QDPT label for the diagonal-Fock single-state PT2: identical
    numbers to method=caspt2 (default h0=fock)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    kw = dict(system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
              pt2={"reference": "casci", "ipea_shift": "0.0"})
    e_caspt2 = _run_pt2(tmp_path, "qdpt_ref_caspt2", method="caspt2", **kw).mol.energies
    e_mrmp2 = _run_pt2(tmp_path, "qdpt_mrmp2", method="mrmp2", **kw).mol.energies
    np.testing.assert_allclose(e_mrmp2, e_caspt2, atol=1.0e-10)


def test_xmcqdpt2_equals_xms_caspt2(tmp_path):
    """XMCQDPT2 (XZERO) is the same single-set rotated path as xms-caspt2 here;
    only the labels differ."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    kw = dict(system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
              ci={"nroot": "2", "target_spin": "singlet"},
              pt2={"reference": "casci", "target_roots": "0,1", "nroot": "2",
                   "ipea_shift": "0.0"})
    e_xms = np.sort(np.asarray(
        _run_pt2(tmp_path, "qdpt_xms", method="xms-caspt2", **kw).mol.energies, dtype=float))
    e_xmc = np.sort(np.asarray(
        _run_pt2(tmp_path, "qdpt_xmcqdpt2", method="xmcqdpt2", **kw).mol.energies, dtype=float))
    np.testing.assert_allclose(e_xmc, e_xms, atol=1.0e-10)


def test_mcqdpt2_is_single_set_not_multiset(tmp_path):
    """MCQDPT2 must stay on the single-set SA-orbital path (GAMESS convention),
    differing from the multi-set ms-caspt2 by the bounded single-set/multi-set
    gap for well-separated H4 states -- and its effective Hamiltonian must be
    consistent (symmetric, diagonal = SS energies, eigenvalues = reported)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    kw = dict(system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
              ci={"nroot": "2", "target_spin": "singlet"},
              pt2={"reference": "casci", "target_roots": "0,1", "nroot": "2",
                   "ipea_shift": "0.0"})
    r_mc = _run_pt2(tmp_path, "qdpt_mcqdpt2", method="mcqdpt2", **kw)
    e_mc = np.sort(np.asarray(r_mc.mol.energies, dtype=float))
    e_ms = np.sort(np.asarray(
        _run_pt2(tmp_path, "qdpt_mscaspt2", method="ms-caspt2", **kw).mol.energies, dtype=float))

    # single-set vs multi-set: same physics, different orbital sets -> small,
    # nonzero, bounded gap (identical would mean the multi-set branch was taken)
    gap = np.max(np.abs(e_mc - e_ms))
    assert 1.0e-6 < gap < 5.0e-3

    heff = np.asarray(r_mc.mol.data["OQP::CASPT2_EFFECTIVE_HAMILTONIAN"], dtype=float)
    ss = np.asarray(r_mc.mol.data["OQP::CASPT2_SS_ENERGIES"], dtype=float)
    np.testing.assert_allclose(heff, heff.T, atol=1.0e-12)
    np.testing.assert_allclose(np.sort(np.diag(heff)), np.sort(ss), atol=1.0e-10)
    np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(heff)), e_mc, atol=1.0e-10)


def test_mcqdpt2_full_active_space_zero_correction(tmp_path):
    """A CAS spanning the whole orbital space leaves no external perturbers:
    every MCQDPT2 state correction is exactly zero."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    runner = _run_pt2(
        tmp_path, "qdpt_full_active", method="mcqdpt2",
        system=_H2, basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "0",
             "max_det": "1000"},
        ci={"nroot": "2", "target_spin": "singlet"},
        pt2={"reference": "casci", "target_roots": "0,1", "nroot": "2"},
    )
    np.testing.assert_allclose(_corrections(runner), 0.0, atol=1.0e-12)


# --------------------------------------------------------------------------- ISA (edshft)
def _run_3root(tmp_path, project, method, pt2_extra):
    pt2 = {"reference": "casci", "target_roots": "0,1,2", "nroot": "3",
           "ipea_shift": "0.0"}
    pt2.update(pt2_extra)
    return _run_pt2(
        tmp_path, project, method=method,
        system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
        ci={"nroot": "3", "target_spin": "singlet"},
        pt2=pt2,
    )


def test_edshft_isa_tames_the_single_set_intruder(tmp_path):
    """The high H4 singlet is an intruder on the single-set path; the GAMESS ISA
    shift (edshft) regularizes it like the imaginary shift does, while leaving
    the well-behaved roots essentially unchanged."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    e2_raw = _corrections(_run_3root(tmp_path, "qdpt_isa_raw", "mcqdpt2", {}))
    e2_isa = _corrections(_run_3root(tmp_path, "qdpt_isa_on", "mcqdpt2", {"edshft": "0.02"}))

    assert abs(e2_raw[2]) > 0.1          # bare intruder blows up
    assert abs(e2_isa[2]) < 0.1          # ISA brings it back to PT2 magnitude
    np.testing.assert_allclose(e2_isa[:2], e2_raw[:2], atol=5.0e-3)


def test_edshft_zero_is_a_strict_noop(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    kw = dict(system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
              pt2={"reference": "casci"})
    e_plain = _run_pt2(tmp_path, "qdpt_noop_a", method="mrmp2", **kw).mol.energies
    kw2 = dict(system=_H4_LINEAR, basis="sto-3g", cas=_CAS22,
               pt2={"reference": "casci", "edshft": "0.0"})
    e_zero = _run_pt2(tmp_path, "qdpt_noop_b", method="mrmp2", **kw2).mol.energies
    np.testing.assert_allclose(e_zero, e_plain, atol=1.0e-12)


# --------------------------------------------------------------------------- option-surface units (no backend needed)
def _options(method, pt2):
    from oqp.library.caspt2_dyall import _caspt2_options
    return _caspt2_options({"input": {"method": method}, "pt2": pt2})


def test_qdpt_method_mapping_units():
    o = _options("mrmp2", {})
    assert (o.family, o.variant, o.h0) == ("qdpt", "caspt2", "fock")
    o = _options("mcqdpt2", {})
    assert (o.family, o.variant) == ("qdpt", "ms")
    o = _options("xmcqdpt2", {})
    assert (o.family, o.variant) == ("qdpt", "xms")
    o = _options("xms-caspt2", {})
    assert (o.family, o.variant) == ("caspt2", "xms")


def test_qdpt_rejects_dyall_h0_units():
    with pytest.raises(ValueError, match="diagonal-Fock"):
        _options("mcqdpt2", {"h0": "dyall"})


def test_edshft_validation_units():
    assert _options("mcqdpt2", {"edshft": "0.02"}).edshft == pytest.approx(0.02)
    with pytest.raises(ValueError, match="mutually exclusive"):
        _options("mcqdpt2", {"edshft": "0.02", "imaginary_shift": "0.1"})
    with pytest.raises(ValueError, match=">= 0"):
        _options("mcqdpt2", {"edshft": "-0.01"})


def test_qdpt_rejects_strong_contraction_units():
    with pytest.raises(ValueError):
        _options("mrmp2", {"contraction": "strong"})
