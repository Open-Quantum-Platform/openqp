"""End-to-end tests for the clean uncontracted Dyall CASPT2 module.

These drive ``method=caspt2`` through the full dispatch (oqp.library.caspt2_dyall)
and check the physics rather than diagnostic sidecars:

* the second-order energy lowers the reference (E2 < 0);
* the stored reference equals the CASCI energy (embedding self-check);
* a full active space leaves no external perturbers (E2 == 0);
* the H4/STO-3G result reproduces the OpenMolcas-anchored validated value.

OpenMolcas 24.10 CIonly CASPT2 reference (H4/STO-3G, CAS(2,2), IPEA=0, no shift):
    E(CASCI) = -2.1147334918   E2 = -0.0184952337   E(CASPT2) = -2.1332287255
With the ``fock`` (CASPT2) H0 this all-valence case matches the OpenMolcas E2 to
~0.5 uEh.  ``test_caspt2_h0_fock_scope`` additionally pins the quantitative
all-valence agreement and the known, bounded atomic-core over-correlation
(LiH/HF/H2O); see the module docstring of ``oqp.library.caspt2_dyall``.
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


def _run_caspt2(tmp_path, project, *, system, basis, cas, pt2, ci=None, method="caspt2"):
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


def test_pt2_nproc_zero_keeps_the_automatic_sentinel():
    from oqp.library.caspt2_dyall import _caspt2_options

    assert _caspt2_options({"pt2": {"nproc": "0"}}).nproc == 0
    assert _caspt2_options({"pt2": {"nproc": "-3"}}).nproc == 0


def _caspt2_energy(runner):
    return float(np.asarray(runner.mol.data["OQP::CASPT2_ENERGIES"], dtype=float)[0])


def _reference_energy(runner):
    return float(np.asarray(runner.mol.data["OQP::CASPT2_REFERENCE_ENERGIES"], dtype=float)[0])


def _correction(runner):
    return float(np.asarray(runner.mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"], dtype=float)[0])


def test_caspt2_dyall_h4_sto3g_matches_validated_reference(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    runner = _run_caspt2(
        tmp_path,
        "caspt2_dyall_h4_sto3g",
        system=_H4_LINEAR,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"},
        pt2={"variant": "caspt2", "reference": "casci",
             "ipea_shift": "0.0", "level_shift": "0.0", "imaginary_shift": "0.0"},
    )

    e_ref = _reference_energy(runner)
    e2 = _correction(runner)
    e_caspt2 = _caspt2_energy(runner)

    # CASCI reference matches OpenMolcas CIonly to micro-Hartree.
    assert e_ref == pytest.approx(-2.1147334918, abs=1.0e-6)
    # The Fock-H0 (CASPT2) E2 matches OpenMolcas CIonly CASPT2 to <= 1e-4:
    # OpenMolcas E2 = -0.0184952337; OpenQP ~ -0.0184957 (sub-uHartree agreement).
    assert e2 == pytest.approx(-0.0184952337, abs=1.0e-4)
    # Total is the sum, and the correction lowers the energy.
    assert e_caspt2 == pytest.approx(e_ref + e2, abs=1.0e-10)
    assert e_caspt2 < e_ref
    assert runner.mol.energies[0] == pytest.approx(e_caspt2, abs=1.0e-12)


def test_nevpt2_h0_dyall_matches_pyscf(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    # h0=dyall selects Dyall's H0 (the reference is an exact H0 eigenstate) = NEVPT2.
    runner = _run_caspt2(
        tmp_path,
        "nevpt2_h4_sto3g",
        system=_H4_LINEAR,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"},
        pt2={"reference": "casci", "h0": "dyall", "ipea_shift": "0.0"},
    )
    # PySCF strongly-contracted NEVPT2 E2 = -0.0195923188 for this system.
    assert _correction(runner) == pytest.approx(-0.0195923188, abs=1.0e-4)


def test_caspt2_dyall_full_active_has_zero_correction(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    # STO-3G H2 has two orbitals; CAS(2,2) with no frozen core spans the whole
    # space, so there is no external first-order space and E2 must vanish.
    runner = _run_caspt2(
        tmp_path,
        "caspt2_dyall_h2_full_active",
        system=_H2,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "0", "max_det": "1000"},
        pt2={"variant": "caspt2", "reference": "casci"},
    )

    assert _correction(runner) == pytest.approx(0.0, abs=1.0e-12)
    assert _caspt2_energy(runner) == pytest.approx(_reference_energy(runner), abs=1.0e-12)


def test_caspt2_dyall_level_shift_reduces_correction_magnitude(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    base = _run_caspt2(
        tmp_path,
        "caspt2_dyall_h4_noshift",
        system=_H4_LINEAR,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"},
        pt2={"variant": "caspt2", "reference": "casci", "level_shift": "0.0"},
    )
    shifted = _run_caspt2(
        tmp_path,
        "caspt2_dyall_h4_shift",
        system=_H4_LINEAR,
        basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"},
        pt2={"variant": "caspt2", "reference": "casci", "level_shift": "0.1"},
    )

    e2_base = _correction(base)
    e2_shift = _correction(shifted)
    # A positive real level shift damps the (negative) correction toward zero.
    assert e2_base < 0.0
    assert e2_shift < 0.0
    assert abs(e2_shift) < abs(e2_base)


_MS_CAS = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "2000"}


def _run_multistate(tmp_path, project, method):
    return _run_caspt2(
        tmp_path,
        project,
        method=method,
        system=_H4_LINEAR,
        basis="sto-3g",
        cas=_MS_CAS,
        ci={"nroot": "2", "target_spin": "singlet"},
        pt2={"reference": "casci", "target_roots": "0,1", "nroot": "2",
             "ipea_shift": "0.0", "level_shift": "0.0"},
    )


# OpenMolcas 24.10 CIonly MS-CASPT2 reference, H4/STO-3G CAS(2,2), Inactive=1,
# Ras2=2, IPEA=0, no shift (the multi-set state-specific-Fock construction):
#   CASCI     -2.11473349 / -1.40952700 / -0.94839445
#   SS-CASPT2 -2.13322873 / -1.43832941 / -1.11662642
#   MS-CASPT2 -2.13341747 / -1.43832941 / -1.11643767  (roots 1&3 mix, root 2 inert)
# ``method=ms-caspt2`` with the default ``h0=fock`` uses the multi-set (per-state
# orbitals) construction and reproduces these to <= 1.8 uEh.
_H4_MS3_MOLCAS_CASCI = np.array([-2.11473349, -1.40952700, -0.94839445])
_H4_MS3_MOLCAS_SS = np.array([-2.13322873, -1.43832941, -1.11662642])
_H4_MS3_MOLCAS_MS = np.array([-2.13341747, -1.43832941, -1.11643767])


def test_ms_caspt2_effective_hamiltonian_is_consistent(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    runner = _run_multistate(tmp_path, "ms_caspt2_h4_2state", "ms-caspt2")
    ms = np.asarray(runner.mol.data["OQP::CASPT2_ENERGIES"], dtype=float)
    ss = np.asarray(runner.mol.data["OQP::CASPT2_SS_ENERGIES"], dtype=float)
    heff = np.asarray(runner.mol.data["OQP::CASPT2_EFFECTIVE_HAMILTONIAN"], dtype=float)

    assert ms.shape == (2,) and ss.shape == (2,) and heff.shape == (2, 2)
    # effective Hamiltonian is symmetric and its diagonal is the per-state CASPT2
    np.testing.assert_allclose(heff, heff.T, atol=1.0e-12)
    np.testing.assert_allclose(np.sort(np.diag(heff)), np.sort(ss), atol=1.0e-10)
    # diagonalizing H_eff reproduces the reported MS energies (ascending)
    np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(heff)), np.sort(ms), atol=1.0e-10)
    # two well-separated states barely mix: MS ~= SS
    np.testing.assert_allclose(np.sort(ms), np.sort(ss), atol=1.0e-3)
    # the ground MS root lies at/below the lowest SS root (variational mixing)
    assert ms.min() <= ss.min() + 1.0e-10
    # the multi-set construction reproduces the OpenMolcas MS-CASPT2 roots (the two
    # lowest H4 states; well-separated, so MS == SS == the OpenMolcas SS-CASPT2 values)
    np.testing.assert_allclose(np.sort(ms), np.sort(_H4_MS3_MOLCAS_SS[:2]), atol=1.0e-4)


def test_xms_caspt2_differs_from_ms_by_the_multiset_single_set_gap(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    # Read the canonical energies from mol.energies (the value the log prints).  The
    # OQP::CASPT2_ENERGIES tagarray export from the multistate path has a known,
    # intermittent uninitialized-readback artifact on the post-2026-06-29 merged
    # backend (the method itself is correct -- mol.energies is always right); see the
    # WF_methods notes for that follow-up.
    ms = np.sort(np.asarray(
        _run_multistate(tmp_path, "xms_cmp_ms", "ms-caspt2").mol.energies, dtype=float))
    xms = np.sort(np.asarray(
        _run_multistate(tmp_path, "xms_cmp_xms", "xms-caspt2").mol.energies, dtype=float))
    # MS-CASPT2 (h0=fock) now uses the multi-set state-specific-Fock construction
    # (per-state orbitals), which reproduces OpenMolcas; XMS-CASPT2 keeps the
    # single-set (state-averaged-orbital) path.  For these well-separated H4 states
    # the two references therefore differ by a small, bounded multi-set/single-set
    # gap (~1.5 mEh), not by the mixing the old single-set MS shared with XMS.
    gap = np.max(np.abs(np.sort(xms) - np.sort(ms)))
    assert 1.0e-4 < gap < 5.0e-3
    # (the multi-set MS-vs-OpenMolcas-multiset validation is in
    # test_ms_caspt2_h4_3root_matches_molcas_multiset; not re-checked here against the
    # SS reference, which is a different quantity from the multi-set MS roots.)


def test_ms_caspt2_h4_3root_matches_molcas_multiset(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    # H4/STO-3G CAS(2,2), 3 singlet roots: the multi-set MS-CASPT2 effective
    # Hamiltonian reproduces OpenMolcas MS-CASPT2 to <= 1.8 uEh (roots 0 and 2 mix,
    # root 1 is inert).  This is the headline multi-set validation.
    runner = _run_3state(tmp_path, "ms_caspt2_h4_3state", 0.0)
    ms = np.sort(np.asarray(runner.mol.data["OQP::CASPT2_ENERGIES"], dtype=float))
    ref = np.sort(np.asarray(runner.mol.data["OQP::CASPT2_REFERENCE_ENERGIES"], dtype=float))

    assert ms.shape == (3,)
    # the CASCI reference reproduces OpenMolcas
    np.testing.assert_allclose(ref, np.sort(_H4_MS3_MOLCAS_CASCI), atol=1.0e-6)
    # the MS-CASPT2 roots reproduce OpenMolcas multi-set MS-CASPT2 to <= 0.1 mEh
    np.testing.assert_allclose(ms, np.sort(_H4_MS3_MOLCAS_MS), atol=1.0e-4)
    # the mixing genuinely splits the outer roots away from the SS-CASPT2 values
    assert abs(ms[0] - np.sort(_H4_MS3_MOLCAS_SS)[0]) > 1.0e-5


def _run_3state(tmp_path, project, imaginary_shift, method="ms-caspt2"):
    return _run_caspt2(
        tmp_path,
        project,
        method=method,
        system=_H4_LINEAR,
        basis="sto-3g",
        cas=_MS_CAS,
        ci={"nroot": "3", "target_spin": "singlet"},
        pt2={"reference": "casci", "variant": method, "target_roots": "0,1,2", "nroot": "3",
             "ipea_shift": "0.0", "imaginary_shift": str(imaginary_shift)},
    )


def test_imaginary_shift_removes_high_root_intruder(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    # The high singlet (root 2) of H4/STO-3G is an intruder in the uncontracted
    # external space: its E(2) blows up without a shift.  This is exercised on the
    # single-set path (XMS-CASPT2) -- the multi-set MS-CASPT2 path uses per-state
    # orbitals that already tame this root, so the bare intruder is only visible in
    # the single-set construction.
    e2_raw = np.asarray(
        _run_3state(tmp_path, "ms3_noshift", 0.0, method="xms-caspt2")
        .mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"],
        dtype=float,
    )
    e2_shift = np.asarray(
        _run_3state(tmp_path, "ms3_shift", 0.1, method="xms-caspt2")
        .mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"],
        dtype=float,
    )

    # without a shift the intruder state over-corrects badly (|E2| > 0.1 Eh)
    assert abs(e2_raw[2]) > 0.1
    # the imaginary shift brings it back to a physical CASPT2 magnitude
    assert abs(e2_shift[2]) < 0.1
    # the well-behaved low roots are barely affected by the shift
    np.testing.assert_allclose(e2_shift[:2], e2_raw[:2], atol=5.0e-3)


# --------------------------------------------------------------------------- H0=fock accuracy scope
# The ``fock`` (CASPT2) H0 is quantitative versus OpenMolcas CIonly CASPT2
# (IPEA=0, no shift).  With the default PT2 frozen core (``[pt2] frozen=auto`` --
# the deep atomic cores OpenMolcas also freezes) EVERY case -- all-valence and
# atomic-core -- matches OpenMolcas to a few uEh.  The earlier 0.03-0.23 mEh
# "atomic-core over-correlation" was an all-electron-vs-frozen-core mismatch:
# ``frozen=0`` correlates the deep core and reproduces the old all-electron value.
_LIH = "\nLi 0 0 0\nH 0 0 1.600"
_HF = "\nH 0 0 0\nF 0 0 0.917"
_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"

# (name, system, basis, frozen_core, OpenMolcas E2) -- default frozen core matches all.
_SCOPE_CASES = [
    ("h4", _H4_LINEAR, "sto-3g", "1", -0.0184952337),
    ("lih", _LIH, "sto-3g", "1", -0.0125985108),
    ("hf", _HF, "sto-3g", "4", -0.0165325858),
    ("h2o", _H2O, "sto-3g", "4", -0.0343680407),
]


@pytest.mark.parametrize(
    "name,system,basis,fc,molcas_e2",
    _SCOPE_CASES,
    ids=[c[0] for c in _SCOPE_CASES],
)
def test_caspt2_h0_fock_scope(tmp_path, name, system, basis, fc, molcas_e2):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    runner = _run_caspt2(
        tmp_path,
        f"caspt2_scope_{name}",
        system=system,
        basis=basis,
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": fc, "max_det": "5000"},
        pt2={"variant": "caspt2", "reference": "casci", "ipea_shift": "0.0"},
    )
    # With the default frozen core, both all-valence and atomic-core agree with
    # OpenMolcas CIonly CASPT2 to <= 22 uEh.
    assert abs(_correction(runner) - molcas_e2) <= 22.0e-6


def test_caspt2_frozen_core_is_the_atomic_core_fix(tmp_path):
    """The default PT2 frozen core (auto = deep cores) closes the atomic-core case to
    OpenMolcas; ``[pt2] frozen=0`` correlates the deep core and reproduces the old
    all-electron over-correlated value -- proving the residual was a frozen-core
    convention, not an H0 error."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    cas = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "4", "max_det": "5000"}
    frozen = _run_caspt2(tmp_path, "hf_frozen_auto", system=_HF, basis="sto-3g", cas=cas,
                         pt2={"variant": "caspt2", "reference": "casci"})
    allcorr = _run_caspt2(tmp_path, "hf_frozen_none", system=_HF, basis="sto-3g", cas=cas,
                          pt2={"variant": "caspt2", "reference": "casci", "frozen": "0"})

    # default (deep core frozen) matches OpenMolcas CIonly CASPT2
    assert abs(_correction(frozen) - (-0.0165325858)) <= 1.0e-5
    # all-electron (frozen=0) reproduces the old over-correlated value, ~30 uEh more bound
    assert _correction(allcorr) == pytest.approx(-0.0165622899, abs=3.0e-6)
    assert _correction(allcorr) < _correction(frozen) - 2.0e-5


def test_ipea_shift_no_op_at_zero_and_monotone(tmp_path):
    """The IPEA shift (Ghigo-Roos-Malmqvist 2004) is a strict no-op at 0 (the default)
    and biases the active H0 diagonal so E2 becomes monotonically less negative for
    IPEA > 0."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    cas = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "4", "max_det": "5000"}

    def hf(ip):
        return _correction(_run_caspt2(
            tmp_path, f"hf_ipea_{int(ip * 100)}", system=_HF, basis="sto-3g", cas=cas,
            pt2={"variant": "caspt2", "reference": "casci", "ipea_shift": str(ip)}))

    e0, e25, e50 = hf(0.0), hf(0.25), hf(0.5)
    # strict no-op at 0: equals the validated (frozen-core) CASPT2 value
    assert e0 == pytest.approx(-0.0165325587, abs=2.0e-6)
    # monotone less negative with increasing IPEA
    assert e0 < e25 < e50


def test_sc_nevpt2_strong_contraction_matches_pyscf(tmp_path):
    """SC-NEVPT2 (`[pt2] contraction=strong`, RDM-based strong contraction) closes
    the CAS(4,4) flavor gap: it matches PySCF SC-NEVPT2 (on CASCI orbitals), whereas
    the uncontracted Dyall path sits ~0.4 mEh below it."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    cas = {"active_electrons": "4", "active_orbitals": "4", "frozen_core": "0", "max_det": "20000"}
    strong = _run_caspt2(
        tmp_path, "sc_nevpt2_h4_631g_cas44", system=_H4_LINEAR, basis="6-31g", cas=cas,
        pt2={"reference": "casci", "h0": "dyall", "contraction": "strong"},
    )
    uncontracted = _run_caspt2(
        tmp_path, "nevpt2_h4_631g_cas44_unc", system=_H4_LINEAR, basis="6-31g", cas=cas,
        pt2={"reference": "casci", "h0": "dyall", "contraction": "none"},
    )
    # PySCF SC-NEVPT2 on CASCI orbitals = -0.0275992508 (independent reference).
    assert _correction(strong) == pytest.approx(-0.0275992508, abs=1.0e-5)
    # The uncontracted Dyall path is the validated all-valence engine, unchanged.
    assert _correction(uncontracted) == pytest.approx(-0.0280053140, abs=1.0e-6)
    # Strong contraction genuinely differs from uncontracted (the ~0.4 mEh flavor gap).
    assert _correction(strong) > _correction(uncontracted) + 3.0e-4


# --------------------------------------------------------------- H0 representations
# The zeroth-order Hamiltonian is stored in whichever exact representation is
# cheapest (diagonal for Fock H0, block-diagonal over the core+virtual occupation
# patterns for Dyall H0).  The explicit dense ndet x ndet matrix with one full
# eigendecomposition per root remains in the module as the fallback; these tests
# pin the cheap representations against it.
def _fast_and_dense(tmp_path, project, **kwargs):
    """Run a case as shipped, then again with the H0 fast paths disabled."""
    import oqp.library.caspt2_dyall as cd
    from oqp.library.fci import _symmetric_eigh

    fast = _run_caspt2(tmp_path, project + "_fast", **kwargs)

    def _full_eigh(self):
        blk = self.dense[np.ix_(self.external, self.external)]
        w, U = _symmetric_eigh(0.5 * (blk + blk.T))
        return w, ((slice(None), U),)

    saved_diag = cd._diagonal_zeroth_order
    saved_eig = cd._ZerothOrder.external_eigenbasis
    cd._diagonal_zeroth_order = lambda *a, **k: None
    cd._ZerothOrder.external_eigenbasis = _full_eigh
    try:
        dense = _run_caspt2(tmp_path, project + "_dense", **kwargs)
    finally:
        cd._diagonal_zeroth_order = saved_diag
        cd._ZerothOrder.external_eigenbasis = saved_eig
    return fast, dense


_CAS44 = {"active_electrons": "4", "active_orbitals": "4", "frozen_core": "0",
          "max_det": "20000"}


def test_diagonal_h0_matches_the_dense_zeroth_order_matrix(tmp_path):
    """``h0=fock`` H0 is a one-electron operator that is already diagonal in the MO
    basis, hence exactly diagonal in the determinant basis: it needs neither the
    dense matrix nor a denominator diagonalization.  Same energy either way."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    fast, dense = _fast_and_dense(
        tmp_path, "caspt2_diag_h0", system=_H4_LINEAR, basis="6-31g", cas=_CAS44,
        pt2={"reference": "casci", "h0": "fock"},
    )
    assert _correction(fast) == pytest.approx(_correction(dense), abs=1.0e-12)
    assert _caspt2_energy(fast) == pytest.approx(_caspt2_energy(dense), abs=1.0e-12)


def test_dyall_h0_block_partition_matches_the_full_eigendecomposition(tmp_path):
    """Dyall H0 cannot move an electron in or out of a core or virtual orbital, so
    ``H0_ext`` is exactly block-diagonal over the core+virtual occupation patterns.
    Diagonalizing those blocks must reproduce the single full solve."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    fast, dense = _fast_and_dense(
        tmp_path, "caspt2_dyall_blocks", system=_H4_LINEAR, basis="6-31g", cas=_CAS44,
        pt2={"reference": "casci", "h0": "dyall", "contraction": "none"},
    )
    assert _correction(fast) == pytest.approx(_correction(dense), abs=1.0e-12)
    assert _caspt2_energy(fast) == pytest.approx(_caspt2_energy(dense), abs=1.0e-12)


def test_occupation_blocks_close_the_dyall_zeroth_order_operator(tmp_path):
    """The partition is only sound if nothing at all lives outside the blocks."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASPT2 tests")

    import oqp.library.caspt2_dyall as cd

    seen = {}
    original = cd._build_operators

    def _spy(*args, **kwargs):
        out = original(*args, **kwargs)
        h0op = out[3]
        if h0op.dense is not None and h0op.blocks is not None:
            blk = h0op.dense[np.ix_(h0op.external, h0op.external)]
            off = 0.0
            for g in h0op.blocks:
                rows = np.abs(blk[g, :])
                rows[:, g] = 0.0
                off = max(off, float(rows.max()))
            seen.update(off=off, nblocks=len(h0op.blocks),
                        scale=float(np.max(np.abs(blk))))
        return out

    cd._build_operators = _spy
    try:
        _run_caspt2(
            tmp_path, "caspt2_dyall_blockcheck", system=_H4_LINEAR, basis="6-31g",
            cas=_CAS44, pt2={"reference": "casci", "h0": "dyall", "contraction": "none"},
        )
    finally:
        cd._build_operators = original

    assert seen, "the Dyall path did not build a dense H0 with a block partition"
    assert seen["nblocks"] > 1
    # Exactly 0.0 in practice; the tolerance only guards a LAPACK-level denormal.
    assert seen["off"] <= 1.0e-12 * seen["scale"]
