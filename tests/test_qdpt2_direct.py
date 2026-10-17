"""Equivalence and capability tests for the matrix-free direct QDPT2 engine.

The direct engine (:mod:`oqp.library.qdpt2_direct`) must reproduce the dense
reference engine bit-for-bit physics-wise (same semicanonical orbitals, same
diagonal H0, same shift semantics): every energy agrees to <= 1e-10 Eh.  The
dense engine is forced with ``[pt2] engine=dense``; ``auto`` routes the QDPT
family to the direct engine.

The capability test runs a system the dense engine cannot represent at all
(C2H4/STO-3G: the correlated determinant space after freezing is ~8.5e5
determinants, far beyond ``max_det``) and checks physical sanity.
"""
import time

import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


def _run(tmp_path, project, *, system, basis, cas, pt2, ci=None, method="mcqdpt2"):
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
            "input": {"system": system, "charge": "0", "basis": basis,
                      "method": method, "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "60",
                    "forced_attempt": "3", "save_molden": "False"},
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


_H4 = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"
_LIH = "\nLi 0 0 0\nH 0 0 1.600"
_C2H4 = ("\nC 0.000000 0.000000 0.669500"
         "\nC 0.000000 0.000000 -0.669500"
         "\nH 0.000000 0.927035 1.235318"
         "\nH 0.000000 -0.927035 1.235318"
         "\nH 0.000000 0.927035 -1.235318"
         "\nH 0.000000 -0.927035 -1.235318")
_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
          "max_det": "2000"}


def _pt2_on_shared_reference(tmp_path, project, pt2_variants, **kw):
    """Run the RHF reference once, then every ``[pt2]`` variant on it.

    Two independent runs do NOT share a reference once OpenMP is on: the
    two-electron Fock build (``int2.F90``, ``schedule(dynamic,1)`` over shell
    pairs into thread-private Fock copies) sums in a run-dependent order, and
    the SCF stops at ``[scf] conv`` along a slightly different path.  On
    threaded OpenBLAS CI the converged orbitals then differ by up to ~2e-7 and
    the PT2 energies by ~1e-7 between two otherwise identical runs, so an
    engine comparison over two SCFs measures the SCF, not the engines.  Each
    variant here re-enters the PT2 driver (integrals, CASCI, semicanonical
    orbitals, PT2) from the identical reference orbitals.

    Returns the sorted PT2 energies and CASCI reference energies per variant.
    """
    from oqp.library.caspt2_dyall import native_caspt2_energy

    # Variants overwrite one resolved [pt2] block in place; with identical key
    # sets nothing one variant sets can survive into the next.
    assert all(set(v) == set(pt2_variants[0]) for v in pt2_variants), \
        "shared-reference PT2 variants must set the same [pt2] keys"
    runner = _run(tmp_path, project, pt2=pt2_variants[0], **kw)
    mol = runner.mol
    # The driver overwrites VEC_MO_A with semicanonical orbitals; it keeps the
    # reference orbitals it started from in CASPT2_REFERENCE_MO.
    ref_mo = np.array(mol.data["OQP::CASPT2_REFERENCE_MO"], dtype=float, copy=True)
    shape = np.asarray(mol.data["OQP::VEC_MO_A"]).shape

    def collect():
        return (np.sort(np.asarray(mol.energies, dtype=float)),
                np.sort(np.asarray(mol.data["OQP::CASPT2_REFERENCE_ENERGIES"],
                                   dtype=float)),
                np.array(mol.data["OQP::VEC_MO_A"], dtype=float, copy=True))

    results = [collect()]
    for pt2 in pt2_variants[1:]:
        # Poison the stash: only a driver that consumed the restored orbitals
        # and republished them can turn it back into ref_mo.
        mol.data["OQP::CASPT2_REFERENCE_MO"] = np.full_like(ref_mo, np.nan)
        mol.data["OQP::VEC_MO_A"][...] = np.ascontiguousarray(ref_mo.T.reshape(shape))
        mol.config["pt2"].update(pt2)
        native_caspt2_energy(mol)
        np.testing.assert_array_equal(
            np.asarray(mol.data["OQP::CASPT2_REFERENCE_MO"], dtype=float), ref_mo)
        results.append(collect())
    # Same orbitals in -> same semicanonical orbitals and CASCI reference out;
    # anything else means the variants did not run on one reference and the
    # comparison is void.
    for _, ref, mo_sc in results[1:]:
        np.testing.assert_allclose(ref, results[0][1], rtol=0, atol=1.0e-12)
        np.testing.assert_allclose(mo_sc, results[0][2], rtol=0, atol=1.0e-12)
    results = [(e, ref) for e, ref, _ in results]
    return results


def _both_engines(tmp_path, tag, method, ci, pt2_common, **kw):
    (e_dense, _), (e_direct, _) = _pt2_on_shared_reference(
        tmp_path, tag, [{**pt2_common, "engine": "dense"},
                        {**pt2_common, "engine": "direct"}],
        method=method, ci=ci, **kw)
    return e_dense, e_direct


@pytest.mark.parametrize("method,nroot,extra", [
    ("mrmp2", 1, {}),
    ("mcqdpt2", 2, {}),
    ("xmcqdpt2", 2, {}),
    ("mcqdpt2", 3, {"edshft": "0.02"}),
    ("xmcqdpt2", 3, {"imaginary_shift": "0.1"}),
])
def test_direct_engine_matches_dense_h4(tmp_path, method, nroot, extra):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    ci = {"nroot": str(nroot), "target_spin": "singlet"} if nroot > 1 else None
    pt2 = {"reference": "casci", "ipea_shift": "0.0"}
    if nroot > 1:
        pt2.update({"target_roots": ",".join(str(i) for i in range(nroot)),
                    "nroot": str(nroot)})
    pt2.update(extra)
    e_dense, e_direct = _both_engines(
        tmp_path, f"eq_{method}_{nroot}", method, ci, pt2,
        system=_H4, basis="sto-3g", cas=_CAS22)
    np.testing.assert_allclose(e_direct, e_dense, rtol=0, atol=1.0e-10)


def test_direct_engine_matches_dense_frozen_core_lih(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    e_dense, e_direct = _both_engines(
        tmp_path, "eq_lih", "mrmp2", None, {"reference": "casci"},
        system=_LIH, basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
             "max_det": "2000"})
    np.testing.assert_allclose(e_direct, e_dense, rtol=0, atol=1.0e-10)


def test_direct_engine_full_active_zero_correction(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    runner = _run(
        tmp_path, "direct_full_active", method="mrmp2",
        system="\nH 0 0 0\nH 0 0 0.740", basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "0",
             "max_det": "1000"},
        pt2={"reference": "casci", "engine": "direct"},
    )
    e2 = float(np.asarray(
        runner.mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"], dtype=float)[0])
    assert e2 == pytest.approx(0.0, abs=1.0e-12)


def test_direct_engine_c2h4_beyond_dense_reach(tmp_path):
    """C2H4/STO-3G CAS(2,2), SA-3 XMCQDPT2: the frozen-reduced correlated space
    is ~8.5e5 determinants -- the dense engine cannot even allocate it.  The
    direct engine runs it in seconds; check physical sanity of the result."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    t0 = time.time()
    runner = _run(
        tmp_path, "direct_c2h4", method="xmcqdpt2",
        system=_C2H4, basis="sto-3g",
        cas={"active_electrons": "2", "active_orbitals": "2", "frozen_core": "7",
             "max_det": "5000"},
        ci={"nroot": "3", "target_spin": "singlet"},
        pt2={"reference": "casci", "target_roots": "0,1,2", "nroot": "3",
             "edshft": "0.02"},
    )
    wall = time.time() - t0
    ms = np.sort(np.asarray(runner.mol.energies, dtype=float))
    ref = np.sort(np.asarray(
        runner.mol.data["OQP::CASPT2_REFERENCE_ENERGIES"], dtype=float))
    e2 = np.asarray(
        runner.mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"], dtype=float)

    assert ms.shape == (3,)
    assert np.all(np.isfinite(ms))
    # dynamic correlation lowers the ground state; magnitudes stay physical
    assert e2[0] < 0.0
    assert np.all(np.abs(e2) < 1.0)
    # ground state gains correlation energy relative to the CASCI reference
    assert ms[0] < ref[0]
    # capability marker for the log: seconds, not hours
    print(f"\nC2H4 SA-3 XMCQDPT2 direct engine wall time: {wall:.2f} s")


def test_direct_engine_nproc_parallel_matches_serial(tmp_path):
    """[pt2] nproc chunks the reference determinants over worker processes;
    results must be identical to the serial stream.  (With only 4 reference
    determinants the pool path is exercised via the small-chunk threshold, so
    force it by nproc=1-vs-2 on the H4 3-root case.)"""
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    ci = {"nroot": "3", "target_spin": "singlet"}
    pt2 = {"reference": "casci", "target_roots": "0,1,2", "nroot": "3",
           "edshft": "0.02"}
    # nproc=2 with 3 singlet references falls below the 4*nproc chunking
    # threshold and must transparently run serial -- same numbers, no crash.
    (e_serial, _), (e_par, _) = _pt2_on_shared_reference(
        tmp_path, "np", [{**pt2, "nproc": "1"}, {**pt2, "nproc": "2"}],
        method="mcqdpt2", system=_H4, basis="sto-3g", cas=_CAS22, ci=ci)
    np.testing.assert_allclose(e_par, e_serial, rtol=0, atol=1.0e-12)

    # exercise the real pool path directly on the streaming layer
    from oqp.library.qdpt2_direct import _stream_all, _stream_chunk
    rng = np.random.default_rng(7)
    norb, ncore, nact = 6, 1, 2
    h1e = rng.standard_normal((norb, norb)); h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4)
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    eps = rng.standard_normal(norb)
    sup_a = np.array([0b000111, 0b001011, 0b001101, 0b001110,
                      0b010011, 0b010101, 0b010110, 0b011001], dtype=np.uint64)
    sup_b = sup_a.copy()
    serial = _stream_chunk((h1e, eri, eps, norb, ncore, nact, sup_a, sup_b, 0))
    par = _stream_all(h1e, eri, eps, norb, ncore, nact, sup_a, sup_b, 2)
    # compare as merged multisets: sort rows of (ka, kb, src, e0, val) jointly
    def canon(arrs):
        ka, kb, val, e0, src = arrs
        order = np.lexsort((val, e0, src, kb, ka))
        return [a[order] for a in (ka, kb, val, e0, src)]
    for s_arr, p_arr in zip(canon(serial), canon(par)):
        np.testing.assert_array_equal(s_arr, p_arr)


def test_dense_engine_rejects_caspt2_direct_dyall_request(tmp_path):
    """engine=direct with a non-diagonal H0 family must be refused."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    from oqp.library.caspt2_dyall import _caspt2_options
    with pytest.raises(ValueError, match="diagonal"):
        _caspt2_options({"input": {"method": "caspt2"},
                         "pt2": {"engine": "direct", "h0": "dyall"}})


# --------------------------------------------------------------------------- Fortran kernel
def _kernel_available() -> bool:
    if not _backend_available():
        return False
    from oqp import lib
    return hasattr(lib, "qdpt2_stream_kernel")


@pytest.mark.parametrize("method,nroot,extra", [
    ("mrmp2", 1, {}),
    ("mcqdpt2", 3, {"edshft": "0.02"}),
    ("xmcqdpt2", 2, {}),
])
def test_fortran_kernel_matches_python_streaming(tmp_path, method, nroot, extra):
    """The liboqp OpenMP kernel must reproduce the NumPy streaming engine."""
    if not _kernel_available():
        pytest.skip("liboqp lacks qdpt2_stream_kernel (rebuild the core)")

    ci = {"nroot": str(nroot), "target_spin": "singlet"} if nroot > 1 else None
    pt2 = {"reference": "casci", "ipea_shift": "0.0"}
    if nroot > 1:
        pt2.update({"target_roots": ",".join(str(i) for i in range(nroot)),
                    "nroot": str(nroot)})
    pt2.update(extra)
    (e_py, _), (e_f, _) = _pt2_on_shared_reference(
        tmp_path, f"fk_{method}{nroot}",
        [{**pt2, "engine": "direct"}, {**pt2, "engine": "fortran"}],
        method=method, system=_H4, basis="sto-3g", cas=_CAS22, ci=ci)
    np.testing.assert_allclose(e_f, e_py, rtol=0, atol=1.0e-10)


def test_fortran_kernel_threads_match_serial(tmp_path):
    if not _kernel_available():
        pytest.skip("liboqp lacks qdpt2_stream_kernel (rebuild the core)")

    ci = {"nroot": "3", "target_spin": "singlet"}
    pt2 = {"reference": "casci", "target_roots": "0,1,2", "nroot": "3",
           "edshft": "0.02", "engine": "fortran"}
    kw = dict(method="xmcqdpt2", system=_C2H4, basis="sto-3g",
              cas={"active_electrons": "2", "active_orbitals": "2",
                   "frozen_core": "7", "max_det": "5000"}, ci=ci)
    (e1, _), (e4, _) = _pt2_on_shared_reference(
        tmp_path, "fk_t", [{**pt2, "nproc": "1"}, {**pt2, "nproc": "4"}], **kw)
    np.testing.assert_allclose(e4, e1, rtol=0, atol=1.0e-12)
