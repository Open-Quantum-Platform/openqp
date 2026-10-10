"""Regression test for the tagarray dangling-view readback artifact.

``OQPData.__getitem__`` returns zero-copy numpy views into Fortran tagarray
storage.  The storage is freed by ``ffi.gc(..., oqp_clean)`` when the OQPData
handle is garbage-collected, so reading through a *temporary* Molecule --
``Runner(...).mol.data[key]`` -- used to hand back a view whose memory was
freed as soon as the expression finished, yielding intermittent garbage
(denormals like 1e-323).  This was the "uninitialized readback artifact" noted
in the multistate CASPT2 tests.  The fix pins the owning handle on the
returned array (``_TagArrayView._oqp_owner``); this test locks that behavior.
"""
import gc

import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


def _fresh_runner(tmp_path, project):
    from oqp.pyoqp import Runner

    return Runner(
        project=project,
        input_file=None,
        log=str(tmp_path / f"{project}.log"),
        input_dict={
            "input": {"system": "\nH 0 0 0\nH 0 0 0.740", "charge": "0",
                      "basis": "sto-3g", "method": "hf", "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "50",
                    "save_molden": "False"},
            "properties": {"scf_prop": ""},
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )


def test_view_from_temporary_molecule_survives_gc(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    def make_view():
        runner = _fresh_runner(tmp_path, "ta_lifetime")
        runner.run(test_mod=True)
        runner.mol.data["OQP::PROBE_LIFETIME"] = np.array([3.14, 2.71, 1.61])
        # returned object is a zero-copy view; the Runner/mol become garbage
        # as soon as this function returns
        return runner.mol.data["OQP::PROBE_LIFETIME"]

    view = make_view()
    gc.collect()
    # touch unrelated allocations to encourage reuse of any freed pages
    _ = [np.random.default_rng(1).standard_normal(4096) for _ in range(8)]
    gc.collect()
    np.testing.assert_allclose(np.asarray(view, dtype=float),
                               [3.14, 2.71, 1.61], atol=1.0e-15)


def test_asarray_and_slices_keep_the_owner_alive(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    runner = _fresh_runner(tmp_path, "ta_lifetime2")
    runner.run(test_mod=True)
    runner.mol.data["OQP::PROBE_LIFETIME2"] = np.arange(6, dtype=np.float64)
    sliced = np.asarray(runner.mol.data["OQP::PROBE_LIFETIME2"], dtype=float)[2:5]
    del runner
    gc.collect()
    np.testing.assert_allclose(sliced, [2.0, 3.0, 4.0], atol=1.0e-15)


def test_inplace_write_through_view_still_works(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    runner = _fresh_runner(tmp_path, "ta_lifetime3")
    runner.run(test_mod=True)
    runner.mol.data["OQP::PROBE_INPLACE"] = np.zeros(4)
    runner.mol.data["OQP::PROBE_INPLACE"][...] = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(
        np.asarray(runner.mol.data["OQP::PROBE_INPLACE"], dtype=float),
        [1.0, 2.0, 3.0, 4.0], atol=1.0e-15)
