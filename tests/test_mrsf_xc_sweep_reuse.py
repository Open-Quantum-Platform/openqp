"""Pin fused XC derivatives against separate sweeps on relaxed MRSF densities."""
import os

import numpy as np
import pytest


@pytest.mark.parametrize("functional,multiplicity,state", [
    ("bhhlyp", 1, 1), ("bhhlyp", 1, 3), ("bhhlyp", 3, 1), ("slater", 1, 1),
])
def test_fused_xc_sweep_uses_production_relaxed_density(tmp_path, functional, multiplicity, state):
    import oqp
    from oqp.pyoqp import Runner

    cfg = {
        "input": {"system": "O .07 .02 .01\nH .95 .06 -.03\nH -.25 .93 .10",
                  "basis": "6-31g*", "functional": functional, "method": "tdhf",
                  "runtype": "grad", "omp_threads": os.environ.get("OMP_NUM_THREADS", "1")},
        "scf": {"type": "rohf", "multiplicity": "3", "conv": "1e-10", "save_molden": "false"},
        "guess": {"type": "huckel", "save_mol": "false"},
        "tdhf": {"type": "mrsf", "nstate": "3", "multiplicity": str(multiplicity),
                 "conv": "1e-9", "zvconv": "1e-9"},
        "properties": {"grad": str(state)},
    }
    runner = Runner(project="xc_sweep", input_file=None, input_dict=cfg,
                    log=str(tmp_path / "run.log"), silent=1, usempi=False)
    runner.run(test_mod=True)
    assert np.all(np.isfinite(runner.mol.grads[state]))
    error = oqp.ffi.new("double *")
    moving_norm = oqp.ffi.new("double *")
    oqp.lib.mrsf_xc_sweep_selftest(runner.mol.data._data, error, moving_norm)
    assert np.isfinite(error[0]) and error[0] < 1e-11
    assert moving_norm[0] > 1e-8  # the omitted-sweep error must be observable
