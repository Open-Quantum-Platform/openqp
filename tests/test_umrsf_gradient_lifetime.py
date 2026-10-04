"""Same-process lifetime of the UMRSF gradient's persistent native buffers.

``tdhf_umrsf_gradient.F90`` keeps module-level ``save`` scratch for the
mean-field contraction (``umrsf_mf_dens``/``umrsf_mf_fout`` sized by the packed
basis, ``umrsf_mf_fxa``.. sized by the basis) and the XC kernel state
(``xc_refa``/``xc_moa``/``xc_molgrid``).  Those survive between calls, so a
fresh-process example proves nothing about a second calculation with a
different basis.  This test drives several gradients through one interpreter
with increasing and decreasing basis sizes, an HF run (no XC state) between two
DFT runs, an energy-only run that never enters the gradient, a different target
state, and an immediate repeat, and requires every gradient to reproduce the
committed example reference (AGENTS.md, resource-lifetime rules).

The reference tolerance is the one used by ``openqp --run_tests`` for the same
examples; the repeat tolerance is far tighter because the same binary on the
same input must reuse (or correctly resize) the buffers without contamination.
"""

import configparser
import json
import os
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "UMRSF-TDDFT"

REF_TOL = 1.0e-6      # Ha/Bohr against the committed example references
REPEAT_TOL = 1.0e-10  # same binary, same input, later in the process


def _drop_stubbed_oqp_modules():
    """tests/test_umrsf_energy_regression.py loads the input checker against stub
    ``oqp``/``oqp.utils``/``oqp.utils.mpi_utils`` modules; in one unforked pytest
    session those stubs would shadow the installed package here."""
    import sys
    for name in [n for n in list(sys.modules) if n == "oqp" or n.startswith("oqp.")]:
        module = sys.modules[name]
        if getattr(module, "__file__", None) is None and getattr(module, "__spec__", None) is None:
            del sys.modules[name]


def _backend_available() -> bool:
    _drop_stubbed_oqp_modules()
    try:
        import oqp
        from oqp import lib
    except Exception:
        return False
    if not bool(getattr(oqp, "BACKEND_AVAILABLE", False)):
        return False
    return hasattr(lib, "tdhf_umrsf_gradient")


requires_backend = pytest.mark.skipif(
    not _backend_available(), reason="liboqp has no tdhf_umrsf_gradient (rebuild needed)"
)


def _example(name, **overrides):
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(EXAMPLES / f"{name}.inp")
    config = {section: dict(parser[section]) for section in parser.sections()}
    for section, values in overrides.items():
        config.setdefault(section, {}).update(values)
    return config


def _reference(name):
    with open(EXAMPLES / f"{name}.json", encoding="utf-8") as handle:
        data = json.load(handle)
    return float(data["energy"]), np.asarray(data["grad"], dtype=float)


def _run(workdir, tag, config):
    _drop_stubbed_oqp_modules()
    from oqp.pyoqp import Runner

    runner = Runner(project=tag, input_file=None, log=str(workdir / f"{tag}.log"),
                    input_dict=config, silent=1, usempi=False)
    runner.run(test_mod=True)
    state = int(config["properties"].get("grad", 0)) if config["input"]["runtype"] == "grad" else None
    grad = None if state is None else np.asarray(runner.mol.grads[state], dtype=float).ravel()
    # energies[0] is the UHF reference energy that the example JSON stores as "energy".
    return float(runner.mol.energies[0]), grad, runner


@requires_backend
def test_umrsf_gradient_buffers_survive_basis_changes_in_one_process(tmp_path):
    h2co = "H2CO_BHHLYP_UMRSFTDDFT_GRAD"          # 6-31G*, 34 bf, BHHLYP (XC state on)
    h2co_big = "H2CO_BHHLYP_UMRSFTDDFT_GRAD_CCPVDZ"  # cc-pVDZ, 38 bf, spherical d
    h2 = "H2_HF_UMRSFTDDFT_GRAD"                   # HF, no XC state, 2 atoms

    e_ref, g_ref = _reference(h2co)
    e_big_ref, g_big_ref = _reference(h2co_big)
    e_h2_ref, g_h2_ref = _reference(h2)

    # 1. fresh allocation
    e1, g1, _ = _run(tmp_path, "step1_h2co", _example(h2co))
    assert abs(e1 - e_ref) < 1e-7 and np.max(np.abs(g1 - g_ref)) < REF_TOL

    # 2. larger basis: packed and square scratch must grow
    # The tighter request exercises the bounded fallback on this build. A reduced
    # solve that passes directly on another platform is accepted below.
    e2, g2, _ = _run(tmp_path, "step2_ccpvdz", _example(h2co_big, tdhf={"zvconv": "5e-10"}))
    assert abs(e2 - e_big_ref) < 1e-7 and np.max(np.abs(g2 - g_big_ref)) < REF_TOL
    z_log = (tmp_path / "step2_ccpvdz.log").read_text()
    # The gate bounds the absolute residual by sqrt(zvconv), as the other
    # OpenQP Z-vector solvers do.
    full_residuals = [float(value) for value in re.findall(
        r"full coupled Z residual =\s*([\d.E+-]+)", z_log
    )]
    z_bound = 5e-10 ** 0.5
    assert full_residuals and full_residuals[-1] <= z_bound
    if full_residuals[0] > z_bound:
        assert "UMRSF Z auto fallback: dense full-block solve" in z_log
        assert len(full_residuals) == 2

    # 3. much smaller system and HF: scratch shrinks, XC state must be released
    e3, g3, _ = _run(tmp_path, "step3_h2_hf", _example(h2))
    assert abs(e3 - e_h2_ref) < 1e-7 and np.max(np.abs(g3 - g_h2_ref)) < REF_TOL

    # 4. energy-only run never reaches the gradient scratch (early exit before the Z-vector)
    e4, _, _ = _run(tmp_path, "step4_energy_only",
                    _example(h2co, input={"runtype": "energy"}, properties={"grad": "0"}))
    assert abs(e4 - e_ref) < 1e-7

    # 5. back to the first basis after HF: XC state and scratch are rebuilt at the old size
    e5, g5, _ = _run(tmp_path, "step5_h2co_again", _example(h2co))
    assert abs(e5 - e1) < REPEAT_TOL and np.max(np.abs(g5 - g1)) < REPEAT_TOL

    # 6. different target state on the same basis (state-count identity, not basis identity)
    e6, g6, _ = _run(tmp_path, "step6_h2co_s1", _example(h2co, properties={"grad": "2"}))
    assert abs(e6 - e1) < REPEAT_TOL                    # same reference SCF
    assert np.max(np.abs(g6 - g1)) > 1e-4               # a genuinely different state

    # 7. immediate repeat with unchanged dimensions reuses the buffers
    e7, g7, _ = _run(tmp_path, "step7_repeat", _example(h2co))
    assert abs(e7 - e1) < REPEAT_TOL and np.max(np.abs(g7 - g1)) < REPEAT_TOL
