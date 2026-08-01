"""End-to-end tests for the PT2-family central-difference numerical gradients.

The native PT2 family (mrmp2/mcqdpt2/xmcqdpt2/caspt2/ms-/xms-caspt2) is
energy-only in the kernels; ``oqp.library.pt2_numgrad`` supplies
central-difference gradients over Cartesian displacements and
``Gradient.gradient`` (single_point.py) dispatches PT2 methods there, which
makes runtype=grad and the gradient-driven optimizers (optimize, penalty
MECI) work for PT2.

One wiring note: these tests construct the Runner with runtype=energy and
flip ``mol.config['input']['runtype']`` afterwards -- ``Runner.run`` reads
the runtype from the live config, so this exercises the dispatch without
re-parsing an input file per case.
"""
import os
import time

import numpy as np
import pytest

BOHR_TO_ANGSTROM = 0.529177210903

_H2_TEMPLATE = "\nH 0 0 0\nH 0 0 %.10f"
_H4_LINEAR = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


def _make_runner(tmp_path, project, system, *, method="mrmp2", nroot=1,
                 frozen_core="0"):
    """PT2 Runner in the proven test configuration (RHF/CASCI reference).

    Constructed with runtype=energy so the input checker passes; callers flip
    the runtype on the live config for grad/optimize/meci runs.
    """
    from oqp.pyoqp import Runner

    ci = {"nroot": str(max(nroot, 1)), "solver": "dense", "eig_tol": "1.0e-10",
          "integral_backend": "native", "integral_cutoff": "5.0e-11"}
    pt2 = {"reference": "casci"}
    if nroot > 1:
        ci["target_spin"] = "singlet"
        pt2.update({"nroot": str(nroot),
                    "target_roots": ",".join(str(i) for i in range(nroot))})
    return Runner(
        project=project,
        input_file=None,
        log=str(tmp_path / f"{project}.log"),
        input_dict={
            "input": {"system": system, "charge": "0", "basis": "sto-3g",
                      "method": method, "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "60",
                    "save_molden": "False"},
            "properties": {"scf_prop": ""},
            "cas": {"active_electrons": "2", "active_orbitals": "2",
                    "frozen_core": frozen_core, "max_det": "2000"},
            "ci": ci,
            "pt2": pt2,
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )


def _h2_energy(tmp_path, tag, r_angstrom):
    """Independent single-point H2 mrmp2 energy at bond length r (Angstrom)."""
    runner = _make_runner(tmp_path, f"e_{tag}", _H2_TEMPLATE % r_angstrom)
    runner.run(test_mod=True)
    return float(runner.mol.energies[0])


def _h2_energy_at_coords(tmp_path, tag, x_bohr):
    """Independent single-point H2 mrmp2 energy at exact Bohr coordinates.

    A fresh Runner (fresh Molecule/config, the full Runner.run ->
    compute_energy pipeline) is used for every point, so this is a genuinely
    different code path from the in-process loop inside pt2_numgrad.
    """
    runner = _make_runner(tmp_path, f"fd_{tag}", _H2_TEMPLATE % 0.74)
    runner.mol.update_system(np.asarray(x_bohr, dtype=float))
    runner.run(test_mod=True)
    return float(runner.mol.energies[0])


def _fit_h2_minimum(tmp_path):
    """Locate the H2 mrmp2/STO-3G minimum by a symmetric 3-point parabola fit."""
    r_c, d = 0.735, 0.015
    e_m = _h2_energy(tmp_path, "scan_m", r_c - d)
    e_c = _h2_energy(tmp_path, "scan_c", r_c)
    e_p = _h2_energy(tmp_path, "scan_p", r_c + d)
    assert e_c < e_m and e_c < e_p, "scan grid does not bracket the minimum"
    r_fit = r_c + 0.5 * d * (e_m - e_p) / (e_m - 2.0 * e_c + e_p)
    return r_fit


def test_h2_mrmp2_numgrad_zero_at_minimum_and_matches_independent_fd(tmp_path):
    """(a) H2 mrmp2: production numerical gradient ~ 0 at the 3-point-fit
    minimum, and matches an independently computed central difference of
    single-point energies (same formula, different code path) to 1e-8."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    r_fit = _fit_h2_minimum(tmp_path)
    assert 0.72 < r_fit < 0.75

    # --- production numerical gradient (runtype=grad through the dispatch) ---
    runner = _make_runner(tmp_path, "h2_grad", _H2_TEMPLATE % r_fit)
    runner.mol.config['input']['runtype'] = 'grad'
    runner.run(test_mod=True)
    grads = np.asarray(runner.mol.grads)
    assert grads.shape == (1, 2, 3)
    assert np.all(np.isfinite(grads))
    grad = grads[0].reshape(-1)

    # Near-zero at the fitted minimum.  The parabola fit on a d=0.015 A grid
    # carries an O(|E'''| d^2) ~ 2e-4 Eh/Bohr cubic residual, so 1e-3 is the
    # honest bound; transverse components vanish by symmetry.
    axial = grad[[2, 5]]
    transverse = grad[[0, 1, 3, 4]]
    assert np.max(np.abs(axial)) < 1.0e-3
    assert np.max(np.abs(transverse)) < 1.0e-6

    # --- independent central differences over the same geometry ------------
    x0 = np.asarray(runner.mol.get_system(), dtype=float).copy()
    step = 1.0e-3  # default [pt2] grad_step (Bohr)
    indep = np.zeros(6)
    for i in range(6):
        xp = x0.copy()
        xp[i] += step
        xm = x0.copy()
        xm[i] -= step
        e_p = _h2_energy_at_coords(tmp_path, f"p{i}", xp)
        e_m = _h2_energy_at_coords(tmp_path, f"m{i}", xm)
        indep[i] = (e_p - e_m) / (2.0 * step)

    np.testing.assert_allclose(grad, indep, atol=1.0e-8)


def test_h2_mrmp2_optimize_converges_to_scan_minimum(tmp_path):
    """(b) runtype=optimize (default lib=oqp) on ground-state H2 mrmp2
    reaches the scan minimum bond length within 1e-3 Angstrom."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    r_fit = _fit_h2_minimum(tmp_path)

    runner = _make_runner(tmp_path, "h2_opt", _H2_TEMPLATE % 0.760)
    runner.mol.config['input']['runtype'] = 'optimize'
    # PT2 state indices index mol.energies; single-state mrmp2 has only
    # state 0 (the [optimize] istate default of 1 addresses TDHF excited
    # states and is out of range here).
    runner.mol.config['optimize']['istate'] = 0

    t0 = time.time()
    runner.run(test_mod=True)
    wall = time.time() - t0

    x = np.asarray(runner.mol.get_system(), dtype=float)
    r_final = np.linalg.norm(x[3:] - x[:3]) * BOHR_TO_ANGSTROM
    print(f"\nPT2 optimize wall time: {wall:.2f} s; "
          f"final r = {r_final:.6f} A vs fit {r_fit:.6f} A")
    assert abs(r_final - r_fit) < 1.0e-3


def test_h4_mcqdpt2_two_state_gradient_executes(tmp_path):
    """(c) A 2-state H4 mcqdpt2 gradient run executes and returns finite,
    state-distinct gradients for the selected (excited) state."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    runner = _make_runner(tmp_path, "h4_grad", _H4_LINEAR,
                          method="mcqdpt2", nroot=2, frozen_core="1")
    runner.mol.config['input']['runtype'] = 'grad'
    runner.mol.config['properties']['grad'] = [1]  # gradient of the S1 surface
    runner.run(test_mod=True)

    energies = np.asarray(runner.mol.energies, dtype=float)
    assert energies.shape == (2,)
    assert np.all(np.isfinite(energies))
    assert energies[0] < energies[1]  # ascending multistate convention

    grads = np.asarray(runner.mol.grads)
    assert grads.shape == (2, 4, 3)
    assert np.all(np.isfinite(grads))
    # The selected state's gradient is a real array of the S1 surface and must
    # differ from S0 (the two surfaces have different slopes at this geometry).
    assert np.max(np.abs(grads[1])) > 1.0e-3
    assert np.max(np.abs(grads[1] - grads[0])) > 1.0e-3


def test_h4_mcqdpt2_meci_penalty_smoke(tmp_path):
    """(d) One/two iterations of the native penalty-function MECI driver
    (gradient-only; no NACs required) execute on H4 mcqdpt2 S0/S1."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end tests")

    runner = _make_runner(tmp_path, "h4_meci", _H4_LINEAR,
                          method="mcqdpt2", nroot=2, frozen_core="1")
    runner.mol.config['input']['runtype'] = 'meci'
    runner.mol.config['optimize']['istate'] = 0   # PT2 state indices, see above
    runner.mol.config['optimize']['jstate'] = 1
    runner.mol.config['optimize']['maxit'] = 2    # smoke: stop after 2 steps
    # meci_search default is 'penalty' (energies+gradients only)

    runner.run(test_mod=True)

    energies = np.asarray(runner.mol.energies, dtype=float)
    assert energies.shape == (2,)
    assert np.all(np.isfinite(energies))
    assert energies[1] >= energies[0]

    log_text = open(runner.mol.log, encoding="utf-8", errors="replace").read()
    assert "Geometry Optimization Step 1" in log_text
    assert "PT2 numerical gradient" in log_text
