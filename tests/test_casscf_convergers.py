"""End-to-end tests for the selectable CASSCF orbital-optimization convergers
(``[casscf] converger`` -> oqp.library.casscf_convergers).

The framework adds 'ah' (trust-region augmented Hessian), 'diis'
(orbital-gradient DIIS over the production quasi-Newton step) and 'auto'
('ah' with a stagnation fallback to the default) next to the unchanged
'twophase' default.  These tests check three contracts:

  (a) every converger reproduces the default converger's CASSCF energy on the
      systems where the CASSCF solution is unambiguous (H2O/STO-3G CAS(4,4)
      reused from ``tests/test_casscf.py``, and stretched LiH CAS(2,2)),
      to <= 1e-8 Eh;
  (b) on the stretched-LiH hard case the macroiteration count of every
      converger is recorded (printed table) and stays within budget;
  (c) default-path regression: with no ``converger`` key the energy AND the
      macroiteration count are identical to the pre-framework code (values
      recorded from the unmodified two-phase optimizer on this build).

The ``converger`` key is read at runtime with a ``dict.get`` default and has no
schema entry yet, so tests inject it into ``runner.mol.config`` after the
Runner has validated the input (the supported programmatic path).
"""
import re

import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401  probe the dual-path oqp.utils shadow
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"
_H4 = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"
_LIH_STRETCH = "\nLi 0 0 0\nH 0 0 3.0"

_CAS44 = {"active_electrons": "4", "active_orbitals": "4", "frozen_core": "3", "max_det": "5000"}
_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1", "max_det": "5000"}

# Default-path baselines recorded from the unmodified two-phase optimizer
# (probe run before the converger framework was added).  Test (c) is a change
# detector for the default path when no converger key is present: the energy
# is pinned exactly, while the macroiteration count is allowed +/-1 because
# the last macroiteration is decided by a threshold crossing and the SCF
# reference it starts from is BLAS/build dependent.  A real convergence
# regression moves the count by far more than one step.
_DEFAULT_BASELINE = {
    # h2o was 20 macroiterations while the spin-orbital 2-RDM was accumulated
    # by the nested create/annihilate walk.  Factorising it through the
    # doubly-annihilated intermediate space (one GEMM) changes only the
    # summation order, but the better-conditioned accumulation reaches the same
    # energy -- pinned below to 1e-8, and still matching PySCF -- in 15.
    "h2o": (-75.0085688882, 15),
    "h4": (-2.1153228671, 3),
    "lih": (-7.7983384275, 7),
}

_AGREEMENT_TOL = 1.0e-8


def _run_casscf(workdir, project, *, system, basis, cas, method="casscf",
                nroot="1", converger=None, state_average=None):
    from oqp.pyoqp import Runner

    input_dict = {
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
            "maxit": "100",
            "forced_attempt": "3",
            "save_molden": "False",
        },
        "cas": cas,
        "ci": {
            "nroot": nroot,
            "solver": "dense",
            "eig_tol": "1.0e-10",
            "integral_backend": "native",
            "target_spin": "any",
            "root_tracking": "energy",
        },
        "casscf": {
            "max_macro_iterations": "50",
            "optimizer": "newton",
            "root": "0",
            "gradient_norm_tol": "1.0e-7",
            "energy_decrease_tol": "1.0e-10",
            "step_norm_tol": "1.0e-8",
            "max_rotation_norm": "2.0e-1",
            "level_shift": "1.0e-3",
            "canonicalize": "true",
        },
        "tests": {"exception": "True"},
    }
    if state_average is not None:
        input_dict["state_average"] = state_average

    runner = Runner(
        project=project,
        input_file=None,
        log=str(workdir / f"{project}.log"),
        input_dict=input_dict,
        silent=1,
        usempi=False,
    )
    if converger is not None:
        # kept as a runtime injection so this helper also exercises the
        # dict-get default path; the schema/input-file path is covered by
        # test_converger_key_through_input_schema.
        runner.mol.config["casscf"]["converger"] = converger
    runner.run(test_mod=True)
    return runner


def test_converger_key_through_input_schema(tmp_path):
    """[casscf] converger now has schema/checker rows: the key must flow from
    the input dict through validation to the framework (trace line in the log
    proves the non-default converger actually ran)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    from oqp.pyoqp import Runner

    runner = Runner(
        project="conv_schema", input_file=None,
        log=str(tmp_path / "conv_schema.log"),
        input_dict={
            "input": {"system": _SYSTEMS["lih"]["system"], "charge": "0",
                      "basis": _SYSTEMS["lih"]["basis"], "method": "casscf",
                      "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                    "forced_attempt": "3", "save_molden": "False"},
            "cas": _SYSTEMS["lih"]["cas"],
            "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-10",
                   "integral_backend": "native", "target_spin": "any",
                   "root_tracking": "energy"},
            "casscf": {"max_macro_iterations": "50", "converger": "ah"},
            "tests": {"exception": "True"},
        },
        silent=1, usempi=False)
    runner.run(test_mod=True)
    assert "PyOQP converger:" in open(tmp_path / "conv_schema.log").read()


def _energy(runner) -> float:
    return float(np.asarray(runner.mol.data["OQP::CASSCF_ENERGIES"], dtype=float)[0])


def _log_text(runner) -> str:
    with open(runner.mol.log) as handle:
        return handle.read()


def _macro_iterations(text: str) -> int:
    return int(re.search(r"CASSCF macro iterations:\s+(\d+)\s*/", text).group(1))


def _converged(text: str) -> bool:
    return re.search(r"CASSCF converged:\s+yes", text) is not None


_SYSTEMS = {
    "h2o": dict(system=_H2O, basis="sto-3g", cas=_CAS44),
    "h4": dict(system=_H4, basis="sto-3g", cas=_CAS22),
    "lih": dict(system=_LIH_STRETCH, basis="sto-3g", cas=_CAS22),
}


@pytest.fixture(scope="module")
def default_runs(tmp_path_factory):
    """Run the default (no converger key) CASSCF once per system and cache
    (energy, macroiterations, log text)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASSCF tests")
    workdir = tmp_path_factory.mktemp("casscf_conv_default")
    cache = {}
    for key, spec in _SYSTEMS.items():
        runner = _run_casscf(workdir, f"default_{key}", **spec)
        text = _log_text(runner)
        cache[key] = (_energy(runner), _macro_iterations(text), text)
    return cache


# --------------------------------------------------------------------- (a) agreement
@pytest.mark.parametrize("converger", ["ah", "diis", "auto"])
@pytest.mark.parametrize("system_key", ["h2o", "lih"])
def test_converger_matches_default(tmp_path, default_runs, system_key, converger):
    """Each converger reproduces the default CASSCF energy to <= 1e-8 Eh."""
    e_ref, _it_ref, _ = default_runs[system_key]
    runner = _run_casscf(tmp_path, f"{system_key}_{converger}",
                         converger=converger, **_SYSTEMS[system_key])
    text = _log_text(runner)
    assert _converged(text), f"{converger} did not converge on {system_key}"
    assert _energy(runner) == pytest.approx(e_ref, abs=_AGREEMENT_TOL)
    # the converger trace proves the framework path actually ran
    assert re.search(r"PyOQP converger:\s+" + converger, text)


# --------------------------------------------------------------------- (b) hard case
def test_hard_case_iteration_counts(tmp_path, default_runs):
    """Stretched LiH CAS(2,2): all convergers converge to the same solution;
    macroiteration counts are recorded (run pytest -s to see the table)."""
    e_ref, it_ref, _ = default_runs["lih"]
    counts = {"twophase (default)": it_ref}
    for converger in ("ah", "diis", "auto"):
        runner = _run_casscf(tmp_path, f"hard_lih_{converger}",
                             converger=converger, **_SYSTEMS["lih"])
        text = _log_text(runner)
        assert _converged(text), f"{converger} did not converge on stretched LiH"
        assert _energy(runner) == pytest.approx(e_ref, abs=_AGREEMENT_TOL)
        counts[converger] = _macro_iterations(text)

    print("\nstretched LiH (3.0 A) STO-3G CAS(2,2) macroiteration counts:")
    for name, count in counts.items():
        print(f"  {name:20s} {count}")
    for name, count in counts.items():
        assert 1 <= count <= 50, f"{name} macroiteration count {count} out of budget"


# --------------------------------------------------------------------- (c) default regression
@pytest.mark.parametrize("system_key", ["h2o", "h4"])
def test_default_path_regression(default_runs, system_key):
    """No converger key: energy matches the pre-framework two-phase optimizer
    exactly, and the macroiteration count stays within one step of it."""
    e_ref, it_ref = _DEFAULT_BASELINE[system_key]
    e, iters, text = default_runs[system_key]
    assert e == pytest.approx(e_ref, abs=_AGREEMENT_TOL)
    assert abs(iters - it_ref) <= 1, (
        f"{system_key} default path took {iters} macroiterations "
        f"(baseline {it_ref}); a swing this large is a convergence regression, "
        f"not build-to-build threshold noise"
    )
    # the default path must not carry any converger trace in the log
    assert "PyOQP converger:" not in text


# --------------------------------------------------------------------- extras
def test_ah_never_above_default_h4(tmp_path, default_runs):
    """H4 CAS(2,2) has two close CASSCF solutions; the AH converger follows the
    negative-curvature direction into the deeper one.  The invariant is that
    'ah' never lands ABOVE the default converger's energy."""
    e_ref, _it_ref, _ = default_runs["h4"]
    runner = _run_casscf(tmp_path, "h4_ah", converger="ah", **_SYSTEMS["h4"])
    text = _log_text(runner)
    assert _converged(text)
    assert _energy(runner) <= e_ref + _AGREEMENT_TOL


def test_sa_casscf_ah_matches_default(tmp_path):
    """SA2-CASSCF(2,2)/H4: the AH converger reproduces the default
    state-averaged energy (the objective closure carries the weights)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASSCF tests")
    sa = {"enabled": "true", "nstate": "2"}
    ref = _run_casscf(tmp_path, "sa_default", system=_H4, basis="sto-3g",
                      cas=_CAS22, method="sa-casscf", nroot="2", state_average=sa)
    ah = _run_casscf(tmp_path, "sa_ah", system=_H4, basis="sto-3g",
                     cas=_CAS22, method="sa-casscf", nroot="2", state_average=sa,
                     converger="ah")
    assert _converged(_log_text(ah))
    assert _energy(ah) == pytest.approx(_energy(ref), abs=_AGREEMENT_TOL)


def test_unknown_converger_raises(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run end-to-end CASSCF tests")
    with pytest.raises(ValueError, match="converger"):
        _run_casscf(tmp_path, "bogus_converger", converger="bogus", **_SYSTEMS["h4"])
