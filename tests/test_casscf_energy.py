"""Equivalence tests for the one-call CASSCF driver (``casscf_energy``,
casscf_driver.F90).

The Python optimizer in ``oqp.library.casscf`` remains the numerical pin:
every native run must reproduce its CASSCF energies in the *same build*, from
the *same* SCF reference -- a hard-coded constant would only pin the build
against itself.  The comparison is therefore made in-process, by disabling the
one-call driver alone (``_casscf_energy_backend``) and leaving every other
native kernel in place, so both arms use identical integrals, identical CI
solves and identical Fock builds and differ only in who runs the loop.  An
out-of-process "Python arm" cannot do this: ``usercustomize`` is not imported
in a venv, so it would quietly be the native arm again.

What is compared is the ENERGY.  The macroiteration count is accumulation-order
sensitive -- mathematically equivalent kernel combinations reach an identical
energy in different numbers of steps -- so it is only checked for a generous
budget, exactly as ``test_casscf_convergers.py`` does.

Also asserts that the fixed-layout option schema in ``casscf_driver.F90``, the
copy in ``include/oqp.h`` and the mirror in ``casscf.py`` are the same list, so
the three cannot drift apart silently.
"""
import os
import re

import numpy as np
import pytest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_FORTRAN = os.path.join(_REPO, "source", "modules", "casscf_driver.F90")
_HEADER = os.path.join(_REPO, "include", "oqp.h")


def _driver_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    if not bool(getattr(oqp, "BACKEND_AVAILABLE", False)):
        return False
    from oqp import lib
    return hasattr(lib, "casscf_energy") and hasattr(oqp, "fci_ao_integrals")


requires_driver = pytest.mark.skipif(
    not _driver_available(), reason="liboqp has no casscf_energy (rebuild needed)"
)


# ------------------------------------------------------------------- schema
def _fortran_schema(prefix):
    """(name, index) pairs for the CAS_<prefix>_* parameters, in file order."""
    with open(_FORTRAN, encoding="utf-8") as handle:
        text = handle.read()
    pattern = re.compile(
        rf"^\s*integer, parameter :: CAS_{prefix}_(\w+)\s*=\s*(\d+)", re.M)
    return [(m.group(1).lower(), int(m.group(2))) for m in pattern.finditer(text)]


def _header_schema(prefix):
    with open(_HEADER, encoding="utf-8") as handle:
        text = handle.read()
    pattern = re.compile(rf"^\s*CAS_{prefix}_(\w+)\s*=\s*(\d+)", re.M)
    return [(m.group(1).lower(), int(m.group(2))) for m in pattern.finditer(text)]


def test_option_schema_matches_fortran():
    """The Python mirror of the iopt/dopt layout is the Fortran one."""
    from oqp.library.casscf import _CAS_DOPT, _CAS_IOPT

    for prefix, mirror in (("I", _CAS_IOPT), ("D", _CAS_DOPT)):
        fortran = _fortran_schema(prefix)
        assert [n for n, _ in fortran] == list(mirror), (
            f"casscf.py _CAS_{prefix}OPT does not match casscf_driver.F90")
        assert [i for _, i in fortran] == list(range(len(mirror)))


def test_option_schema_matches_header():
    """include/oqp.h documents exactly the Fortran layout."""
    from oqp.library.casscf import _CAS_DOPT, _CAS_IOPT

    assert [n for n, _ in _header_schema("I")] == list(_CAS_IOPT)
    assert [n for n, _ in _header_schema("D")] == list(_CAS_DOPT)


def test_schema_counts_are_consistent():
    from oqp.library.casscf import _CAS_DOPT, _CAS_IOPT

    text = open(_FORTRAN, encoding="utf-8").read()
    assert int(re.search(r"CAS_NIOPT\s*=\s*(\d+)", text).group(1)) == len(_CAS_IOPT)
    assert int(re.search(r"CAS_NDOPT\s*=\s*(\d+)", text).group(1)) == len(_CAS_DOPT)


def test_ci_schema_is_imported_not_redeclared():
    """The CI half reuses fci_solve's schema, so it cannot drift here.

    A copy of the FCI_I_* / FCI_D_* constants in this file would be a second
    authority for the same layout; the driver imports them instead."""
    text = open(_FORTRAN, encoding="utf-8").read()
    assert not re.search(r"^\s*integer, parameter :: FCI_[ID]_", text, re.M)
    assert "use fci_driver_mod" in text


# -------------------------------------------------------------- dispatch
@pytest.mark.parametrize("cfg,expected", [
    # An unset `converger` -- and the explicit spelling `default` -- now mean
    # `trah`, so these go through the converger framework and print a trace.
    ({}, (4, 0, True)),
    ({"converger": "twophase"}, (0, 0, False)),
    ({"converger": "default"}, (4, 0, True)),
    ({"hessian": "fd"}, (4, 0, True)),
    ({"converger": "ah"}, (1, 0, True)),
    # `trah` is its own converger now (the shared trust-region core over a
    # matrix-free Hessian-vector product), not an alias for `ah`.
    ({"converger": "trah"}, (4, 0, True)),
    ({"converger": "diis"}, (2, 0, True)),
    ({"converger": "auto"}, (3, 0, True)),
    ({"hessian": "analytic"}, (4, 1, True)),
    ({"converger": "twophase", "hessian": "analytic"}, (0, 1, False)),
    ({"converger": "ah", "hessian": "analytic"}, (1, 1, True)),
    # `2phase` resolves to the two-phase code but DOES go through the Python
    # converger framework, so it prints a trace; the flag follows the raw
    # spelling, not the resolved code.
    ({"converger": "2phase"}, (0, 0, True)),
])
def test_native_dispatch_resolves_every_supported_converger(cfg, expected):
    """The one-call driver claims all four convergers and both Hessians."""
    from oqp.library.casscf import _native_converger_codes

    assert _native_converger_codes(cfg) == expected


@pytest.mark.parametrize("cfg", [
    {"converger": "bogus"},
    {"hessian": "bogus"},
    # accepted by the input checker but rejected by `_hessian_provider`, which
    # owns that message; the driver must not quietly accept it instead
    {"hessian": "analytical"},
])
def test_native_dispatch_declines_unknown_spellings(cfg):
    """Unknown spellings go to Python, which owns the error message."""
    from oqp.library.casscf import _native_converger_codes

    assert _native_converger_codes(cfg) is None


# -------------------------------------------------------------- numerics
_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"
_H4 = "\nH 0 0 0\nH 0 0 0.740\nH 0 0 1.480\nH 0 0 2.220"
_LIH_STRETCH = "\nLi 0 0 3.0"

_CAS44 = {"active_electrons": "4", "active_orbitals": "4",
          "frozen_core": "3", "max_det": "5000"}
_CAS22 = {"active_electrons": "2", "active_orbitals": "2",
          "frozen_core": "1", "max_det": "5000"}


def _case(system, cas, *, method="casscf", basis="sto-3g", nroot="1",
          state_average=None, ci=None, casscf=None, budget=60):
    return dict(system=system, basis=basis, method=method, nroot=nroot,
                cas=cas, state_average=state_average, ci=ci or {},
                casscf=casscf or {}, budget=budget)


_CASES = {
    # the shipped shapes
    "h2o": _case(_H2O, _CAS44),
    "h4": _case(_H4, {"active_electrons": "4", "active_orbitals": "4",
                      "frozen_core": "0", "max_det": "5000"}),
    "lih": _case("\nLi 0 0 0\nH 0 0 3.0", _CAS22),
    "lih_sa": _case("\nLi 0 0 0\nH 0 0 1.6", _CAS22, method="sa-casscf",
                    nroot="2",
                    state_average={"enabled": "true", "weights": "0.5,0.5",
                                   "nstate": "2", "target_roots": "0,1",
                                   "equal_weights": "true"}),
    # option branches the examples do not reach
    # `powell` is a scaled-gradient fallback: it reaches the right CASSCF
    # solution but needs hundreds of macroiterations to do it.  Budgeted and
    # toleranced so it genuinely CONVERGES (759 iterations, |g_orb| 9.1e-06,
    # about a second on this system) rather than being cut off -- an
    # unconverged CASSCF is now an error, and comparing two runs that both
    # stopped early tested less anyway.  It lands on -7.88104526, the same
    # solution as examples/WF_methods/LiH_CASSCF.
    "powell": _case("\nLi 0 0 0\nH 0 0 1.6", _CAS22, budget=800,
                    casscf={"optimizer": "powell",
                            "max_macro_iterations": "800",
                            "gradient_norm_tol": "1.0e-5"}),
    # the two-phase path, now reachable only by asking for it: it is the one
    # converger that does NOT go through the framework, so it prints no trace
    "twophase": _case(_H2O, _CAS44, casscf={"converger": "twophase"}),
    "no_canonicalize": _case(_H2O, _CAS44, casscf={"canonicalize": "false"}),
    "excited_root": _case("\nLi 0 0 0\nH 0 0 1.6", _CAS22, nroot="2",
                          casscf={"root": "1"}),
    "davidson": _case(_H2O, _CAS44, ci={"solver": "davidson"}),
    "spin_filter": _case(_H2O, _CAS44, ci={"target_spin": "singlet"}),
    # the convergers and the analytic Hessian, each against its NumPy pin
    "analytic": _case(_H2O, _CAS44, casscf={"hessian": "analytic"}),
    "ah": _case(_H2O, _CAS44, casscf={"converger": "ah"}),
    "ah_analytic": _case(_H2O, _CAS44,
                         casscf={"converger": "ah", "hessian": "analytic"}),
    "diis": _case(_H2O, _CAS44, casscf={"converger": "diis"}),
    "auto": _case(_H2O, _CAS44, casscf={"converger": "auto"}),
    "ah_lih": _case("\nLi 0 0 0\nH 0 0 3.0", _CAS22, casscf={"converger": "ah"}),
    "diis_analytic_lih": _case("\nLi 0 0 0\nH 0 0 3.0", _CAS22,
                               casscf={"converger": "diis",
                                       "hessian": "analytic"}),
    "auto_analytic_lih": _case("\nLi 0 0 0\nH 0 0 3.0", _CAS22,
                               casscf={"converger": "auto",
                                       "hessian": "analytic"}),
    "sa_ah_analytic": _case("\nLi 0 0 0\nH 0 0 1.6", _CAS22,
                            method="sa-casscf", nroot="2",
                            state_average={"enabled": "true",
                                           "weights": "0.5,0.5",
                                           "nstate": "2",
                                           "target_roots": "0,1",
                                           "equal_weights": "true"},
                            casscf={"converger": "ah",
                                    "hessian": "analytic"}),
}

# How closely the native and the Python arm must land on the same energy.
#
# Both stop when the orbital gradient falls below `gradient_norm_tol` (1e-6 by
# default), so they stop at slightly DIFFERENT points inside that ball.  Near a
# stationary point the energy is quadratic, so the two energies differ by about
# |g|^2 / lambda -- with the smallest orbital curvatures of these active spaces
# (lambda ~ 5e-2) that is a few times 1e-8, and lih_sa measures 2e-8.
#
# The old 1e-9 was not a principled figure: it held because the two-phase
# converger, which used to be the default, happened to stop at very nearly the
# same point on both arms.  It is not the right yardstick for a converger that
# stops on a gradient norm.  1e-7 is the quadratic estimate with room, and is
# still two orders below chemical accuracy.
_AGREEMENT_TOL = 1.0e-7


def _run(workdir, project, case):
    from oqp.pyoqp import Runner

    ci = {"nroot": case["nroot"], "solver": "dense", "eig_tol": "1.0e-10",
          "integral_backend": "native", "target_spin": "any",
          "root_tracking": "energy"}
    ci.update(case.get("ci", {}))
    casscf = {"max_macro_iterations": "50", "optimizer": "newton",
              "root": "0", "gradient_norm_tol": "1.0e-7",
              "energy_decrease_tol": "1.0e-10", "step_norm_tol": "1.0e-8",
              "max_rotation_norm": "2.0e-1", "level_shift": "1.0e-3",
              "canonicalize": "true"}
    casscf.update(case.get("casscf", {}))
    input_dict = {
        "input": {"system": case["system"], "charge": "0",
                  "basis": case["basis"], "method": case["method"],
                  "runtype": "energy"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "forced_attempt": "3", "save_molden": "False",
                "conv": "1.0e-10"},
        "cas": case["cas"],
        "ci": ci,
        "casscf": casscf,
        "tests": {"exception": "True"},
    }
    if case["state_average"] is not None:
        input_dict["state_average"] = case["state_average"]

    runner = Runner(project=project, input_file=None,
                    log=str(workdir / f"{project}.log"),
                    input_dict=input_dict, silent=1, usempi=False)
    runner.run(test_mod=True)
    text = open(runner.mol.log, encoding="utf-8").read()
    return (np.asarray(runner.mol.data["OQP::CASSCF_ENERGIES"], dtype=float),
            int(re.search(r"CASSCF macro iterations:\s+(\d+)\s*/", text).group(1)),
            re.search(r"CASSCF converged:\s+yes", text) is not None,
            text)


@requires_driver
@pytest.mark.parametrize("name", sorted(_CASES))
def test_native_run_reproduces_the_python_optimizer(tmp_path, monkeypatch, name):
    """Same build, same process: only the one-call driver is switched off."""
    import oqp.library.casscf as casscf

    case = _CASES[name]
    e_native, it_native, conv_native, _ = _run(tmp_path, f"{name}_native", case)

    monkeypatch.setattr(casscf, "_casscf_energy_backend", lambda: None)
    e_python, it_python, conv_python, _ = _run(tmp_path, f"{name}_python", case)

    assert conv_native == conv_python
    assert e_native.shape == e_python.shape
    assert np.max(np.abs(e_native - e_python)) < _AGREEMENT_TOL, (
        f"{name}: native {e_native} vs python {e_python}")
    # Iteration counts are accumulation-order sensitive; only budgeted.
    assert it_native <= case["budget"]
    assert it_python <= case["budget"]


@requires_driver
@pytest.mark.parametrize("name", ["h2o", "ah_analytic", "diis", "auto"])
def test_native_path_is_actually_taken(tmp_path, name):
    """One crossing per run, not one per microiteration.

    The whole point of the driver is that the boundary is crossed once; a
    silent fallback would still produce the right energy, so the crossing
    count is what proves the native path ran.  ``ah_analytic`` is the
    configuration that used to cross it hundreds of times: every AH model
    step, every Hessian column build, every CI-relaxation amplitude and every
    RDM was its own call."""
    import collections

    import oqp.library.fci as fci

    counter = collections.Counter()
    lib, ffi = fci._lib_backend()

    class _Counting:
        def __getattr__(self, attr):
            fn = getattr(lib, attr)
            if not callable(fn):
                return fn

            def wrapped(*args, **kwargs):
                counter[attr] += 1
                return fn(*args, **kwargs)
            return wrapped

    saved = fci._BACKEND_CACHE
    fci._BACKEND_CACHE = (_Counting(), ffi)
    try:
        _run(tmp_path, f"{name}_crossings", _CASES[name])
    finally:
        fci._BACKEND_CACHE = saved

    assert counter["casscf_energy"] == 1
    # every per-microiteration kernel is now called from inside liboqp
    for kernel in ("fci_solve", "casscf_gfock_grad", "casscf_orbital_rotate",
                   "mo_transform_eri", "rdm2_spatial", "casscf_ah_model_step",
                   "casscf_diis_coeffs", "casscf_hess_bmat", "casscf_hess_amp",
                   "casscf_hess_wmat", "casscf_hess_relax",
                   "casscf_excitation_stack", "fci_dense_hamiltonian"):
        assert counter[kernel] == 0, f"{kernel} still crosses the boundary"


@requires_driver
def test_converger_trace_survives_the_round_trip(tmp_path):
    """The driver returns counters; Python formats the same trace it used to.

    The trace is the only user-visible evidence that a non-default converger
    ran, so its lines have to come out of the native path unchanged."""
    _e, _n, _c, text = _run(tmp_path, "h2o_trace", _CASES["ah_analytic"])
    assert re.search(r"PyOQP converger:\s+ah", text)
    assert re.search(r"PyOQP orbital hessian:\s+analytic", text)
    assert re.search(r"PyOQP converger CI evaluations:\s+\d+", text)
    builds = int(re.search(r"PyOQP analytic hessian builds:\s+(\d+)",
                           text).group(1))
    assert builds >= 1

    # ... and the two-phase path, which does not go through the framework,
    # still carries no trace at all.  (It is no longer the default: an unset
    # `converger` means `trah`, which does print one.)
    _e, _n, _c, plain = _run(tmp_path, "h2o_notrace", _CASES["twophase"])
    assert "PyOQP converger:" not in plain


@requires_driver
def test_macroiteration_table_survives_the_round_trip(tmp_path):
    """The driver hands the history back and Python prints it unchanged.

    The log is the one thing that must not change shape when the loop moves
    into Fortran, so the table is parsed back out and checked against the
    summary lines it has to agree with."""
    _e, niter, converged, text = _run(tmp_path, "h2o_table", _CASES["h2o"])
    assert converged

    block = text.split("--- macro iterations ---", 1)[1]
    rows = []
    for line in block.splitlines():
        m = re.match(r"\s+(\d+)\s+(-?\d+\.\d+)\s+(\S*)\s*(\S+)\s*(\S*)\s*$", line)
        if m is None:
            if rows:
                break
            continue
        rows.append((int(m.group(1)), float(m.group(2)), float(m.group(4))))

    assert len(rows) >= 2
    assert rows[0][0] == 0                      # the seed row carries no dE/|step|
    # The reported count is the macroiteration the loop STOPPED in; when it
    # stops on the gradient test at the top of that iteration no row is
    # appended for it, so the table ends one short.  Both are the Python
    # behaviour and the driver reproduces it.
    assert rows[-1][0] in (niter, niter - 1)
    final_g = float(re.search(r"CASSCF final \|g_orb\|:\s+(\S+)", text).group(1))
    assert rows[-1][2] == pytest.approx(final_g, rel=1e-3)
    # a converged run ends at (or below) its own starting energy
    assert rows[-1][1] <= rows[0][1] + 1.0e-12
