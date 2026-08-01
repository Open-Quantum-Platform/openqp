"""Correctness gate for the analytic MCSCF orbital-rotation Hessian
(``pyoqp/oqp/library/casscf_hessian.py``) and its ``[casscf] hessian``
wiring into the converger framework.

Contracts:

  (1) 2-RDM source validation: the full-space spin-summed 1-/2-RDMs built
      from the determinant CI vectors reproduce the CASCI energy through
      ``E = sum h*D1 + 0.5 sum (pq|rs)*D2 + Enuc`` to 1e-10, per root.
  (2) Machinery anchors: the Z directional-derivative intermediate
      reproduces the production generalized-Fock orbital gradient, and the
      core-folded active Hamiltonian built from the excitation-matrix stack
      reproduces the CI solver's energies.
  (3) The analytic Hessian matches ``casscf._fd_orbital_hessian`` to
      atol 1e-6 on stretched LiH CAS(2,2), H2O CAS(4,4) and
      SA2-CASSCF(2,2)/H4 -- at BOTH the initial (RHF) orbitals and the
      converged CASSCF orbitals.  The residual is FD truncation
      (O(h^2), h = 1e-4); the synthetic h^2-scaling check during development
      confirmed exactness (error ratio 9.00 for step ratio 3).
  (4) Converger wiring: 'ah' with ``hessian = analytic`` reproduces the
      fd-'ah' energies to <= 1e-8 Eh with a large drop in CI evaluations
      (trace counter), the default path stays untouched when the key is
      absent (or explicitly 'fd'), and an unknown value raises.

The ``hessian`` key is read with a ``dict.get`` default and has no schema
entry yet (proposed rows are in WF_CONVERGER_NOTES.md), so tests inject it
into ``runner.mol.config`` after validation -- the same supported
programmatic path the converger tests use.
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

_SYSTEMS = {
    "lih": dict(system=_LIH_STRETCH, basis="sto-3g", cas=_CAS22),
    "h2o": dict(system=_H2O, basis="sto-3g", cas=_CAS44),
    "h4_sa": dict(system=_H4, basis="sto-3g", cas=_CAS22, method="sa-casscf",
                  nroot="2", state_average={"enabled": "true", "nstate": "2"}),
}

_HESS_TOL = 1.0e-6          # analytic-vs-FD gate (FD truncation is ~1e-7 here)
_ENERGY_TOL = 1.0e-8        # cross-run energy agreement
_IDENTITY_TOL = 1.0e-10     # 2-RDM energy identity


def _run(workdir, project, *, system, basis, cas=None, method="casscf",
         nroot="1", state_average=None, inject=None):
    from oqp.pyoqp import Runner

    input_dict = {
        "input": {"system": system, "charge": "0", "basis": basis,
                  "method": method, "runtype": "energy"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "forced_attempt": "3", "save_molden": "False"},
        "tests": {"exception": "True"},
    }
    if cas is not None:
        input_dict["cas"] = cas
        input_dict["ci"] = {"nroot": nroot, "solver": "dense", "eig_tol": "1.0e-10",
                            "integral_backend": "native", "target_spin": "any",
                            "root_tracking": "energy"}
        input_dict["casscf"] = {"max_macro_iterations": "50", "optimizer": "newton",
                                "root": "0", "gradient_norm_tol": "1.0e-7",
                                "energy_decrease_tol": "1.0e-10",
                                "step_norm_tol": "1.0e-8",
                                "max_rotation_norm": "2.0e-1",
                                "level_shift": "1.0e-3", "canonicalize": "true"}
    if state_average is not None:
        input_dict["state_average"] = state_average

    runner = Runner(project=project, input_file=None,
                    log=str(workdir / f"{project}.log"),
                    input_dict=input_dict, silent=1, usempi=False)
    if inject:
        # runtime injection through the dict-get path (see module docstring)
        runner.mol.config["casscf"].update(inject)
    runner.run(test_mod=True)
    return runner


def _context(runner_cas, runner_hf):
    """Rebuild the ``casscf._optimize`` evaluation context for a finished run.

    Mirrors the header of ``CASSCF.energy`` (same buffers, same closures) so
    the FD and the analytic Hessian are probed on exactly the production
    objective.  ``C_conv`` is the committed (canonicalized) CASSCF solution;
    ``C_init`` is the RHF reference from an identical plain-HF run -- the
    orbitals every CASSCF macroiteration starts from."""
    import oqp
    from oqp.library.casscf import (
        CASSCF, _casscf_options, _generalized_fock, _nonredundant_pairs,
        _solve_active,
    )
    from oqp.library.fci import (
        settings_from_casci_config, _transform_integrals, _unpack_lower_triangle,
    )

    mol = runner_cas.mol
    settings = settings_from_casci_config(mol.config)
    options = _casscf_options(mol.config)
    nbf = int(mol.data.get_basis()["nbf"])
    oqp.fci_ao_integrals(mol)
    hcore_ao = _unpack_lower_triangle(np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
    eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape((nbf,) * 4, order="F")
    enuc = float(mol.mol_energy.nenergy)
    ncore = int(settings.frozen_core)
    nact = int(settings.active_orbitals)
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    active_nelec = (nelec[0] - ncore, nelec[1] - ncore)
    weights, roots, _sa = CASSCF(mol)._state_average_plan(settings, options)
    pairs = _nonredundant_pairs(ncore, nact, nbf)

    C_conv = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
    molh = runner_hf.mol
    oqp.fci_ao_integrals(molh)
    C_init = np.asarray(molh.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T

    def evaluate(C):
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, C)
        energies, coeffs, dets, D, G = _solve_active(
            h1e, eri, ncore, nact, active_nelec, enuc, settings, weights, roots)
        objective = float(np.dot(weights, energies[roots]))
        F = _generalized_fock(D, G, h1e, eri)
        grad = np.array([2.0 * (F[q, p] - F[p, q]) for (p, q) in pairs])
        return objective, grad, energies, coeffs

    def transform(C):
        return _transform_integrals(hcore_ao, eri_ao, C)

    return dict(mol=mol, settings=settings, nbf=nbf, enuc=enuc, ncore=ncore,
                nact=nact, active_nelec=active_nelec, weights=weights,
                roots=roots, pairs=pairs, C_init=C_init, C_conv=C_conv,
                evaluate=evaluate, transform=transform)


@pytest.fixture(scope="module")
def contexts(tmp_path_factory):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    workdir = tmp_path_factory.mktemp("casscf_hessian")
    cache = {}
    for key, spec in _SYSTEMS.items():
        spec = dict(spec)
        hf = _run(workdir, f"{key}_hf", system=spec["system"], basis=spec["basis"],
                  method="hf")
        cas = _run(workdir, f"{key}_cas", **spec)
        cache[key] = _context(cas, hf)
    return cache


# ------------------------------------------------------------------ (1) 2-RDM identity
@pytest.mark.parametrize("case", list(_SYSTEMS))
def test_2rdm_energy_identity(contexts, case):
    """E = sum h*D1 + 0.5 sum (pq|rs)*D2 + Enuc reproduces every CASCI root
    energy to 1e-10 (full-space spin-summed RDMs from the CI vectors)."""
    from oqp.library.casscf import _full_rdms
    from oqp.library.fci import _determinants
    from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

    ctx = contexts[case]
    h1e, eri = ctx["transform"](ctx["C_conv"])
    _obj, _grad, energies, coeffs = ctx["evaluate"](ctx["C_conv"])
    dets = _determinants(ctx["nact"], ctx["active_nelec"])
    for r in ctx["roots"]:
        gamma = make_rdm1_spatial(coeffs[:, r], dets, ctx["nact"])
        Gamma = make_rdm2_spatial(coeffs[:, r], dets, ctx["nact"])
        D1, D2 = _full_rdms(gamma, Gamma, ctx["ncore"], ctx["nact"], ctx["nbf"])
        e_rdm = float(np.einsum("pq,pq", h1e, D1)
                      + 0.5 * np.einsum("pqrs,pqrs", eri, D2)
                      + ctx["enuc"])
        assert e_rdm == pytest.approx(float(energies[r]), abs=_IDENTITY_TOL)


# ------------------------------------------------------------------ (2) machinery anchors
@pytest.mark.parametrize("case", list(_SYSTEMS))
def test_z_intermediate_reproduces_gradient(contexts, case):
    """The Z directional-derivative intermediate evaluated on the plain
    integrals reproduces the production gradient 2(F_qp - F_pq)."""
    from oqp.library.casscf import _full_rdms
    from oqp.library.casscf_hessian import _pair_differences, _z_matrix
    from oqp.library.fci import _determinants
    from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

    ctx = contexts[case]
    for C in (ctx["C_init"], ctx["C_conv"]):
        h1e, eri = ctx["transform"](C)
        _obj, grad, _energies, coeffs = ctx["evaluate"](C)
        dets = _determinants(ctx["nact"], ctx["active_nelec"])
        D = np.zeros((ctx["nbf"], ctx["nbf"]))
        G = np.zeros((ctx["nbf"],) * 4)
        for w, r in zip(ctx["weights"], ctx["roots"]):
            gamma = make_rdm1_spatial(coeffs[:, r], dets, ctx["nact"])
            Gamma = make_rdm2_spatial(coeffs[:, r], dets, ctx["nact"])
            Dr, Gr = _full_rdms(gamma, Gamma, ctx["ncore"], ctx["nact"], ctx["nbf"])
            D += w * Dr
            G += w * Gr
        gz = _pair_differences(_z_matrix(D, G, h1e, eri), ctx["pairs"])
        assert np.max(np.abs(gz - grad)) < 1.0e-10


@pytest.mark.parametrize("case", list(_SYSTEMS))
def test_active_hamiltonian_reproduces_ci_energies(contexts, case):
    """Core folding + excitation-matrix operator assembly reproduce the CI
    solver: folded integrals match fci._active_space exactly and the dense
    active Hamiltonian's lowest eigenvalues match the solved root energies."""
    from oqp.library.casscf_hessian import (
        _active_operator_matrix, _excitation_matrices, _fold_active,
    )
    from oqp.library.fci import _active_space, _symmetric_eigh

    ctx = contexts[case]
    h1e, eri = ctx["transform"](ctx["C_conv"])
    _obj, _grad, energies, _coeffs = ctx["evaluate"](ctx["C_conv"])
    nel_tot = (ctx["ncore"] + ctx["active_nelec"][0], ctx["ncore"] + ctx["active_nelec"][1])
    h_ref, g_ref, _an, ecore, _meta = _active_space(
        h1e, eri, nel_tot, ctx["enuc"], ctx["settings"])
    f0, g0 = _fold_active(h1e, eri, ctx["ncore"], ctx["nact"])
    assert np.max(np.abs(f0 - h_ref)) < 1.0e-12
    assert np.max(np.abs(g0 - g_ref)) < 1.0e-12

    _dets, stack = _excitation_matrices(ctx["nact"], *ctx["active_nelec"])
    eps, _vecs = _symmetric_eigh(_active_operator_matrix(f0, g0, stack))
    nroot = max(ctx["roots"]) + 1
    for r in range(nroot):
        assert float(eps[r] + ecore) == pytest.approx(float(energies[r]), abs=1.0e-9)


# ------------------------------------------------------------------ (3) FD agreement gate
@pytest.mark.parametrize("point", ["initial", "converged"])
@pytest.mark.parametrize("case", list(_SYSTEMS))
def test_analytic_hessian_matches_fd(contexts, case, point):
    """The arbiter: analytic vs symmetrized-FD orbital Hessian, atol 1e-6."""
    from oqp.library.casscf import _fd_orbital_hessian
    from oqp.library.casscf_hessian import analytic_orbital_hessian

    ctx = contexts[case]
    C = ctx["C_init"] if point == "initial" else ctx["C_conv"]
    _obj, grad, _energies, coeffs = ctx["evaluate"](C)
    h1e, eri = ctx["transform"](C)
    H_fd = _fd_orbital_hessian(C, ctx["evaluate"], ctx["pairs"], ctx["nbf"])
    H_an = analytic_orbital_hessian(
        h1e, eri, ctx["ncore"], ctx["nact"], ctx["active_nelec"], ctx["pairs"],
        ctx["weights"], ctx["roots"], coeffs)
    diff = float(np.max(np.abs(H_an - H_fd)))
    print(f"\n{case}/{point}: n_par={len(ctx['pairs'])} |g|={np.linalg.norm(grad):.2e} "
          f"|H_fd|max={np.max(np.abs(H_fd)):.3f} max|analytic-fd|={diff:.3e}")
    assert H_an == pytest.approx(H_an.T), "analytic Hessian must be symmetric"
    assert diff < _HESS_TOL


# ------------------------------------------------------------------ (4) converger wiring
def _energy(runner) -> float:
    return float(np.asarray(runner.mol.data["OQP::CASSCF_ENERGIES"], dtype=float)[0])


def _log_text(runner) -> str:
    with open(runner.mol.log) as handle:
        return handle.read()


def _ci_evaluations(text: str) -> int:
    return int(re.search(r"converger CI evaluations:\s+(\d+)", text).group(1))


def test_ah_analytic_matches_fd_ah_with_fewer_ci_solves(tmp_path):
    """'ah' + hessian=analytic reproduces fd-'ah' energies (<= 1e-8 Eh) while
    the CI-evaluation counter collapses (no 2*n_par FD solves per iteration)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    rows = []
    for key in ("lih", "h2o"):
        spec = dict(_SYSTEMS[key])
        fd = _run(tmp_path, f"{key}_ah_fd", inject={"converger": "ah"}, **spec)
        an = _run(tmp_path, f"{key}_ah_analytic",
                  inject={"converger": "ah", "hessian": "analytic"}, **spec)
        fd_text, an_text = _log_text(fd), _log_text(an)
        assert re.search(r"CASSCF converged:\s+yes", fd_text)
        assert re.search(r"CASSCF converged:\s+yes", an_text)
        assert re.search(r"orbital hessian:\s+analytic", an_text)
        assert re.search(r"analytic hessian builds:\s+\d+", an_text)
        assert "orbital hessian:" not in fd_text
        assert _energy(an) == pytest.approx(_energy(fd), abs=_ENERGY_TOL)
        n_fd, n_an = _ci_evaluations(fd_text), _ci_evaluations(an_text)
        n_builds = int(re.search(r"analytic hessian builds:\s+(\d+)", an_text).group(1))
        rows.append((key, n_fd, n_an, n_builds, _energy(fd), _energy(an)))
        # "far fewer": each macroiteration drops from ~2*n_par+1 to ~1 CI solve
        assert n_an * 3 < n_fd, f"{key}: analytic-ah used {n_an} vs fd-ah {n_fd}"

    print("\n'ah' CI evaluations, fd vs analytic hessian:")
    print(f"  {'case':8s} {'fd':>6s} {'analytic':>9s} {'hess builds':>12s}  E(fd) == E(analytic)")
    for key, n_fd, n_an, n_builds, e_fd, e_an in rows:
        print(f"  {key:8s} {n_fd:6d} {n_an:9d} {n_builds:12d}  {e_fd:.10f} == {e_an:.10f}")


def test_sa_ah_analytic_matches_fd_ah(tmp_path):
    """SA2-CASSCF: the analytic Hessian carries the state-average weights
    (weighted RDMs + weighted CI-relaxation term)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    spec = dict(_SYSTEMS["h4_sa"])
    fd = _run(tmp_path, "h4_sa_ah_fd", inject={"converger": "ah"}, **spec)
    an = _run(tmp_path, "h4_sa_ah_analytic",
              inject={"converger": "ah", "hessian": "analytic"}, **spec)
    assert re.search(r"CASSCF converged:\s+yes", _log_text(an))
    assert _energy(an) == pytest.approx(_energy(fd), abs=_ENERGY_TOL)
    assert _ci_evaluations(_log_text(an)) * 3 < _ci_evaluations(_log_text(fd))


def test_default_newton_path_accepts_analytic_hessian(tmp_path):
    """No converger key (two-phase default) + hessian=analytic: the default
    Newton loop consumes the analytic Hessian and lands on the same solution."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    spec = dict(_SYSTEMS["lih"])
    ref = _run(tmp_path, "lih_default", **spec)
    ana = _run(tmp_path, "lih_default_analytic", inject={"hessian": "analytic"}, **spec)
    assert re.search(r"CASSCF converged:\s+yes", _log_text(ana))
    assert _energy(ana) == pytest.approx(_energy(ref), abs=_ENERGY_TOL)
    # the default two-phase log carries no converger/hessian trace either way
    assert "PyOQP converger:" not in _log_text(ana)


def test_hessian_fd_explicit_is_default_path(tmp_path):
    """hessian=fd (explicit) is the identical FD path: same energy as the key
    being absent and no hessian trace in the log."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    spec = dict(_SYSTEMS["lih"])
    ref = _run(tmp_path, "lih_ah_ref", inject={"converger": "ah"}, **spec)
    exp = _run(tmp_path, "lih_ah_fd_explicit",
               inject={"converger": "ah", "hessian": "fd"}, **spec)
    ref_text, exp_text = _log_text(ref), _log_text(exp)
    assert _energy(exp) == pytest.approx(_energy(ref), abs=1.0e-12)
    assert _ci_evaluations(exp_text) == _ci_evaluations(ref_text)
    assert "orbital hessian:" not in exp_text


def test_unknown_hessian_value_raises(tmp_path):
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run CASSCF Hessian tests")
    with pytest.raises(ValueError, match="hessian"):
        _run(tmp_path, "lih_bad_hessian", inject={"hessian": "bogus"},
             **dict(_SYSTEMS["lih"]))
