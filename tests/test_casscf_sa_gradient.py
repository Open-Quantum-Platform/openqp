"""SA-CASSCF nuclear gradients: the weighted objective, and individual roots.

Two derivatives are under test and they are not variants of one another
(``docs/sa_casscf_gradients.md``):

* the WEIGHTED OBJECTIVE ``L = sum_I w_I E_I`` is stationary in every
  wavefunction parameter, so its gradient carries no response term;
* an INDIVIDUAL averaged root is stationary against CI variation and *not*
  against orbital rotation, so it needs a coupled orbital+CI Z-vector.

Layered on purpose, cheapest and sharpest first:

1. **The abstract protocol.**  A basis-free SA-CASSCF model whose "geometry" is
   a scalar parameter inside the integrals.  There is no AO basis, no
   derivative integral and no Pulay term, so

       dE_J/dx = sum D^eff (dh/dx) + 1/2 sum Gamma^eff (dg/dx)

   isolates the Lagrangian/Z-vector algebra completely: a failure there is a
   defect in the response and in nothing else.  This runs before any molecule.
2. **Identities with no finite differences at all.**  ``dL/dx = sum_I w_I
   dE_I/dx`` is exact, and the effective generalized Fock must be symmetric --
   which in the active-active block holds only because the CI multipliers were
   solved.  Both are free of truncation error, so they fail sharply.
3. **Negative controls.**  The same gradient with the response removed must
   disagree with the reference by orders of magnitude.  Without this, a test
   suite cannot distinguish "the response term is applied" from "the response
   term is computed and discarded".
4. **Molecules**, against a five-point O(h^4) finite difference, plus
   translational invariance -- the sharpest cheap check on the Pulay and
   two-particle terms, which the Hellmann-Feynman part alone would not satisfy.
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False))


def _ao_gradient_available() -> bool:
    if not _backend_available():
        return False
    try:
        from oqp.library.casscf_sa_gradient import _ao_gradient_backend
    except Exception:
        return False
    return _ao_gradient_backend() is not None


# =========================================================== the abstract model
class _Model:
    """A smooth one-parameter family of orthonormal-basis integrals.

    ``g`` is ``sum_a L^a_pq L^a_rs`` from symmetric factors, which gives it the
    full eight-fold permutational symmetry of a real ERI tensor and makes it
    positive semidefinite, so the active spectrum behaves like a real one.  Both
    ``h`` and ``L`` are affine in ``x``, so the derivative integrals are exact
    rather than differenced -- the finite difference is spent on the CONVERGED
    energy, which is the only place it is needed.
    """

    def __init__(self, nbf, nfac=4, seed=20260816, coupling=0.35):
        rng = np.random.default_rng(seed)

        def sym(scale):
            m = rng.normal(scale=scale, size=(nbf, nbf))
            return 0.5 * (m + m.T)

        self.nbf = nbf
        self.h0 = sym(1.0) - np.diag(np.arange(nbf, dtype=float) * 1.5 + 3.0)
        self.h1 = sym(coupling)
        self.l0 = np.array([sym(0.45) for _ in range(nfac)])
        self.l1 = np.array([sym(coupling * 0.30) for _ in range(nfac)])

    def integrals(self, x):
        lmat = self.l0 + x * self.l1
        return self.h0 + x * self.h1, np.einsum("apq,ars->pqrs", lmat, lmat)

    def derivatives(self, x):
        lmat = self.l0 + x * self.l1
        return self.h1, (np.einsum("apq,ars->pqrs", self.l1, lmat)
                         + np.einsum("apq,ars->pqrs", lmat, self.l1))


def _expm_antisym(kmat):
    n = max(int(np.ceil(np.log2(max(np.max(np.abs(kmat)), 1e-300)) + 6)), 0)
    a = kmat / (2.0 ** n)
    out = np.eye(kmat.shape[0])
    term = np.eye(kmat.shape[0])
    for k in range(1, 18):
        term = term @ a / k
        out = out + term
    for _ in range(n):
        out = out @ out
    return out


def _solve_ci(h, g, ncore, nact, nelec, nroot):
    from oqp.library.casscf_hessian import (_active_hamiltonian,
                                            _excitation_matrices, _fold_active)
    from oqp.library.fci import _symmetric_eigh

    f, g_act = _fold_active(h, g, ncore, nact)
    dets, stack = _excitation_matrices(nact, *nelec)
    eps, vecs = _symmetric_eigh(_active_hamiltonian(f, g_act, stack, dets, nact))
    ecore = 0.0
    if ncore:
        c = slice(0, ncore)
        ecore = (2.0 * float(np.trace(h[c, c]))
                 + 2.0 * float(np.einsum("iijj->", g[c, c, c, c]))
                 - float(np.einsum("ijji->", g[c, c, c, c])))
    return eps[:nroot] + ecore, vecs[:, :nroot], dets


def _optimize(model, x, ncore, nact, nelec, weights, roots, coeff=None,
              tol=1.0e-12, maxit=300):
    """Drive the SA orbital gradient to zero with an |eigenvalue| Newton step.

    Taking the absolute value of the Hessian eigenvalues targets a STATIONARY
    point of any index, which is what a state average needs: the converged
    solution is routinely a saddle of the weighted objective in the orbital
    space and a descent method would walk away from it.
    """
    from oqp.library.casscf import _kappa_matrix, _nonredundant_pairs
    from oqp.library.casscf_hessian import analytic_orbital_hessian
    from oqp.library.casscf_sa_gradient import _gfock_and_grad
    from oqp.library.fci import _symmetric_eigh, _transform_integrals
    from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

    nbf = model.nbf
    h_ref, g_ref = model.integrals(x)
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    coeff = np.eye(nbf) if coeff is None else np.array(coeff, copy=True)
    for _ in range(maxit):
        h1e, eri = _transform_integrals(h_ref, g_ref, coeff)
        energies, ci, dets = _solve_ci(h1e, eri, ncore, nact, nelec,
                                       max(roots) + 1)
        gammas = np.array([make_rdm1_spatial(ci[:, r], dets, nact) for r in roots])
        gammas2 = np.array([make_rdm2_spatial(ci[:, r], dets, nact) for r in roots])
        _f, grad = _gfock_and_grad(nbf, ncore, nact, weights, gammas, gammas2,
                                   h1e, eri, pairs)
        if np.linalg.norm(grad) < tol:
            return dict(coeff=coeff, h1e=h1e, eri=eri, energies=energies,
                        ci=ci, dets=dets, gammas=gammas, gammas2=gammas2,
                        pairs=pairs, gnorm=float(np.linalg.norm(grad)))
        hess = analytic_orbital_hessian(h1e, eri, ncore, nact, nelec, pairs,
                                        weights, roots, ci)
        ev, evec = _symmetric_eigh(hess)
        step = -evec @ ((evec.T @ grad) / np.maximum(np.abs(ev), 1.0e-8))
        norm = float(np.linalg.norm(step))
        if norm > 0.25:
            step *= 0.25 / norm
        coeff = coeff @ _expm_antisym(_kappa_matrix(step, pairs, nbf))
    raise AssertionError(f"abstract SA-CASSCF did not converge, |g| = "
                         f"{np.linalg.norm(grad):.3e}")


def _abstract_analytic(sol, model, x, ncore, nact, nelec, weights, roots,
                       target, response="full"):
    from oqp.library.casscf_sa_gradient import (_gfock_and_grad,
                                                averaged_objective_densities,
                                                check_densities,
                                                relaxed_state_densities)
    from oqp.library.fci import _transform_integrals

    dh, dg = model.derivatives(x)
    dh_c, dg_c = _transform_integrals(dh, dg, sol["coeff"])
    if target is None:
        built = averaged_objective_densities(ncore, nact, model.nbf, weights,
                                             sol["gammas"], sol["gammas2"])
        built["fock"] = _gfock_and_grad(model.nbf, ncore, nact, weights,
                                        sol["gammas"], sol["gammas2"],
                                        sol["h1e"], sol["eri"], sol["pairs"])[0]
    else:
        built = relaxed_state_densities(
            sol["h1e"], sol["eri"], ncore, nact, nelec, sol["pairs"], weights,
            roots, target, sol["ci"], sol["dets"], sol["gammas"],
            sol["gammas2"], response=response)
    value = (float(np.einsum("pq,pq->", built["d_eff"], dh_c))
             + 0.5 * float(np.einsum("pqrs,pqrs->", built["g_eff"], dg_c)))
    return value, built, check_densities


def _abstract_numeric(model, x, ncore, nact, nelec, weights, roots, target,
                      coeff, step=2.0e-3):
    """Five-point O(h^4) difference of the RE-CONVERGED energy.

    Each displaced point restarts from the reference orbitals, which is the
    cheapest way to stay on one solution branch: a state average has many
    stationary points and a cold start is free to find a different one, with
    nothing in the number saying that it did.
    """
    total = 0.0
    for off, c in zip((-2, -1, 1, 2),
                      (1 / 12.0, -2 / 3.0, 2 / 3.0, -1 / 12.0)):
        sol = _optimize(model, x + off * step, ncore, nact, nelec, weights,
                        roots, coeff=coeff)
        e = (float(np.dot(weights, sol["energies"][roots])) if target is None
             else float(sol["energies"][target]))
        total += c * e
    return total / step


_ABSTRACT_CASES = [
    pytest.param(6, 1, 3, (2, 2), [0.5, 0.5], [0, 1], 20260816,
                 id="equal-weights"),
    pytest.param(6, 1, 3, (2, 2), [0.7, 0.3], [0, 1], 20260816,
                 id="unequal-weights"),
    pytest.param(6, 1, 3, (2, 2), [0.5, 0.3, 0.2], [0, 1, 2], 90210,
                 id="three-states"),
    pytest.param(5, 0, 3, (2, 2), [0.6, 0.4], [0, 1], 4242,
                 id="no-inactive"),
    pytest.param(6, 1, 3, (2, 2), [0.5, 0.5], [0, 2], 20260816,
                 id="non-contiguous-roots"),
]


@pytest.mark.parametrize("nbf,ncore,nact,nelec,wts,roots,seed", _ABSTRACT_CASES)
def test_abstract_sa_gradients_match_finite_difference(
        nbf, ncore, nact, nelec, wts, roots, seed):
    """Both paths, on an abstract model with no basis set and no Pulay term."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    model = _Model(nbf, seed=seed)
    weights = np.asarray(wts, dtype=float)
    weights /= weights.sum()
    sol = _optimize(model, 0.0, ncore, nact, nelec, weights, roots)
    assert sol["gnorm"] < 1.0e-11

    for target in [None] + list(roots):
        ana, built, check = _abstract_analytic(
            sol, model, 0.0, ncore, nact, nelec, weights, roots, target)
        num = _abstract_numeric(model, 0.0, ncore, nact, nelec, weights, roots,
                                target, sol["coeff"])
        assert ana == pytest.approx(num, abs=5.0e-8), (
            f"target={target}: analytic {ana:.12f} vs numeric {num:.12f}")
        check(built["d_eff"], built["g_eff"], built["fock"],
              built["sep_terms"], built["low_terms"], sol["eri"])


@pytest.mark.parametrize("nbf,ncore,nact,nelec,wts,roots,seed", _ABSTRACT_CASES)
def test_abstract_weighted_sum_rule(nbf, ncore, nact, nelec, wts, roots, seed):
    """``dL/dx = sum_I w_I dE_I/dx`` -- exact, with no finite differences.

    The weights are constants, so this identity has no truncation error at all
    and holds to machine precision.  It is the sharpest available statement
    that the two paths are derivatives of the same object, and it would break
    immediately if either the averaged densities or the Z-vector were wrong."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    model = _Model(nbf, seed=seed)
    weights = np.asarray(wts, dtype=float)
    weights /= weights.sum()
    sol = _optimize(model, 0.0, ncore, nact, nelec, weights, roots)
    args = (sol, model, 0.0, ncore, nact, nelec, weights, roots)
    avg = _abstract_analytic(*args, None)[0]
    parts = [_abstract_analytic(*args, r)[0] for r in roots]
    assert avg == pytest.approx(float(np.dot(weights, parts)), abs=1.0e-11)


def test_abstract_response_term_is_not_negligible():
    """Negative control: remove the response and the answer must break.

    Both removals are tested because they fail in different ways.  Dropping the
    whole Z-vector leaves the state-specific formula, which is wrong by a first-
    order term.  Dropping only the CI multipliers -- the plausible-looking
    mistake -- leaves a gradient that is still wrong AND breaks the symmetry of
    the effective generalized Fock, because the active-active block of that
    matrix is annihilated by CI stationarity and by nothing else.  The guard
    must catch the second case on its own, without a finite difference."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    from oqp.library.casscf_sa_gradient import (check_densities,
                                                relaxed_state_densities)

    nbf, ncore, nact, nelec = 6, 1, 3, (2, 2)
    roots = [0, 1]
    weights = np.array([0.7, 0.3])
    model = _Model(nbf, seed=20260816)
    sol = _optimize(model, 0.0, ncore, nact, nelec, weights, roots)

    for target in roots:
        num = _abstract_numeric(model, 0.0, ncore, nact, nelec, weights, roots,
                                target, sol["coeff"])
        full = _abstract_analytic(sol, model, 0.0, ncore, nact, nelec, weights,
                                  roots, target)[0]
        bare = _abstract_analytic(sol, model, 0.0, ncore, nact, nelec, weights,
                                  roots, target, response="none")[0]
        assert abs(full - num) < 5.0e-8
        assert abs(bare - num) > 1.0e-3, (
            "dropping the Z-vector changed nothing: the response term is being "
            "computed but not applied")

        built = relaxed_state_densities(
            sol["h1e"], sol["eri"], ncore, nact, nelec, sol["pairs"], weights,
            roots, target, sol["ci"], sol["dets"], sol["gammas"],
            sol["gammas2"], response="orbital")
        with pytest.raises(ValueError, match="not symmetric"):
            check_densities(built["d_eff"], built["g_eff"], built["fock"],
                            built["sep_terms"], built["low_terms"], sol["eri"])


def test_abstract_near_degenerate_roots():
    """Unequal weights with the two averaged roots close together.

    This is the regime state averaging exists for and the one the formulation
    is most fragile in: with unequal weights the CI response denominators ARE
    the root gap.  The identities that carry no truncation error -- the sum
    rule and the Fock symmetry -- must survive it unchanged, and they are what
    is asserted here; the finite difference itself becomes the limiting factor
    long before the analytic gradient does."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    nbf, ncore, nact, nelec = 6, 1, 3, (2, 2)
    roots = [0, 1]
    weights = np.array([0.65, 0.35])

    # An avoided crossing in a model with no symmetry has a nonzero minimum
    # gap, so a small gap is SEARCHED for rather than assumed.
    best = None
    for seed in (777, 4711):
        model = _Model(nbf, seed=seed, coupling=0.9)
        for x in np.linspace(-1.5, 1.5, 25):
            try:
                sol = _optimize(model, float(x), ncore, nact, nelec, weights,
                                roots)
            except (AssertionError, ValueError):
                continue
            gap = float(sol["energies"][1] - sol["energies"][0])
            if best is None or gap < best[0]:
                best = (gap, float(x), sol, model)
    assert best is not None
    gap, x_star, sol, model = best
    assert gap < 0.1, (
        f"the search did not find a near-degeneracy (best gap {gap:.3e}); the "
        "test is not exercising the regime it claims to")

    args = (sol, model, x_star, ncore, nact, nelec, weights, roots)
    avg, built, check = _abstract_analytic(*args, None)
    parts = []
    for r in roots:
        val, b, _c = _abstract_analytic(*args, r)
        parts.append(val)
        check(b["d_eff"], b["g_eff"], b["fock"], b["sep_terms"],
              b["low_terms"], sol["eri"])
    assert avg == pytest.approx(float(np.dot(weights, parts)), abs=1.0e-10)

    num = _abstract_numeric(model, x_star, ncore, nact, nelec, weights, roots,
                            1, sol["coeff"], step=min(2.0e-3, 0.05 * gap))
    assert parts[1] == pytest.approx(num, abs=max(1.0e-7, 1.0e-5 * abs(num)))

    # With EQUAL weights the averaged-pair coupling cancels, so the SA orbital
    # Hessian is perfectly well defined at the same near-crossing and says
    # nothing about it.  The individual state is still not differentiable at a
    # crossing, and only the gradient's own root-gap guard catches that -- which
    # is why the guard exists as a separate layer.  Here the nearest CI root to
    # the target IS the other averaged root, so a tolerance just above their gap
    # reaches that branch without disturbing the external couplings.
    from oqp.library.casscf_sa_gradient import relaxed_state_densities
    equal = np.array([0.5, 0.5])
    sol_eq = _optimize(model, x_star, ncore, nact, nelec, equal, roots,
                       coeff=sol["coeff"])
    gap_eq = float(sol_eq["energies"][1] - sol_eq["energies"][0])
    common = (sol_eq["h1e"], sol_eq["eri"], ncore, nact, nelec, sol_eq["pairs"],
              equal, roots, 1, sol_eq["ci"], sol_eq["dets"], sol_eq["gammas"],
              sol_eq["gammas2"])
    relaxed_state_densities(*common, degen_tol=0.1 * gap_eq)    # admitted
    with pytest.raises(ValueError, match="degenerate with another CI root"):
        relaxed_state_densities(*common, degen_tol=2.0 * gap_eq)


def test_abstract_degenerate_roots_are_refused():
    """A crossing is refused, not regularized.

    The adiabatic energy is not differentiable at an exact crossing whatever
    the weights, and the CI response denominator is the gap itself.  Forcing
    the tolerance above the real gap is the portable way to reach that branch."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    from oqp.library.casscf_sa_gradient import relaxed_state_densities

    nbf, ncore, nact, nelec = 6, 1, 3, (2, 2)
    roots = [0, 1]
    model = _Model(nbf, seed=20260816)

    # Unequal weights: the CI response denominator between the two averaged
    # roots is live, so the SA orbital Hessian itself is undefined and refuses
    # first -- the coupled system has no solution to hand back.
    weights = np.array([0.7, 0.3])
    sol = _optimize(model, 0.0, ncore, nact, nelec, weights, roots)
    with pytest.raises(ValueError, match="degener"):
        relaxed_state_densities(
            sol["h1e"], sol["eri"], ncore, nact, nelec, sol["pairs"], weights,
            roots, 1, sol["ci"], sol["dets"], sol["gammas"], sol["gammas2"],
            degen_tol=10.0)


# ================================================================== the molecule
_LIH = "\nLi 0.0 0.0 0.0\nH 0.0 0.0 1.6"
_H4 = ("\nH 0.0  0.041 0.0"
       "\nH 0.017 0.0   0.795"
       "\nH -0.023 0.062 1.523"
       "\nH 0.009 -0.031 2.310")


def _sa_config(system, cas, nroot, roots, weights=None, gradient_state="averaged",
               runtype="grad", basis="sto-3g"):
    sa = {"enabled": "true", "nstate": str(len(roots)),
          "target_roots": ",".join(str(r) for r in roots)}
    if weights is None:
        sa["equal_weights"] = "true"
    else:
        sa["equal_weights"] = "false"
        sa["weights"] = ",".join(f"{w}" for w in weights)
    return {
        "input": {"system": system, "charge": "0", "basis": basis,
                  "method": "sa-casscf", "runtype": runtype},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "conv": "1.0e-10",
                "maxit": "80", "forced_attempt": "3", "save_molden": "False"},
        "properties": {"grad": "0"},
        "cas": cas,
        "ci": {"nroot": str(nroot), "solver": "dense", "eig_tol": "1.0e-10",
               "integral_backend": "native", "target_spin": "singlet",
               "root_tracking": "energy"},
        "casscf": {"max_macro_iterations": "150", "optimizer": "newton",
                   "gradient_norm_tol": "1.0e-9",
                   "energy_decrease_tol": "1.0e-12",
                   "step_norm_tol": "1.0e-10", "max_rotation_norm": "2.0e-1",
                   "level_shift": "1.0e-3", "canonicalize": "true",
                   "gradient_state": str(gradient_state)},
        "state_average": sa,
        "tests": {"exception": "True"},
    }


_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
          "max_det": "5000", "orbital_source": "rhf", "sort_orbitals": "energy"}


def _runner(tmp_path, project, config):
    from oqp.pyoqp import Runner
    return Runner(project=project, input_file=None,
                  log=str(tmp_path / f"{project}.log"), input_dict=config,
                  silent=1, usempi=False)


_RESTART_TAGS = ('OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::DM_A', 'OQP::DM_B',
                 'OQP::E_MO_A', 'OQP::E_MO_B')


def _sa_energy(mol, target, roots):
    """The differentiated scalar: the objective, or one averaged root."""
    if target is None:
        return float(mol.mol_energy.energy)
    return float(mol.energies[[int(r) for r in roots].index(int(target))])


def _fd_along(mol, sp, target, roots, k, x0, e0, step):
    """O(h^4) derivative along coordinate ``k``, and the energy smoothness.

    A state average has many stationary points and a displaced run has an
    independent chance of landing on a different one, so the displaced solves
    restart from the reference orbitals and the sampled energies are checked
    for smoothness.  Without that, a branch change reads as a gradient defect.
    """
    snap = {}
    for tag in _RESTART_TAGS:
        try:
            snap[tag] = np.array(mol.data[tag], dtype=float, copy=True)
        except Exception:
            pass

    energies = {}
    for m in (-2, -1, 1, 2):
        x = x0.copy()
        x[k] += m * step
        for tag, val in snap.items():
            try:
                mol.data[tag][...] = val
            except Exception:
                pass
        mol.update_system(x)
        sp.energy()
        energies[m] = _sa_energy(mol, target, roots)

    for tag, val in snap.items():
        try:
            mol.data[tag][...] = val
        except Exception:
            pass
    mol.update_system(x0)
    sp.energy()

    deriv = (energies[-2] - 8 * energies[-1]
             + 8 * energies[1] - energies[2]) / (12 * step)
    xs = np.array([-2, -1, 0, 1, 2], dtype=float) * step
    ys = np.array([energies[-2], energies[-1], e0, energies[1], energies[2]])
    smooth = float(np.max(np.abs(ys - np.polyval(np.polyfit(xs, ys, 3), xs))))
    return deriv, smooth


def _sa_gradient(mol, target):
    """Gradient for one selector, reusing the converged wavefunction."""
    from oqp.library.casscf_sa_gradient import sa_casscf_gradient
    mol.config["casscf"]["gradient_state"] = (
        "averaged" if target is None else str(target))
    return sa_casscf_gradient(mol)


@pytest.mark.parametrize("weights,ids", [(None, "equal"), ([0.7, 0.3], "uneq")])
def test_sa_casscf_gradients_match_finite_difference(tmp_path, weights, ids):
    """Both paths on LiH/STO-3G CAS(2,2), against an O(h^4) difference."""
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    from oqp.library.single_point import SinglePoint

    config = _sa_config(_LIH, _CAS22, 2, [0, 1], weights=weights)
    runner = _runner(tmp_path, f"sagrad_{ids}", config)
    mol = runner.mol
    sp = SinglePoint(mol)
    sp.energy()
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()

    for target in (None, 0, 1):
        grad, diag = _sa_gradient(mol, target)
        assert diag["gnorm_sa"] < 1.0e-7
        assert diag["fock_asymmetry"] < 1.0e-7
        assert diag["factorization_mismatch"] < 1.0e-8
        e0 = _sa_energy(mol, target, [0, 1])
        flat = np.asarray(grad, dtype=float).reshape(-1)
        k = int(np.argmax(np.abs(flat)))
        deriv, smooth = _fd_along(mol, sp, target, [0, 1], k, x0, e0, 1.0e-3)
        assert smooth < 1.0e-8, (
            f"target={target}: the displaced SA-CASSCF energies are not smooth "
            f"(cubic-fit residual {smooth:.2e}), so the finite difference is "
            "not a derivative -- a displaced point converged to a different "
            "state-averaged solution")
        assert flat[k] == pytest.approx(deriv, abs=5.0e-7), (
            f"target={target}, coordinate {k}")


def test_sa_casscf_weighted_sum_rule_molecular(tmp_path):
    """``dL/dx = sum_I w_I dE_I/dx`` on a real molecule, exactly.

    No finite differences: the weights are constants.  This ties the two
    production paths -- including the AO transform, the factorized two-particle
    density and the derivative-integral contraction, which the abstract
    protocol cannot reach -- to one another."""
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    from oqp.library.single_point import SinglePoint

    config = _sa_config(_LIH, _CAS22, 2, [0, 1], weights=[0.7, 0.3])
    runner = _runner(tmp_path, "sagrad_sumrule", config)
    mol = runner.mol
    SinglePoint(mol).energy()

    avg = _sa_gradient(mol, None)[0]
    parts = [_sa_gradient(mol, r)[0] for r in (0, 1)]
    combined = 0.7 * parts[0] + 0.3 * parts[1]
    assert np.max(np.abs(avg - combined)) < 1.0e-9


def test_sa_casscf_gradients_are_translationally_invariant(tmp_path):
    """A rigid shift cannot change any energy, averaged or individual.

    The sharpest cheap check on the Pulay and two-particle terms: a missing or
    mis-signed derivative-integral contribution breaks it while the
    Hellmann-Feynman part alone would not.  It also covers the relaxed
    densities, which enter the same contraction."""
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    from oqp.library.single_point import SinglePoint

    config = _sa_config(_H4, _CAS22, 2, [0, 1])
    runner = _runner(tmp_path, "sagrad_trans", config)
    mol = runner.mol
    SinglePoint(mol).energy()

    for target in (None, 0, 1):
        grad = np.asarray(_sa_gradient(mol, target)[0], dtype=float)
        assert np.max(np.abs(grad.sum(axis=0))) < 1.0e-9, f"target={target}"


def test_sa_casscf_individual_roots_differ(tmp_path):
    """Root selection is real: the two averaged roots have different gradients.

    Pins the part that is easiest to get silently wrong -- returning the same
    quantity whatever ``gradient_state`` says -- and shows the averaged
    objective is neither of them."""
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    from oqp.library.single_point import SinglePoint

    config = _sa_config(_LIH, _CAS22, 2, [0, 1])
    runner = _runner(tmp_path, "sagrad_roots", config)
    mol = runner.mol
    SinglePoint(mol).energy()

    g_avg = _sa_gradient(mol, None)[0]
    g0 = _sa_gradient(mol, 0)[0]
    g1 = _sa_gradient(mol, 1)[0]
    assert np.max(np.abs(g0 - g1)) > 1.0e-3
    assert np.max(np.abs(g_avg - g0)) > 1.0e-4
    assert np.max(np.abs(g_avg - g1)) > 1.0e-4


def test_one_state_average_reproduces_the_state_specific_gradient(tmp_path):
    """A one-state "average" must give the validated state-specific gradient.

    This is the sharpest available check on everything the abstract protocol
    cannot reach -- the AO transform, the factorized two-particle density, the
    Pulay term and the new derivative-integral entry point -- because the
    reference is an INDEPENDENT implementation of the same number
    (``casscf_gradient.F90``, which builds its own AO density and runs its own
    contraction), not a rearrangement of this one.

    Both paths are checked.  Path A is the algebraic identity: with one weight
    of 1.0 the averaged objective IS that state.  Path B is stronger: it runs
    the whole Z-vector machinery, which must then find ``g^J = g^SA = 0``,
    solve for ``z = 0``, and reduce to the same answer.  A response term with a
    sign or scale error would generally survive Path A and fail here.
    """
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    import oqp
    if not hasattr(oqp, "casscf_gradient"):
        pytest.skip("liboqp has no casscf_gradient reference entry point")

    from oqp.library.casscf_gradient import casscf_analytic_gradient
    from oqp.library.single_point import SinglePoint

    ss = {
        "input": {"system": _LIH, "charge": "0", "basis": "sto-3g",
                  "method": "casscf", "runtype": "grad"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "conv": "1.0e-10",
                "maxit": "80", "forced_attempt": "3", "save_molden": "False"},
        "properties": {"grad": "0"},
        "cas": _CAS22,
        "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-10",
               "integral_backend": "native", "target_spin": "singlet"},
        "casscf": {"max_macro_iterations": "150", "optimizer": "newton",
                   "root": "0", "gradient_norm_tol": "1.0e-9",
                   "energy_decrease_tol": "1.0e-12", "step_norm_tol": "1.0e-10",
                   "max_rotation_norm": "2.0e-1", "level_shift": "1.0e-3",
                   "canonicalize": "true"},
        "tests": {"exception": "True"},
    }
    runner = _runner(tmp_path, "sagrad_ss", ss)
    SinglePoint(runner.mol).energy()
    reference = np.asarray(casscf_analytic_gradient(runner.mol),
                           dtype=float).reshape(-1, 3)

    sa = _sa_config(_LIH, _CAS22, 1, [0])
    runner = _runner(tmp_path, "sagrad_one", sa)
    mol = runner.mol
    SinglePoint(mol).energy()

    g_avg, diag_avg = _sa_gradient(mol, None)
    assert np.max(np.abs(np.asarray(g_avg) - reference)) < 1.0e-9

    g_root, diag_root = _sa_gradient(mol, 0)
    assert diag_root["zvector"]["z_norm"] < 1.0e-7, (
        "with a single averaged state the individual gradient IS the objective "
        "gradient, so the Z-vector must come out zero")
    assert np.max(np.abs(np.asarray(g_root) - reference)) < 1.0e-9
    assert diag_root["fock_asymmetry"] < 1.0e-8
    assert diag_avg["factorization_mismatch"] < 1.0e-8


def test_sa_response_term_is_not_negligible_molecular(tmp_path):
    """Negative control on a molecule: remove the response, break the answer.

    The comparison is against the finite difference of the SA-CASSCF energy
    itself, so this states that the Z-vector is not a small correction to the
    published number but the difference between right and wrong."""
    if not _ao_gradient_available():
        pytest.skip("liboqp has no casscf_ao_gradient; rebuild to run this test")

    from oqp.library.casscf_sa_gradient import sa_casscf_gradient
    from oqp.library.single_point import SinglePoint

    config = _sa_config(_LIH, _CAS22, 2, [0, 1], weights=[0.7, 0.3])
    runner = _runner(tmp_path, "sagrad_neg", config)
    mol = runner.mol
    sp = SinglePoint(mol)
    sp.energy()
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()

    mol.config["casscf"]["gradient_state"] = "1"
    full = np.asarray(sa_casscf_gradient(mol)[0], dtype=float).reshape(-1)
    bare = np.asarray(sa_casscf_gradient(mol, response="none")[0],
                      dtype=float).reshape(-1)

    e0 = _sa_energy(mol, 1, [0, 1])
    k = int(np.argmax(np.abs(full)))
    deriv, smooth = _fd_along(mol, sp, 1, [0, 1], k, x0, e0, 1.0e-3)
    assert smooth < 1.0e-8
    assert full[k] == pytest.approx(deriv, abs=5.0e-7)
    assert abs(bare[k] - deriv) > 1.0e-4, (
        "dropping the Z-vector changed nothing: the response term is being "
        "computed but not applied")


# ==================================================================== selection
def test_gradient_state_selector_rejects_unaveraged_root():
    """A root the objective never averaged over has no Z-vector to build."""
    from oqp.library.casscf_sa_gradient import resolve_gradient_state

    assert resolve_gradient_state({"casscf": {}}, [0, 1]) is None
    assert resolve_gradient_state({"casscf": {"gradient_state": "averaged"}},
                                  [0, 1]) is None
    assert resolve_gradient_state({"casscf": {"gradient_state": "1"}},
                                  [0, 1]) == 1
    with pytest.raises(ValueError, match="not one of the"):
        resolve_gradient_state({"casscf": {"gradient_state": "3"}}, [0, 1])
    with pytest.raises(ValueError, match="must be"):
        resolve_gradient_state({"casscf": {"gradient_state": "lowest"}}, [0, 1])


def test_state_specific_entry_point_refuses_a_state_average():
    """The state-specific formula must not be reused on an averaged run."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    from oqp.library.casscf_gradient import _state_average_requested

    class _Mol:
        config = {"input": {"method": "sa-casscf"}}

    class _Settings:
        state_average_enabled = False

    assert _state_average_requested(_Mol(), _Settings())

# ===================================================================== preflight
def _check(config):
    from oqp.utils.input_checker import check_input_values
    return check_input_values(config, raise_error=False, emit=False)


def _preflight_config(runtype="grad", gradient_state="averaged", grad="0",
                      istate=None, roots="0,1"):
    cfg = {
        "input": {"method": "sa-casscf", "runtype": runtype, "system": _LIH,
                  "basis": "sto-3g", "charge": 0},
        "scf": {"type": "rhf", "multiplicity": 1},
        "properties": {"grad": grad},
        "cas": {"active_electrons": 2, "active_orbitals": 2, "frozen_core": 1},
        "ci": {"nroot": 2, "solver": "dense", "target_spin": "singlet"},
        "casscf": {"gradient_state": gradient_state},
        "state_average": {"enabled": "true", "nstate": 2,
                          "target_roots": roots, "equal_weights": "true"},
    }
    if istate is not None:
        cfg["optimize"] = {"istate": istate}
    return cfg


def _errors(report, field):
    return [e for e in report.to_text().splitlines() if field in e]


def test_preflight_admits_both_sa_gradient_paths():
    """A state-averaged gradient run is no longer rejected outright."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")
    for state in ("averaged", "0", "1"):
        report = _check(_preflight_config(gradient_state=state))
        assert report.ok, report.to_text()


def test_preflight_rejects_a_gradient_state_outside_the_average():
    """A root the objective never averaged over has no response to build."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")
    report = _check(_preflight_config(gradient_state="2"))
    assert not report.ok
    assert _errors(report, "casscf.gradient_state")


def test_preflight_rejects_a_properties_grad_that_contradicts_the_selector():
    """The two selectors must not disagree; the run must not pick one silently."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")
    # the objective is not a state, so no state index is consistent with it
    assert not _check(_preflight_config(gradient_state="averaged", grad="1")).ok
    # naming a different root than the one being differentiated
    assert not _check(_preflight_config(gradient_state="1", grad="0,2")).ok
    # the conventional 0 and the differentiated root itself are both accepted
    assert _check(_preflight_config(gradient_state="1", grad="0,1")).ok
    assert _check(_preflight_config(gradient_state="1", grad="1")).ok


def test_preflight_refuses_the_objective_for_optimizer_runtypes():
    """An optimizer pairs energies[istate] with grads[istate].

    The weighted objective is not in that energy list, so a geometry
    optimization asked to follow it would minimize a state energy while
    stepping along the objective's gradient -- and nothing in the output would
    say so.  For an individual root the pairing is legal, but only when istate
    is the root's POSITION in target_roots."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    report = _check(_preflight_config(runtype="optimize",
                                      gradient_state="averaged", istate=0))
    assert not report.ok
    assert _errors(report, "casscf.gradient_state")

    assert _check(_preflight_config(runtype="optimize", gradient_state="1",
                                    istate=1)).ok
    bad = _check(_preflight_config(runtype="optimize", gradient_state="1",
                                   istate=0))
    assert not bad.ok
    assert _errors(bad, "optimize.istate")
