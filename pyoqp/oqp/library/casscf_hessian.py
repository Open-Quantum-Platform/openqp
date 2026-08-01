"""Analytic MCSCF orbital-rotation Hessian for the native CASSCF.

This module removes the finite-difference orbital Hessian bottleneck of
:mod:`oqp.library.casscf` (``_fd_orbital_hessian``: ``2*n_par`` gradient
evaluations = ``2*n_par`` active-space CI solves per build) by assembling the
exact second derivative of the (state-averaged) CASSCF energy with respect to
the non-redundant orbital rotations, using only the current MO integrals, the
CI vectors of the averaged roots, and one dense diagonalization of the active
Hamiltonian.  No calls to the macroiteration ``evaluate`` closure (i.e. no
extra CI solves) are made.

Conventions (identical to ``casscf.py``)
----------------------------------------
Orbitals are rotated ``C(x) = C exp(K(x))`` with

    K(x) = sum_k x_k K^(k),      K^(k) = e_{p_k q_k} - e_{q_k p_k},

over the non-redundant pair list produced by ``casscf._nonredundant_pairs``
(the exact ``pairs`` sequence is passed in, so the Hessian rows/columns use
the same ordering as ``_fd_orbital_hessian``): (active, inactive) then
(virtual, inactive) then (virtual, active), with the first pair index the
higher space.  MO integrals transform as ``h(x) = e^{-K} h e^{K}`` and the
chemist ERIs ``(pq|rs)(x)`` by the same one-body rotation on all four indices.

What the finite-difference Hessian actually measures
----------------------------------------------------
``_fd_orbital_hessian`` differentiates the *variational* orbital gradient: at
every displaced point the active-space CI is re-solved.  Its symmetrized limit
is therefore the second derivative of

    E(x) = sum_I w_I E_I(x),

where ``E_I(x)`` are eigenvalues of the core-folded active Hamiltonian

    Heff(x) = Ecore(x) + sum_tu f_tu(x) E_tu
              + 1/2 sum_tuvw (tu|vw)(x) (E_tu E_vw - delta_uv E_tw),

with ``f`` and ``Ecore`` the standard inactive-Fock folding of
``fci._active_space``.  By eigenvalue perturbation theory this second
derivative splits exactly into an orbital ("fixed-CI") part and a
CI-relaxation part:

    H_kl = sum_I w_I <I| d2 Heff / dx_k dx_l |I>
         + 2 sum_I w_I sum_{J != I} <I| dHeff/dx_k |J> <J| dHeff/dx_l |I>
                                     / (E_I - E_J),

where ``J`` runs over the COMPLETE eigenspectrum of ``Heff`` in the same
determinant basis the CI solver uses.  Both parts are assembled here; omitting
the second (as a naive fixed-CI orbital Hessian would) does NOT reproduce the
finite-difference Hessian away from trivial cases.

Part 1: fixed-CI orbital Hessian from the SA RDMs
-------------------------------------------------
With the full-space spin-summed state-averaged RDMs ``D`` (1-RDM) and ``G``
(chemist 2-RDM, ``E2 = 1/2 sum (pq|rs) G_pqrs``, built by
``casscf._full_rdms``), the fixed-CI energy is the linear functional

    E_fix(x) = sum_pq D_pq h_pq(x) + 1/2 sum_pqrs G_pqrs (pq|rs)(x),

and (because the CAS core-folding identity is linear in the integrals) its
second derivative equals ``sum_I w_I <I| d2 Heff |I>`` exactly, including the
``d2 Ecore`` scalar.  Directional derivatives of ``E_fix`` along an arbitrary
rotation generator ``K`` are ``sum_mn K_mn Z_mn`` with the intermediate

    Z(t, T)_mn = (D t)_nm - (t D)_nm
               + 1/2 sum_qrs [ G_nqrs T_mqrs + G_qnrs T_qmrs
                             + G_qrns T_qrms + G_qrsn T_qrsm ],

evaluated on integral-like tensors ``(t, T)``.  For ``(t, T) = (h, g)`` the
pair contraction ``Z_{p q} - Z_{q p}`` reproduces the production gradient
``2 (F_qp - F_pq)`` of ``casscf.py`` (asserted in the tests).  The Hessian is
the symmetrized second directional derivative

    H_fix = 1/2 (B + B^T),
    B_kl  = Z(h^(l), g^(l))_{p_k q_k} - Z(h^(l), g^(l))_{q_k p_k},

built from the one-index derivative integrals

    h^(l) = [h, K^(l)],
    g^(l)_pqrs = sum_m [ K^(l)_mp g_mqrs + K^(l)_mq g_pmrs
                       + K^(l)_mr g_pqms + K^(l)_ms g_pqrm ].

The symmetrization matches the explicit ``0.5 (H + H^T)`` of
``_fd_orbital_hessian``; the antisymmetric parts of both are the same
gradient-weighted commutator term and cancel identically.  This is
algebraically the standard MCSCF orbital-orbital Hessian (Siegbahn, Almloef,
Heiberg, Roos, J. Chem. Phys. 74 (1981) 2384; Helgaker, Jorgensen, Olsen,
"Molecular Electronic-Structure Theory", Sec. 10.8), assembled through
derivative-integral / generalized-Fock ("Y intermediate") contractions rather
than per-block closed formulas so that one code path covers all rotation
classes and the SA case.

Part 2: CI relaxation from the active spectrum
----------------------------------------------
``dHeff/dx_k`` is the effective Hamiltonian of the core-folded derivative
integrals (folding is linear): ``f^(k), g^(k)_act`` from ``(h^(k), g^(k))``.
Its matrix elements are evaluated with precomputed spin-summed excitation
matrices ``E_tu`` over the same determinant list ``fci._determinants`` the CI
solver uses; the full spectrum ``(E_J, |J>)`` comes from one dense
``_symmetric_eigh`` of ``Heff`` (the honest cost of this module -- one dense
active-space diagonalization per Hessian, replacing ``2*n_par`` CI solves).
The scalar ``dEcore/dx_k`` drops from all off-diagonal elements.  Degeneracy
handling: couplings between two *averaged* states enter as
``0.5 (w_I - w_J)`` per ordered visit, so the equal-weight SA-CASSCF terms
cancel exactly (the standard invariance of equal-weight state averaging); a
near-degenerate ``E_J ~ E_I`` with a genuinely non-zero coupling and non-zero
weight factor makes the SA objective non-smooth, and this module raises
instead of returning a garbage Hessian (the FD Hessian is equally undefined
there).  Spin-forbidden couplings vanish identically because the ``E_tu`` are
spin-summed and are skipped by the noise threshold.

Scope / cost limits (honest)
----------------------------
Dense-spectrum CI relaxation costs ``O(ndet^3)`` time and
``O(nact^2 ndet^2)`` memory for the excitation-matrix stack; a guard raises
above ``_MAX_STACK_BYTES`` with the advice to fall back to
``[casscf] hessian = fd``.  This matches the validation-grade scope of the
native CASSCF (small/medium active spaces).  Large active spaces need
iterative linear-response solves for part 2 (and a density-direct part 1),
which is the documented upgrade hook, not what is implemented here.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb

import numpy as np

from oqp.library.fci import _determinants, _symmetric_eigh, _transform_integrals
from oqp.library.rdm import annihilate, create, make_rdm1_spatial, make_rdm2_spatial

__all__ = ["analytic_orbital_hessian", "make_hessian_provider"]

# |E_I - E_J| below which a response denominator is treated as degenerate.
_DEGENERACY_TOL = 1.0e-8
# Coupling amplitudes below this are numerical zeros (spin-forbidden partners
# come out at ~1e-14 * integral scale); genuine couplings are orders larger.
_COUPLING_NOISE = 1.0e-8
# Memory guard for the (nact, nact, ndet, ndet) excitation-matrix stack.
_MAX_STACK_BYTES = 2 * 1024**3


# ----------------------------------------------------------------- derivative integrals
def _one_index_derivative_h(h, p, q):
    """``[h, K]`` for the single-pair generator ``K = e_pq - e_qp``."""
    t = np.zeros_like(h)
    t[:, q] += h[:, p]
    t[:, p] -= h[:, q]
    t[p, :] -= h[q, :]
    t[q, :] += h[p, :]
    return t


def _one_index_derivative_g(g, p, q):
    """One-index derivative of the chemist ERI tensor along ``K = e_pq - e_qp``."""
    T = np.zeros_like(g)
    T[q, :, :, :] += g[p, :, :, :]
    T[p, :, :, :] -= g[q, :, :, :]
    T[:, q, :, :] += g[:, p, :, :]
    T[:, p, :, :] -= g[:, q, :, :]
    T[:, :, q, :] += g[:, :, p, :]
    T[:, :, p, :] -= g[:, :, q, :]
    T[:, :, :, q] += g[:, :, :, p]
    T[:, :, :, p] -= g[:, :, :, q]
    return T


def _z_matrix(D, G, t, T):
    """Directional-derivative intermediate: ``dE_fix along K = sum_mn K_mn Z_mn``.

    ``Z(h, g)`` contracted over a pair reproduces the production gradient
    ``2 (F_qp - F_pq)``; evaluated on derivative integrals it yields the
    fixed-CI Hessian columns.  No permutational symmetry of ``t``/``T`` is
    assumed (derivative integrals are not full-ERI-symmetric)."""
    z = (D @ t - t @ D).T
    z += 0.5 * (
        np.einsum("nqrs,mqrs->mn", G, T, optimize=True)
        + np.einsum("qnrs,qmrs->mn", G, T, optimize=True)
        + np.einsum("qrns,qrms->mn", G, T, optimize=True)
        + np.einsum("qrsn,qrsm->mn", G, T, optimize=True)
    )
    return z


def _pair_differences(Z, pairs):
    return np.array([Z[p, q] - Z[q, p] for (p, q) in pairs])


def _fold_active(t, T, ncore, nact):
    """Inactive-Fock core folding of integral-like tensors (sequential spaces).

    Mirrors ``fci._active_space`` exactly:
    ``f[t,u] = t[t,u] + sum_i [2 T[t,u,i,i] - T[t,i,i,u]]``; linear in the
    integrals, hence valid for derivative integrals as well."""
    act = slice(ncore, ncore + nact)
    f = np.array(t[act, act], dtype=float, copy=True)
    if ncore:
        f += 2.0 * np.einsum("pqii->pq", T[act, act, :ncore, :ncore])
        f -= np.einsum("piiq->pq", T[act, :ncore, :ncore, act])
    g = np.array(T[act, act, act, act], dtype=float, copy=True)
    return f, g


# ----------------------------------------------------------------- determinant machinery
@lru_cache(maxsize=8)
def _excitation_matrices(nact, nalpha, nbeta):
    """Spin-summed excitation matrices ``(E_tu)_{row,col} = <row| E_tu |col>``.

    Same determinant ordering and bit convention (alpha bits low, beta bits
    high) as ``fci._determinants`` / ``rdm.make_rdm*``.  Cached per active
    space; the returned stack is read-only."""
    dets = tuple(_determinants(int(nact), (int(nalpha), int(nbeta))))
    index = {det: i for i, det in enumerate(dets)}
    ndet = len(dets)
    stack = np.zeros((nact, nact, ndet, ndet))
    for col, det in enumerate(dets):
        for off in (0, nact):
            for u in range(nact):
                ann = annihilate(det, u + off)
                if ann is None:
                    continue
                det_u, phase_u = ann
                for t in range(nact):
                    cre = create(det_u, t + off)
                    if cre is None:
                        continue
                    det_tu, phase_t = cre
                    row = index.get(det_tu)
                    if row is not None:
                        stack[t, u, row, col] += phase_u * phase_t
    stack.setflags(write=False)
    return dets, stack


def _active_operator_matrix(f, g, stack):
    """Dense determinant-basis matrix of ``sum f_tu E_tu + 1/2 sum g (E E - d E)``."""
    ham = np.tensordot(f, stack, axes=([0, 1], [0, 1]))
    inter = np.tensordot(g, stack, axes=([2, 3], [0, 1]))  # (t, u, x, b)
    ham += 0.5 * np.tensordot(stack, inter, axes=([0, 1, 3], [0, 1, 2]))
    ham -= 0.5 * np.tensordot(np.einsum("tuuw->tw", g), stack, axes=([0, 1], [0, 1]))
    return ham


def _apply_active_operator(f, g, stack, wmat):
    """Apply the active operator of ``(f, g)`` to the vector behind ``wmat``.

    ``wmat[t, u, :] = E_tu |c>`` is precomputed once per CI vector; each call
    is then pure dense tensor algebra (no determinant loops)."""
    sigma = np.tensordot(f, wmat, axes=([0, 1], [0, 1]))
    x = np.tensordot(g, wmat, axes=([2, 3], [0, 1]))  # (t, u, ndet)
    sigma += 0.5 * np.einsum("tuab,tub->a", stack, x, optimize=True)
    sigma -= 0.5 * np.einsum("tw,twa->a", np.einsum("tuuw->tw", g), wmat, optimize=True)
    return sigma


# ----------------------------------------------------------------- main assembly
def analytic_orbital_hessian(h1e, eri, ncore, nact, active_nelec, pairs,
                             weights, roots, ci_coeffs, *,
                             degeneracy_tol=_DEGENERACY_TOL):
    """Exact orbital-rotation Hessian of the (SA-)CASSCF energy.

    Parameters mirror the ``casscf._optimize`` context: MO-basis ``h1e`` /
    chemist ``eri`` at the current orbitals, sequential ``ncore``/``nact``
    spaces, the ``pairs`` list of ``casscf._nonredundant_pairs`` (ordering is
    preserved in the returned matrix), the state-average ``weights`` over the
    ``roots`` columns of ``ci_coeffs`` (determinant-basis CI vectors from the
    solve at these orbitals; for a state-specific run one weight of 1.0).

    Returns the dense ``(n_par, n_par)`` Hessian matching the symmetrized
    ``casscf._fd_orbital_hessian`` limit.  Raises for degenerate roots with
    genuine coupling (non-smooth objective) and for active spaces beyond the
    dense-spectrum memory guard.
    """
    h1e = np.asarray(h1e, dtype=float)
    eri = np.asarray(eri, dtype=float)
    nbf = h1e.shape[0]
    if h1e.shape != (nbf, nbf) or eri.shape != (nbf, nbf, nbf, nbf):
        raise ValueError("analytic hessian needs square h1e and matching 4-index eri")
    ncore = int(ncore)
    nact = int(nact)
    nalpha, nbeta = int(active_nelec[0]), int(active_nelec[1])
    pairs = [(int(p), int(q)) for (p, q) in pairs]
    npar = len(pairs)
    weights = np.asarray(weights, dtype=float)
    roots = [int(r) for r in roots]
    if weights.shape != (len(roots),):
        raise ValueError("weights and roots must have matching lengths")

    ndet = comb(nact, nalpha) * comb(nact, nbeta)
    stack_bytes = 8 * nact * nact * ndet * ndet
    if stack_bytes > _MAX_STACK_BYTES:
        raise ValueError(
            f"analytic CASSCF Hessian needs a {nact}^2 x {ndet}^2 excitation stack "
            f"(~{stack_bytes / 1024**3:.1f} GiB); this active space is beyond the "
            "dense-spectrum path -- use [casscf] hessian = fd")
    ci = np.asarray(ci_coeffs, dtype=float)
    if ci.ndim != 2 or ci.shape[0] != ndet or max(roots) >= ci.shape[1]:
        raise ValueError("ci_coeffs must be (ndet, nroot) covering the averaged roots")

    dets, stack = _excitation_matrices(nact, nalpha, nbeta)

    # --- state-averaged RDMs (same construction as casscf._solve_active)
    from oqp.library.casscf import _full_rdms  # deferred: casscf lazily imports us

    D = np.zeros((nbf, nbf))
    G = np.zeros((nbf, nbf, nbf, nbf))
    for w, r in zip(weights, roots):
        gamma = make_rdm1_spatial(ci[:, r], dets, nact)
        Gamma = make_rdm2_spatial(ci[:, r], dets, nact)
        Dr, Gr = _full_rdms(gamma, Gamma, ncore, nact, nbf)
        D += w * Dr
        G += w * Gr

    # --- part 1: fixed-CI orbital Hessian + per-pair derivative integrals
    B = np.empty((npar, npar))
    f_der = np.empty((npar, nact, nact))
    g_der = np.empty((npar, nact, nact, nact, nact))
    for l, (p, q) in enumerate(pairs):
        t_l = _one_index_derivative_h(h1e, p, q)
        T_l = _one_index_derivative_g(eri, p, q)
        B[:, l] = _pair_differences(_z_matrix(D, G, t_l, T_l), pairs)
        f_der[l], g_der[l] = _fold_active(t_l, T_l, ncore, nact)
    hess = 0.5 * (B + B.T)

    # --- part 2: CI relaxation over the complete active spectrum
    f0, g0 = _fold_active(h1e, eri, ncore, nact)
    hact = _active_operator_matrix(f0, g0, stack)
    eps, vecs = _symmetric_eigh(hact)
    ovl = (ci[:, roots].T @ vecs) ** 2  # (n_avg, ndet) averaged-root overlaps

    for i_avg, (w_i, r) in enumerate(zip(weights, roots)):
        c = np.array(ci[:, r], dtype=float, copy=True)
        nrm = float(np.linalg.norm(c))
        if nrm <= 0.0:
            raise ValueError("CI vector for an averaged root has zero norm")
        c /= nrm
        e_i = float(c @ (hact @ c))
        wmat = np.tensordot(stack, c, axes=([3], [0]))  # (nact, nact, ndet)
        amp = np.empty((npar, ndet))
        for k in range(npar):
            amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

        factors = np.zeros(ndet)
        for j in range(ndet):
            if ovl[i_avg, j] > 0.5:
                continue  # this eigenstate IS the averaged root I
            m = int(np.argmax(ovl[:, j]))
            if ovl[m, j] > 0.5:
                # coupling within the averaged set: (w_I - w_J), split over the
                # two ordered visits -> exact cancellation for equal weights
                coef = 0.5 * (w_i - float(weights[m]))
            else:
                coef = float(w_i)
            if coef == 0.0:
                continue
            denom = e_i - float(eps[j])
            if abs(denom) < degeneracy_tol:
                if float(np.max(np.abs(amp[:, j]))) * abs(coef) < _COUPLING_NOISE:
                    continue  # symmetry/spin-forbidden partner: exact zero coupling
                raise ValueError(
                    "analytic CASSCF Hessian: root degeneracy with non-zero orbital "
                    f"coupling (|E_I - E_J| = {abs(denom):.2e}); the objective is not "
                    "smooth here and no orbital Hessian is defined")
            factors[j] = 2.0 * coef / denom
        nz = factors != 0.0
        if np.any(nz):
            hess += (amp[:, nz] * factors[nz]) @ amp[:, nz].T

    return hess


def make_hessian_provider(hcore_ao, eri_ao, ncore, nact, active_nelec, pairs,
                          weights, roots, settings=None):
    """Bind the AO integrals and active-space plan into ``hess_fn(C, ci_coeffs)``.

    The returned callable matches the pluggable-Hessian slot of
    ``casscf._floored_newton_loop`` / ``casscf_convergers._ah_inner``: it
    transforms the AO integrals at the given orbitals and assembles the
    analytic Hessian from the CI vectors of the *current* solve (which the
    macroiteration loops already hold), so it never triggers a CI solve of its
    own."""
    if settings is not None and (
        tuple(getattr(settings, "active_orbital_indices", ()) or ())
        or tuple(getattr(settings, "core_orbital_indices", ()) or ())
    ):
        raise ValueError(
            "casscf.hessian=analytic supports sequential active spaces only "
            "(unset cas.active_orbital_indices / cas.core_orbital_indices)")
    hcore_ao = np.asarray(hcore_ao, dtype=float)
    eri_ao = np.asarray(eri_ao, dtype=float)
    ncore = int(ncore)
    nact = int(nact)
    active_nelec = (int(active_nelec[0]), int(active_nelec[1]))
    pairs = [(int(p), int(q)) for (p, q) in pairs]
    weights = np.asarray(weights, dtype=float)
    roots = [int(r) for r in roots]

    def hess_fn(C, ci_coeffs):
        h1e, eri = _transform_integrals(hcore_ao, eri_ao, C)
        return analytic_orbital_hessian(
            h1e, eri, ncore, nact, active_nelec, pairs, weights, roots, ci_coeffs)

    return hess_fn
