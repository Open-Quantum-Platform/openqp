"""Analytic nuclear gradient for strongly contracted NEVPT2 (SC-NEVPT2).

This is a genuine analytic derivative of the SC-NEVPT2 total energy: every
term is a closed-form expression, there is no finite difference anywhere on the
path, and the cost is independent of the number of nuclei.

Scope
-----
* **Strongly contracted** NEVPT2 only (``[pt2] h0=dyall, contraction=strong``).
  The uncontracted Dyall-H0 path (``contraction=none``) and the partially
  contracted variant (which OpenQP does not implement at all) are refused
  rather than silently given this formula -- they are different first-order
  spaces and therefore different derivatives.
* **State specific**, on a converged **state-specific CASSCF** reference.
* Closed-shell singlet, Hartree-Fock integrals -- the same conditions the
  SC-NEVPT2 energy already enforces.

Theory
------
SC-NEVPT2 is not variational, so the energy derivative needs the response of
the reference.  The Lagrangian is

    L = E_CAS + E^(2)
        + sum_{pq in nonredundant} z_pq  g^CAS_pq(C, c)
        + sum_K            y_K  r_K(C, c),      r_K = 2[(H_act - E_CAS) c]_K

with ``g^CAS`` the CASSCF orbital-rotation gradient and ``r`` the CI residual;
both vanish at the reference, so ``L = E`` there, and the multipliers are fixed
by making ``L`` stationary against every wavefunction parameter.  That is the
ordinary coupled CASSCF response (Z-vector) system

    [ H_oo  H_oc ] [ z ]     [ dE2/dkappa ]
    [ H_oc^T H_cc ] [ y ]  = -[ dE2/dc     ],

whose blocks are the CASSCF electronic Hessian, reused from
:mod:`oqp.library.casscf_hessian` rather than re-derived.

Three structural facts make the rest exact and compact.

**1. The active gauge is free.**  SC-NEVPT2 is invariant under an active-active
orbital rotation when the CI vector rotates with it (verified numerically to
1e-12).  Intra-active rotations therefore need no multiplier: stationarity
there follows from CI stationarity, and it shows up as a symmetric
active-active block of the final Lagrangian.

**2. The semicanonical basis is not free.**  The strongly contracted
denominators are diagonal in the Dyall zeroth-order Hamiltonian, which is only
correct when the inactive and virtual generalized-Fock blocks are diagonal.  The
energy is consequently NOT invariant under core-core or virtual-virtual
rotations, and the response of the semicanonicalizing rotation is a real
contribution -- on a CAS(4,4) probe it is ~4% of the total derivative, far
above any tolerance one would accept.

The invariant form of the functional carries the full inactive and virtual Fock
blocks and coincides with the implemented one when they are diagonal.  Its
off-diagonal derivative follows from invariance alone, with no per-subspace
derivation: for ``i != j`` inside one block an infinitesimal rotation gives
``dF_ij = (eps_i - eps_j) kappa_ij`` (twice, for ``ij`` and ``ji``) and leaves
every diagonal element stationary, so

    0 = dE2_inv/dkappa_ij = dE2_diag/dkappa_ij + 2 Lambda_ij (eps_i - eps_j)

    =>  Lambda_ij = (F_ij - F_ji) / (eps_i - eps_j),

``F`` being the orbital Lagrangian of the diagonal functional.  The correctness
of this closure is checked at run time: with ``Lambda`` folded in, the
intra-inactive and intra-virtual blocks of the total Lagrangian must be
symmetric, and they are (~1e-14).

**3. The multipliers reduce to one-index transformations.**  The derivative of
``sum z_pq g^CAS_pq`` with respect to the integrals is the one-index
transformation of the CASSCF densities by the antisymmetric ``Z``, and the
derivative of the CI term is the ``y``/``c`` transition density.  Both enter the
same effective one- and two-particle densities as the CASSCF and second-order
pieces, so a single contraction with the derivative integrals
(``nevpt2_gradient.F90``) finishes the gradient.

The energy-weighted density is then the generalized Fock of the total relaxed
densities.  Because ``L`` is stationary against every orbital rotation, that
matrix must come out symmetric -- which is the strongest single self-check
available, and it is enforced below.

Denominator edge cases
----------------------
``Lambda`` divides by ``eps_i - eps_j`` inside one block, and the two ways that
can go wrong are genuinely different.

An EXACT (symmetry) degeneracy -- any linear molecule's pi pair, any atom -- is
a free gauge: a rotation inside the degenerate set preserves the Fock block's
diagonality, so the energy does not depend on it and the numerator vanishes with
the denominator.  The 0/0 limit is zero.  That is verified rather than assumed
(the numerator is measured against the Lagrangian's own scale) and it is what
makes the method usable on symmetric molecules at all.

A NEAR degeneracy is not rescued by any invariance: the ratio is a large finite
number formed from a difference of nearly equal quantities.  That case is
refused with the offending pair named, rather than returning a
plausible-looking number.
"""
from __future__ import annotations

from math import comb

import numpy as np

import oqp
from oqp.library.casscf import (
    _full_rdms,
    _generalized_fock,
    _kappa_matrix,
    _log,
    _nonredundant_pairs,
)
from oqp.library.casscf_hessian import (
    _active_hamiltonian,
    _apply_active_operator,
    _excitation_matrices,
    _fold_active,
    _one_index_derivative_g,
    _one_index_derivative_h,
    _pair_differences,
    _z_matrix,
)
from oqp.library.caspt2_dyall import (
    _caspt2_options,
    _freeze_core,
    _pt2_memory,
    _reference_roots,
    _run_casscf_reference,
    _semicanonicalize,
    _standard_core_count,
)
from oqp.library.fci import (
    _unpack_lower_triangle,
    check_ao_eri_budget,
    contiguous_active_space,
    settings_from_casci_config,
)
from oqp.library.nevpt2_adjoint import sc_nevpt2_energy_adjoints
from oqp.library.nevpt2_sc import make_rdms
from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

#: Default ceiling on the basis size.  The gradient holds several dense nbf^4
#: tensors (the second-order two-particle density, its AO transform and the
#: derivative-integral buffer), so the guard is tighter than the energy's.
#: Overridable through ``[pt2] max_memory`` / ``[cas] max_memory``.
DEFAULT_MAX_NBF = 64

#: Accepted CASSCF orbital-rotation gradient norm at the reference.  The
#: Lagrangian assumes ``g^CAS = 0``; the gradient error from a non-stationary
#: point is first order in this norm.
_STATIONARITY_FLOOR = 1.0e-4
#: Multiple of the run's own declared threshold above which the point is
#: reported as loosely converged (but still differentiated).
_STATIONARITY_WARN = 1.0e2
#: Largest accepted asymmetry of the total orbital Lagrangian, in Hartree.  A
#: violation means a multiplier is wrong, not that the answer is imprecise.
_LAGRANGIAN_ASYMMETRY_LIMIT = 1.0e-6
#: Inactive/virtual orbital-energy gaps below this make the semicanonical
#: multiplier ill-conditioned; between this and _EXACT_DEGENERACY_TOL the run is
#: refused rather than guessed.
_DEGENERACY_TOL = 1.0e-6
#: Gaps below this are a symmetry degeneracy, i.e. exactly degenerate at working
#: precision.  A rotation inside such a set preserves the Fock block's
#: diagonality, so it is a free gauge and the multiplier's 0/0 limit is zero.
_EXACT_DEGENERACY_TOL = 1.0e-10
#: Largest residual energy dependence, relative to the Lagrangian's own scale,
#: still accepted as "this rotation is a free gauge".  Measured at 1.1e-20 on
#: LiH/STO-3G's exactly degenerate 2p-pi virtuals, so this is orders of
#: magnitude of headroom over the case it is meant to admit.
_FREE_GAUGE_COUPLING = 1.0e-10
#: Accepted relative residual of the coupled response solve.
_RESPONSE_RESIDUAL_TOL = 1.0e-8
#: Dyall denominators below this are reported as a likely intruder state.  The
#: same threshold ``caspt2_dyall`` uses for the uncontracted path, so the two
#: report the same physics with the same number.
_INTRUDER_DENOMINATOR = 0.1


# --------------------------------------------------------------- small helpers
def _sym8(g):
    """Project a four-index tensor onto the eight-fold symmetry of a real ERI.

    Only this component of a raw derivative ever contributes, because the
    integrals it multiplies carry that symmetry; the derivative-ERI driver also
    requires it of the density it is handed.
    """
    g = 0.5 * (g + g.transpose(1, 0, 2, 3))
    g = 0.5 * (g + g.transpose(0, 1, 3, 2))
    g = 0.5 * (g + g.transpose(2, 3, 0, 1))
    return g


def _one_index_transform_1(K, D):
    """``sum_m (K_pm D_mq + K_qm D_pm)`` -- the commutator ``[K, D]``."""
    return K @ D - D @ K


def _one_index_transform_2(K, G):
    """One-index transformation of a four-index density by antisymmetric ``K``.

    The derivative of ``sum_pq K_pq F_pq(D, G)`` with respect to the two-electron
    integrals: ``K`` is applied to each slot of ``G`` in turn, the four-index
    analogue of :func:`_one_index_transform_1`.
    """
    out = np.einsum("pm,mqrs->pqrs", K, G, optimize=True)
    out += np.einsum("qm,pmrs->pqrs", K, G, optimize=True)
    out += np.einsum("rm,pqms->pqrs", K, G, optimize=True)
    out += np.einsum("sm,pqrm->pqrs", K, G, optimize=True)
    return out


def _orbital_lagrangian(hbar, gbar, h1e, eri):
    """Generalized Fock of the effective densities behind ``(hbar, gbar)``.

    ``hbar``/``gbar`` are raw derivatives with respect to independent integral
    elements.  Only their symmetric projections can contribute (the integrals
    they multiply are symmetric), and after that projection the pair
    ``(P, Gamma) = (sym(hbar), 2 sym8(gbar))`` is in exactly the CASSCF
    convention ``E = sum P h + 1/2 sum Gamma g``, so the same
    :func:`oqp.library.casscf._generalized_fock` applies and

        dE/dkappa_pq = 2 (X_qp - X_pq)

    over the pair list, the convention ``casscf_gfock_grad`` uses.
    """
    P = 0.5 * (hbar + hbar.T)
    G = 2.0 * _sym8(gbar)
    return P, G, _generalized_fock(P, G, h1e, eri)


# ------------------------------------------------------------ CI-space algebra
def _excitation_vectors(stack, vec):
    """``wmat[t, u, :] = E_tu |vec>`` for the determinant-basis vector ``vec``."""
    return np.tensordot(stack, np.asarray(vec, dtype=float), axes=([3], [0]))


def _apply_dm_bar(bar, order, stack, W1, W2):
    """``sum_I bar[I] E_{p1q1} ... E_{pnqn} |c>`` for one density-matrix adjoint.

    ``bar`` is ``dE2/d dm_n`` in the ``dm_n[p1,q1,...] = <0|E_p1q1 ... |0>``
    convention.  The chain is split in half so the cost stays a pair of GEMMs
    over the composite index ``a = (p, q)`` instead of an explicit product of
    ``n`` sparse operators: ``W2[a, b] = E_a E_b |c>`` is built once and reused
    by every order.
    """
    ndet = stack.shape[2]
    A = stack.shape[0] * stack.shape[1]
    S = stack.reshape(A, ndet, ndet)
    if order == 1:
        return np.einsum("A,Aa->a", bar.reshape(A), W1, optimize=True)
    if order == 2:
        Y = bar.reshape(A, A) @ W1
    elif order == 3:
        Y = bar.reshape(A, A * A) @ W2
    elif order == 4:
        Y2 = (bar.reshape(A * A, A * A) @ W2).reshape(A, A, ndet)
        Y = np.einsum("Bab,ABb->Aa", S, Y2, optimize=True)
    else:  # pragma: no cover - defensive
        raise ValueError(f"unsupported density-matrix order {order}")
    return np.einsum("Aab,Ab->a", S, Y, optimize=True)


def _ci_bar_from_dm_bars(dmbars, stack, ci):
    """``dE2/dc`` from the adjoints of ``dm1..dm4``.

    ``dm_n[I] = <c| O_I |c>`` with ``O_I = E_p1q1 ... E_pnqn``, so

        d dm_n[I] / d c = O_I |c> + O_I^dagger |c>,

    and ``O_I^dagger`` is the same product read backwards with every pair
    transposed -- which is exactly ``bar`` with all its axes reversed.
    """
    ndet = stack.shape[2]
    W1 = _excitation_vectors(stack, ci).reshape(-1, ndet)
    A = W1.shape[0]
    S = stack.reshape(A, ndet, ndet)
    W2 = np.einsum("Aab,Bb->ABa", S, W1, optimize=True).reshape(A * A, ndet)
    out = np.zeros(ndet)
    for order, bar in enumerate(dmbars, start=1):
        if not np.any(bar):
            continue
        out += _apply_dm_bar(bar, order, stack, W1, W2)
        out += _apply_dm_bar(bar.T, order, stack, W1, W2)
    return out


def _transition_densities(stack, bra, ket):
    """Spin-free transition 1-/2-RDM between two determinant-basis vectors.

    Conventions match :func:`oqp.library.casscf_hessian._active_operator_matrix`,
    whose operator is ``sum f_tu E_tu + 1/2 sum g_tuvw (E_tu E_vw - d_uv E_tw)``,
    so these are exactly the densities that reproduce ``<bra|H_act|ket>``:

        gamma_tu   = <bra| E_tu |ket>
        Gamma_tuvw = <bra| E_tu E_vw |ket> - delta_uv <bra| E_tw |ket>

    ``E_tu`` is not self-adjoint (``E_tu^dagger = E_ut``), so the product is
    evaluated as ``(E_ut |bra>) . (E_vw |ket>)``.  Both densities are projected
    onto the symmetry of the integrals they meet -- the folded one-electron
    matrix is symmetric and the active ERI block is eight-fold symmetric -- which
    changes neither contraction.
    """
    nact = stack.shape[0]
    wket = _excitation_vectors(stack, ket)
    wbra = _excitation_vectors(stack, bra)
    gamma_raw = np.einsum("a,tua->tu", np.asarray(bra, dtype=float), wket,
                          optimize=True)
    two = np.einsum("uta,vwa->tuvw", wbra, wket, optimize=True)
    for u in range(nact):
        two[:, u, u, :] -= gamma_raw
    return 0.5 * (gamma_raw + gamma_raw.T), _sym8(two)


# ------------------------------------------------------- semicanonical closure
def _semicanonical_multipliers(F0, eps, ncore, nact, nbf):
    """``Lambda`` on the inactive-inactive and virtual-virtual Fock blocks.

    See the module docstring: invariance of the exact (non-diagonal-denominator)
    functional under an intra-block rotation fixes these multipliers without any
    per-subspace derivation,

        Lambda_ij = (F_ij - F_ji) / (eps_i - eps_j).

    Degeneracy splits into two genuinely different cases, and they get different
    answers rather than one blanket refusal.

    **Exactly degenerate** (a symmetry degeneracy: any linear molecule's pi
    pair, any atom).  A rotation inside such a set preserves the diagonality of
    the Fock block, so the energy is invariant along it -- the direction is a
    free gauge, and the NUMERATOR vanishes with the denominator.  The limit of a
    genuine 0/0 along a direction the energy does not depend on is zero.  That
    is not assumed: the numerator is measured, and only a numerator small on the
    scale of the Lagrangian itself is accepted as the free-gauge case.  Verified
    on LiH/STO-3G, whose 2p-pi virtuals are exactly degenerate -- the numerator
    comes out at 1.1e-20 and the resulting gradient reproduces five-point
    central differences to 3.4e-11 Eh/Bohr.

    **Near degenerate** but not exactly so.  Here the ratio is a large finite
    number computed as a difference of nearly equal quantities, i.e.
    ill-conditioned, and no invariance rescues it.  This is refused with the
    offending pair named, rather than returning a plausible-looking number.
    """
    nocc = ncore + nact
    Lam = np.zeros((nbf, nbf))
    smallest_gap = float("inf")
    scale = max(float(np.max(np.abs(F0))) if F0.size else 1.0, 1.0)
    degenerate_pairs = 0
    for lo, hi in ((0, ncore), (nocc, nbf)):
        for i in range(lo, hi):
            for j in range(lo, hi):
                if i == j:
                    continue
                gap = float(eps[i] - eps[j])
                smallest_gap = min(smallest_gap, abs(gap))
                coupling = float(F0[i, j] - F0[j, i])
                if abs(gap) < _EXACT_DEGENERACY_TOL:
                    if abs(coupling) > _FREE_GAUGE_COUPLING * scale:
                        raise ValueError(
                            "Analytic SC-NEVPT2 gradient refused: "
                            f"semicanonical orbitals {i} and {j} are degenerate "
                            f"to {abs(gap):.3e} Hartree, but the rotation "
                            "between them is NOT a free gauge -- the energy "
                            f"still depends on it at {abs(coupling):.3e} "
                            "Hartree. The semicanonical multiplier is a 0/0 "
                            "with a non-vanishing numerator and is genuinely "
                            "undetermined here. Break the degeneracy, or use a "
                            "numerical gradient.")
                    degenerate_pairs += 1
                    continue                       # free gauge: Lambda_ij = 0
                if abs(gap) < _DEGENERACY_TOL:
                    raise ValueError(
                        "Analytic SC-NEVPT2 gradient refused: semicanonical "
                        f"orbitals {i} and {j} are separated by only "
                        f"{abs(gap):.3e} Hartree (limit "
                        f"{_DEGENERACY_TOL:.1e}). The strongly contracted "
                        "denominators are diagonal only in the semicanonical "
                        "basis, and the response of the rotation that produces "
                        "that basis is ill-conditioned at a near-degeneracy "
                        "this small -- neither zero (the pair is not exactly "
                        "degenerate, so the rotation is not a free gauge) nor "
                        "reliably computable. Break the degeneracy "
                        "(symmetry-distinct geometry or a different active "
                        "space), or use a numerical gradient.")
                Lam[i, j] = coupling / gap
    return Lam, smallest_gap, degenerate_pairs


# ------------------------------------------------------------- response solve
def _response_blocks(h1e, eri, D, G, ncore, nact, pairs, stack, dets, ci):
    """CASSCF electronic-Hessian blocks at the reference.

    Returns ``(H_oo, sigma, hact, e_ci)`` where ``sigma[k]`` is
    ``(dH_act/dkappa_k) |c>`` in the determinant basis; the orbital-CI coupling
    is ``2 sigma`` and the CI-CI block is ``2 (H_act - E)``, both projected onto
    the complement of ``c`` by the caller.
    """
    npar = len(pairs)
    B = np.empty((npar, npar))
    f_der = np.empty((npar, nact, nact))
    g_der = np.empty((npar, nact, nact, nact, nact))
    for l, (p, q) in enumerate(pairs):
        t_l = _one_index_derivative_h(h1e, p, q)
        T_l = _one_index_derivative_g(eri, p, q)
        B[:, l] = _pair_differences(_z_matrix(D, G, t_l, T_l), pairs)
        f_der[l], g_der[l] = _fold_active(t_l, T_l, ncore, nact)
    H_oo = 0.5 * (B + B.T)

    wmat = _excitation_vectors(stack, ci)
    sigma = np.empty((npar, stack.shape[2]))
    for l in range(npar):
        sigma[l] = _apply_active_operator(f_der[l], g_der[l], stack, wmat)

    f0, g0 = _fold_active(h1e, eri, ncore, nact)
    hact = _active_hamiltonian(f0, g0, stack, dets, nact)
    hact = 0.5 * (hact + hact.T)
    e_ci = float(ci @ (hact @ ci))
    return H_oo, sigma, hact, e_ci


def _solve_response(H_oo, sigma, hact, e_ci, ci, rhs_orb, rhs_ci):
    """Solve the coupled orbital/CI response system in the complement of ``c``.

    The CI parameters are the components of ``dc`` orthogonal to ``c`` (the
    norm direction is fixed by the normalization constraint and carries no
    information), so the system is built in an orthonormal basis of that
    complement.  That also removes the exact null vector the unprojected
    CI-CI block would otherwise carry.
    """
    npar = H_oo.shape[0]
    _u, _s, vt = np.linalg.svd(ci.reshape(1, -1))
    Q = vt[1:].T                                   # (ndet, ndet-1), orthonormal
    H_cc = 2.0 * (Q.T @ (hact @ Q) - e_ci * (Q.T @ Q))
    H_oc = 2.0 * (sigma @ Q)
    A = np.zeros((npar + Q.shape[1], npar + Q.shape[1]))
    A[:npar, :npar] = H_oo
    A[:npar, npar:] = H_oc
    A[npar:, :npar] = H_oc.T
    A[npar:, npar:] = H_cc
    b = -np.concatenate([np.asarray(rhs_orb, dtype=float), Q.T @ rhs_ci])
    try:
        x = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        x = np.linalg.lstsq(A, b, rcond=None)[0]
    resid = float(np.linalg.norm(A @ x - b))
    scale = max(float(np.linalg.norm(b)), 1.0)
    if resid > _RESPONSE_RESIDUAL_TOL * scale:
        raise ValueError(
            "Analytic SC-NEVPT2 gradient refused: the coupled CASSCF response "
            f"system did not solve (relative residual {resid / scale:.3e}). "
            "The reference is at (or extremely near) a point where the CASSCF "
            "electronic Hessian is singular -- typically a root crossing. Move "
            "off the crossing or use a numerical gradient.")
    z_orb = x[:npar]
    y_ci = Q @ x[npar:]
    return z_orb, y_ci, resid / scale


# ------------------------------------------------------------ relaxed densities
def _relaxed_densities(hbar, gbar, D_cas, G_cas, K, gamma_ci, Gamma_ci,
                       ncore, nact, nbf):
    """Total relaxed one- and two-particle densities in the MO basis.

    Four additive pieces, in the shared convention
    ``E = sum D h + 1/2 sum Gamma g``:

    * the CASSCF reference densities;
    * the second-order effective densities, i.e. the symmetric projections of
      the SC-NEVPT2 integral adjoints;
    * the orbital multiplier, entering as the one-index transformation of the
      CASSCF densities by the antisymmetric ``K``;
    * the CI multiplier, entering as its transition densities embedded exactly
      the way :func:`oqp.library.casscf._full_rdms` embeds the CI densities --
      the all-active block is the transition 2-RDM and the inactive cross terms
      come from linearizing the separable part.
    """
    P2 = 0.5 * (hbar + hbar.T)
    G2 = 2.0 * _sym8(gbar)

    Pz = _one_index_transform_1(K, D_cas)
    Gz = _one_index_transform_2(K, G_cas)

    active = list(range(ncore, ncore + nact))
    Dci = np.zeros((nbf, nbf))
    Dci[np.ix_(active, active)] = gamma_ci
    Gci = (np.einsum("pq,rs->pqrs", Dci, D_cas, optimize=True)
           + np.einsum("pq,rs->pqrs", D_cas, Dci, optimize=True)
           - 0.5 * np.einsum("ps,rq->pqrs", Dci, D_cas, optimize=True)
           - 0.5 * np.einsum("ps,rq->pqrs", D_cas, Dci, optimize=True))
    Gci[np.ix_(active, active, active, active)] = Gamma_ci

    D_tot = D_cas + P2 + Pz + Dci
    G_tot = G_cas + G2 + Gz + Gci
    return D_tot, G_tot


# --------------------------------------------------------------------- driver
def _validate(mol, options, settings):
    """Refuse every configuration this derivative is not the derivative of."""
    if options.contraction != "strong" or options.h0 != "dyall":
        raise ValueError(
            "The analytic NEVPT2 gradient implements the STRONGLY CONTRACTED "
            "variant ([pt2] h0=dyall, contraction=strong). The uncontracted "
            "Dyall-H0 path ([pt2] contraction=none) expands a different "
            "first-order space and has a different derivative; the partially "
            "contracted variant is not implemented in OpenQP at all. Set "
            "[pt2] h0=dyall and contraction=strong, or use a numerical "
            "gradient.")
    if options.variant != "caspt2":
        raise ValueError(
            "The analytic SC-NEVPT2 gradient is state specific. Multistate "
            "(MS/XMS) strong contraction is not implemented for the energy "
            "either; drop [pt2] variant.")
    if options.family != "caspt2":
        raise ValueError(
            "The analytic SC-NEVPT2 gradient does not apply to the QDPT2 "
            "family (mrmp2/mcqdpt2/xmcqdpt2), which uses a diagonal-Fock "
            "zeroth order.")
    if options.reference != "casscf":
        raise ValueError(
            "The analytic SC-NEVPT2 gradient requires a CASSCF reference "
            "([pt2] reference=casscf). With reference=casci the orbitals are "
            "not stationary for the CASSCF energy, so the orbital multiplier "
            "solves a different response equation than the one implemented "
            "here.")
    unapplied = [name for name, value in (
        ("level_shift", options.level_shift),
        ("imaginary_shift", options.imaginary_shift),
        ("ipea_shift", options.ipea_shift),
        ("edshft", options.edshft)) if value]
    if unapplied:
        raise ValueError(
            "[pt2] %s does not reach the strongly contracted denominators, so "
            "the SC-NEVPT2 energy already refuses it; the gradient does too."
            % ", ".join(unapplied))
    if settings.integral_backend != "native":
        raise ValueError("SC-NEVPT2 gradient integral_backend must be native")
    if mol.config["input"]["functional"]:
        raise ValueError(
            "The analytic SC-NEVPT2 gradient requires Hartree-Fock integrals; "
            "unset input.functional.")
    if mol.config["scf"]["type"] != "rhf":
        raise ValueError(
            "The analytic SC-NEVPT2 gradient supports a closed-shell RHF "
            "reference only.")
    if int(mol.data["nelec_A"]) != int(mol.data["nelec_B"]):
        raise ValueError(
            "The analytic SC-NEVPT2 gradient supports closed-shell singlets "
            "only.")
    method = str(mol.config["input"]["method"]).strip().lower().replace("_", "-")
    if method in {"sa-casscf", "sacasscf"} or str(
            mol.config.get("state_average", {}).get("enabled", False)
    ).strip().lower() in ("true", "1", "yes", "on"):
        raise ValueError(
            "The analytic SC-NEVPT2 gradient is built on a STATE-SPECIFIC "
            "CASSCF reference. A state-averaged reference optimizes "
            "sum_I w_I E_I, whose orbital multiplier solves a different "
            "response equation.")


#: Dense nbf^4 tensors held simultaneously on the gradient path: the MO
#: integrals, the second-order integral adjoint before and after the Fock
#: chain rule, the CASSCF and total two-particle densities, and the AO-basis
#: density handed to the derivative-integral driver.
_LIVE_NBF4_TENSORS = 6
#: Dense nact^8 / nact^6 tensors: the four-particle RDM and its adjoint, the
#: three-particle RDM and its adjoint, and the two eri-folded intermediates
#: with their adjoints.
_LIVE_NACT8_TENSORS = 2


def _check_size(nbf, nact=None, ndet=None, budget_mib=None, budget_label=""):
    """Refuse a system whose dense intermediates exceed the configured budget.

    The energy already guards its own ``nact^8`` four-particle RDM; the gradient
    holds that RDM's ADJOINT alongside it, plus several ``nbf^4`` tensors the
    energy never forms, so it needs its own accounting rather than inheriting
    the energy's.
    """
    if nbf > DEFAULT_MAX_NBF:
        raise ValueError(
            f"The analytic SC-NEVPT2 gradient holds dense nbf^4 tensors and is "
            f"guarded at nbf <= {DEFAULT_MAX_NBF}; this basis has nbf={nbf}. "
            "Reduce the basis, or use a numerical gradient.")
    if nact is None or budget_mib is None:
        return
    nbf4 = 8 * int(nbf) ** 4
    nact8 = 8 * int(nact) ** 8
    nact6 = 8 * int(nact) ** 6
    work = 8 * int(ndet or 0) * int(nact) ** 4
    live = (_LIVE_NBF4_TENSORS * nbf4
            + _LIVE_NACT8_TENSORS * nact8
            + 4 * nact6
            + 2 * work)
    cap = max(1, int(budget_mib)) * 1024 ** 2
    if live > cap:
        raise ValueError(
            "The analytic SC-NEVPT2 gradient needs ~%.2f GiB of dense "
            "intermediates (%d x nbf^4 at nbf=%d, the four-particle RDM and "
            "its adjoint at nact=%d, and the CI-derivative workspace), above "
            "the %s max_memory budget of %d MiB. Reduce the basis or the "
            "active space, raise %s max_memory, or use a numerical gradient."
            % (live / 1024 ** 3, _LIVE_NBF4_TENSORS, nbf, nact,
               budget_label or "[cas]", int(budget_mib),
               budget_label or "[cas]"))


def sc_nevpt2_gradient_route(mol):
    """Which nuclear-gradient route a configured PT2 calculation should take.

    Returns ``("analytic", "")`` when the calculation is exactly the one this
    module differentiates, and ``("numerical", reason)`` otherwise.  With
    ``[pt2] gradient=analytic`` an out-of-scope calculation raises instead of
    falling back, so a run that asked for the analytic derivative never
    silently receives central differences.
    """
    options = _caspt2_options(mol.config)
    if options.gradient == "numerical":
        return "numerical", "[pt2] gradient=numerical"
    settings = settings_from_casci_config(mol.config)
    try:
        _validate(mol, options, settings)
        nbf = int(mol.data.get_basis()["nbf"])
        nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
        _ncore, nact, active_nelec, _plan = contiguous_active_space(
            nbf, nelec, settings, "SC-NEVPT2 gradient")
        mem, mem_label = _pt2_memory(options, settings)
        _check_size(nbf, nact,
                    comb(nact, active_nelec[0]) * comb(nact, active_nelec[1]),
                    mem, mem_label)
        if _gradient_backend() is None:
            raise ValueError(
                "the compiled OpenQP library has no nevpt2_gradient entry "
                "point; rebuild liboqp to use analytic SC-NEVPT2 gradients")
    except ValueError as exc:
        if options.gradient == "analytic":
            raise
        return "numerical", str(exc)
    return "analytic", ""


def _gradient_backend():
    """liboqp ``(lib, ffi)`` exposing the derivative-integral contraction."""
    from oqp.library.fci import _lib_backend
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "nevpt2_gradient"):
        return None
    return lib, ffi


_NEVPT2G_MESSAGES = {
    -1: "the basis size disagreed with the native handle",
    -2: "the native driver could not allocate its working arrays",
    -3: "the relaxed two-particle density was not eight-fold permutation "
        "symmetric",
}


def sc_nevpt2_analytic_gradient(mol, ref_energy=None):
    """Analytic nuclear gradient of the SC-NEVPT2 total energy.

    Returns a ``(1, natom, 3)`` array and leaves the same gradient in the native
    buffer ``mol.get_grad()`` reads, so every gradient-driven runtype (``grad``,
    ``optimize``, ``ts``, ``mep``, ``irc``) consumes it through the ordinary
    path.  ``mol.energies`` is set to the SC-NEVPT2 total energy of the
    differentiated root, so a caller needs no separate energy pass.
    """
    options = _caspt2_options(mol.config)
    settings = settings_from_casci_config(mol.config)
    _validate(mol, options, settings)

    backend = _gradient_backend()
    if backend is None:
        raise ValueError(
            "The compiled OpenQP library has no nevpt2_gradient entry point; "
            "rebuild liboqp to use analytic SC-NEVPT2 gradients.")
    lib, ffi = backend

    nbf = int(mol.data.get_basis()["nbf"])
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    ncore, nact, active_nelec, _plan = contiguous_active_space(
        nbf, nelec, settings, "SC-NEVPT2 gradient")
    mem, mem_label = _pt2_memory(options, settings)
    _check_size(nbf, nact,
                comb(nact, active_nelec[0]) * comb(nact, active_nelec[1]),
                mem, mem_label)
    roots = _reference_roots(options)
    root = int(roots[0])
    weights = np.array([1.0])

    # ---- reference: the same state-specific CASSCF the energy driver runs
    _run_casscf_reference(mol, ref_energy, roots, weights)
    check_ao_eri_budget(nbf, mem, mem_label)
    oqp.fci_ao_integrals(mol)
    hcore_ao = _unpack_lower_triangle(
        np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
    eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
        (nbf, nbf, nbf, nbf), order="F")
    coeff = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape(
        (nbf, nbf)).T
    enuc = float(mol.mol_energy.nenergy)

    coeff_sc, h1e, eri, coeffs, energies, dets, eps, D_sa = _semicanonicalize(
        hcore_ao, eri_ao, coeff, ncore, nact, active_nelec, enuc, settings,
        roots, weights)
    target = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float)
    mol.data["OQP::VEC_MO_A"][...] = np.ascontiguousarray(
        coeff_sc.T.reshape(target.shape))
    ci = np.asarray(coeffs[:, root], dtype=float)
    ci = ci / np.linalg.norm(ci)
    e_casci = float(energies[root])

    # ---- CASSCF densities and the stationarity precondition
    gamma = make_rdm1_spatial(ci, dets, nact)
    Gamma_act = make_rdm2_spatial(ci, dets, nact)
    D_cas, G_cas = _full_rdms(gamma, Gamma_act, ncore, nact, nbf)
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    F_cas = _generalized_fock(D_cas, G_cas, h1e, eri)
    g_orb = np.array([2.0 * (F_cas[q, p] - F_cas[p, q]) for (p, q) in pairs])
    gnorm = float(np.linalg.norm(g_orb)) if g_orb.size else 0.0
    declared = float(mol.config.get("casscf", {}).get("gradient_norm_tol", 1e-6)
                     or 1e-6)
    limit = max(_STATIONARITY_FLOOR, _STATIONARITY_WARN * declared)
    if not np.isfinite(gnorm) or gnorm > limit:
        raise ValueError(
            "Analytic SC-NEVPT2 gradient refused: the CASSCF orbital-rotation "
            f"gradient at the reference is {gnorm:.3e}, above the acceptance "
            f"limit {limit:.3e}. The Lagrangian assumes a stationary CASSCF "
            "reference and its error is first order in this norm. Converge "
            "the orbital optimization ([casscf] gradient_norm_tol / "
            "max_macro_iterations).")

    # ---- PT2 frozen core: the second order is evaluated in the reduced space
    if options.frozen < 0:
        nfrozen = max(0, min(_standard_core_count(mol), ncore))
    else:
        nfrozen = int(options.frozen)
        if nfrozen > ncore:
            raise ValueError(
                f"[pt2] frozen={nfrozen} exceeds the {ncore} inactive "
                "orbital(s) available.")
        nfrozen = max(0, nfrozen)
    if nfrozen:
        h1e_f, eri_f, eps_f, _D_f, ncore_f, _nbf_f, _enuc_f = _freeze_core(
            h1e, eri, eps, D_sa, ncore, nbf, enuc, nfrozen)
    else:
        h1e_f, eri_f, eps_f, ncore_f = h1e, eri, eps, ncore

    # ---- second order: energy and its exact derivatives
    dms = make_rdms(ci, nact, active_nelec, upto=4)
    (e2, comp, hbar_f, gbar_f, epsbar_f, dmbars,
     min_denominator) = sc_nevpt2_energy_adjoints(
        h1e_f, eri_f, eps_f, ncore_f, nact, active_nelec, ci, dms=dms)
    e_total = e_casci + e2

    # unfold the frozen-core dressing back onto the full orbital space
    hbar0 = np.zeros((nbf, nbf))
    gbar0 = np.zeros((nbf,) * 4)
    epsbar = np.zeros(nbf)
    f = nfrozen
    hbar0[f:, f:] = hbar_f
    gbar0[f:, f:, f:, f:] = gbar_f
    epsbar[f:] = epsbar_f
    if f:
        fidx = np.arange(f)
        gbar0[f:, f:, fidx, fidx] += 2.0 * hbar_f[:, :, None]
        gbar0[f:, fidx, fidx, f:] -= hbar_f[:, None, :]

    # ---- semicanonical closure (module docstring, point 2)
    _P0, _G0, F0 = _orbital_lagrangian(hbar0, gbar0, h1e, eri)
    Lam, smallest_gap, degenerate_pairs = _semicanonical_multipliers(
        F0, eps, ncore, nact, nbf)
    Fbar = np.diag(epsbar) + Lam
    hbar = hbar0 + Fbar
    gbar = (gbar0
            + np.einsum("pq,rs->pqrs", Fbar, D_cas, optimize=True)
            - 0.5 * np.einsum("pq,rs->prsq", Fbar, D_cas, optimize=True))
    Dbar = (np.einsum("pq,pqrs->rs", Fbar, eri, optimize=True)
            - 0.5 * np.einsum("pq,prsq->rs", Fbar, eri, optimize=True))

    # ---- right-hand sides of the response system
    P2, G2, F2 = _orbital_lagrangian(hbar, gbar, h1e, eri)
    for label, lo, hi in (("inactive", 0, ncore), ("virtual", ncore + nact, nbf)):
        if hi - lo < 2:
            continue
        blk = F2[lo:hi, lo:hi]
        asym = float(np.max(np.abs(blk - blk.T)))
        if asym > _LAGRANGIAN_ASYMMETRY_LIMIT:
            raise ValueError(
                f"Analytic SC-NEVPT2 gradient refused: the intra-{label} "
                f"Lagrangian is asymmetric by {asym:.3e} Hartree after the "
                "semicanonical closure. That closure is exact, so this "
                "indicates the reference is not semicanonical to working "
                "precision.")
    rhs_orb = np.array([2.0 * (F2[q, p] - F2[p, q]) for (p, q) in pairs])

    _dets_stack, stack = _excitation_matrices(nact, active_nelec[0],
                                              active_nelec[1])
    ci_bar = _ci_bar_from_dm_bars(dmbars, stack, ci)
    active = list(range(ncore, ncore + nact))
    dbar_act = Dbar[np.ix_(active, active)]
    dbar_act = dbar_act + dbar_act.T
    wmat = _excitation_vectors(stack, ci)
    ci_bar = ci_bar + np.einsum("tu,tua->a", dbar_act, wmat, optimize=True)
    ci_bar = ci_bar - float(ci @ ci_bar) * ci        # project onto c-complement

    # ---- coupled orbital/CI response
    H_oo, sigma, hact, e_ci = _response_blocks(
        h1e, eri, D_cas, G_cas, ncore, nact, pairs, stack, dets, ci)
    z_orb, y_ci, resid = _solve_response(
        H_oo, sigma, hact, e_ci, ci, rhs_orb, ci_bar)

    # ---- relaxed densities.  The CI Lagrangian term is sum_K y_K r_K with
    # r = 2 (H - E) c, i.e. 2 y^T H c, hence the factor two on the transition
    # densities of y.
    K = _kappa_matrix(z_orb, pairs, nbf)
    gamma_ci, Gamma_ci = _transition_densities(stack, 2.0 * y_ci, ci)
    D_tot, G_tot = _relaxed_densities(hbar, gbar, D_cas, G_cas, K,
                                      gamma_ci, Gamma_ci, ncore, nact, nbf)
    G_tot = _sym8(G_tot)

    # ---- energy-weighted density.  L is stationary against every orbital
    # rotation, so its generalized Fock must be symmetric; that is the single
    # strongest check available on the whole construction.
    X = _generalized_fock(D_tot, G_tot, h1e, eri)
    asym = float(np.max(np.abs(X - X.T)))
    scale = max(float(np.max(np.abs(X))), 1.0)
    if asym > _LAGRANGIAN_ASYMMETRY_LIMIT * scale:
        raise ValueError(
            "Analytic SC-NEVPT2 gradient refused: the total Lagrangian is "
            f"asymmetric by {asym:.3e} Hartree. A stationary Lagrangian has a "
            "symmetric generalized Fock, so a response multiplier did not "
            "solve; the gradient would be wrong rather than imprecise.")
    W = 0.5 * (X + X.T)

    # ---- MO -> AO and the derivative-integral contraction
    C = np.ascontiguousarray(coeff_sc)
    D_ao = C @ D_tot @ C.T
    W_ao = C @ W @ C.T
    G_ao = np.einsum("pqrs,ap,bq,cr,ds->abcd", G_tot, C, C, C, C, optimize=True)
    G_ao = _sym8(G_ao)

    d_c = np.ascontiguousarray(D_ao, dtype=np.float64)
    w_c = np.ascontiguousarray(W_ao, dtype=np.float64)
    g_c = np.ascontiguousarray(G_ao, dtype=np.float64)
    status = int(lib.nevpt2_gradient(
        mol.data._data, int(nbf),
        ffi.cast("double *", d_c.ctypes.data),
        ffi.cast("double *", w_c.ctypes.data),
        ffi.cast("double *", g_c.ctypes.data)))
    if status < 0:
        reason = _NEVPT2G_MESSAGES.get(
            status, f"the native driver declined (status {status})")
        raise ValueError(f"Analytic SC-NEVPT2 gradient failed: {reason}.")

    mol.energies = [e_total]
    mol.mol_energy.energy = e_total
    mol.data["OQP::CASPT2_ENERGIES"] = np.ascontiguousarray(
        [e_total], dtype=np.float64)
    mol.data["OQP::CASPT2_REFERENCE_ENERGIES"] = np.ascontiguousarray(
        [e_casci], dtype=np.float64)
    mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"] = np.ascontiguousarray(
        [e2], dtype=np.float64)

    _log(mol, "")
    _log(mol, "   ==============================================")
    _log(mol, "   PyOQP: analytic SC-NEVPT2 nuclear gradient")
    _log(mol, "   ==============================================")
    _log(mol, f"   {'differentiated root':<38}{root}")
    _log(mol, f"   {'E(CASSCF) reference':<38}{e_casci:>20.10f}")
    _log(mol, f"   {'E(2) strong contraction':<38}{e2:>20.10f}")
    _log(mol, f"   {'E(SC-NEVPT2) total':<38}{e_total:>20.10f}")
    _log(mol, f"   {'PT2 frozen-core orbitals':<38}{nfrozen:>20d}")
    _log(mol, f"   {'CASSCF orbital gradient |g_orb|':<38}{gnorm:>20.3e}")
    if np.isfinite(smallest_gap):
        _log(mol, f"   {'smallest semicanonical gap':<38}"
                  f"{smallest_gap:>20.3e}")
    if degenerate_pairs:
        _log(mol, f"   {'exactly degenerate pairs (free gauge)':<38}"
                  f"{degenerate_pairs:>20d}")
    intruder = ("" if min_denominator >= _INTRUDER_DENOMINATOR
                else "   <-- WARNING: likely intruder state")
    _log(mol, f"   {'smallest Dyall denominator':<38}"
              f"{min_denominator:>20.3e}{intruder}")
    if min_denominator < _INTRUDER_DENOMINATOR:
        _log(mol, "   NOTE: the strongly contracted energy admits no level "
                  "shift, so an intruder")
        _log(mol, "         cannot be regularized here; the derivative goes "
                  "as the denominator")
        _log(mol, "         SQUARED, so it degrades faster than the energy "
                  "does.")
    _log(mol, f"   {'response solve relative residual':<38}{resid:>20.3e}")
    _log(mol, f"   {'Lagrangian asymmetry max|X-X^T|':<38}{asym:>20.3e}")
    _log(mol, f"   {'orbital multiplier norm |z|':<38}"
              f"{float(np.linalg.norm(z_orb)):>20.3e}")
    _log(mol, f"   {'CI multiplier norm |y|':<38}"
              f"{float(np.linalg.norm(y_ci)):>20.3e}")
    _log(mol, "   Dyall subspace energies:")
    for key in ("Sr", "Si", "Sijrs", "Sijr", "Srsi", "Srs", "Sij", "Sir"):
        _log(mol, f"      {key:7s} {comp[key]:18.10f}")
    _log(mol, "   ==============================================")

    natom = int(mol.data["natom"])
    grad = np.asarray(mol.get_grad(), dtype=float).reshape((natom, 3))
    return np.array([grad.copy()]).reshape((1, natom, 3))
