"""Analytic CASPT2 / MS-CASPT2 / XMS-CASPT2 nuclear gradient.

This is the derivative of the *reported* PT2 total energy with respect to the
nuclear coordinates.  Nothing here is a finite difference: every term is an
explicit derivative of the functional :mod:`oqp.library.caspt2_dyall` evaluates,
and every non-variational parameter carries its own multiplier equation.

Why the determinant formulation makes this tractable
----------------------------------------------------
OpenQP's CASPT2 is *uncontracted*: the first-order interacting space ``Q`` is the
set of external determinants, an orthonormal basis whose labels do not depend on
the geometry.  There is no internally-contracted overlap metric, no linear
dependence removal, and no perturber normalization -- the three things that make
an internally-contracted CASPT2 gradient hard.  What remains is

    E = E_ref + E2 ,   E_ref = <Psi_0|H|Psi_0> + E_nuc ,
    E2 = V . T ,       V = (P_Q H |Psi_0>) ,   T = -G(D)^{-1} V ,
    D = P_Q H0 P_Q - E0 ,                      E0 = <Psi_0|H0|Psi_0>

with ``G`` the (regularized) denominator function of §2 below.

Zeroth-order Hamiltonian in orbital-rotation invariant form
-----------------------------------------------------------
The energy path builds ``H0 = diag(eps)`` with ``eps_p = f_pp`` and ``f`` the
closed+active Fock ``h + J(D_sa) - 1/2 K(D_sa)``.  That representation is tied to
the semicanonical basis.  This module uses instead the *block-diagonal Fock
operator*

    H0^inv = sum_{p,q in the same orbital block} f_pq E_pq

which is numerically IDENTICAL at the semicanonical orbitals (where the
within-block off-diagonals of ``f`` vanish, checked here and refused otherwise)
but is invariant under rotations inside the inactive, active and virtual blocks.
Those rotations are therefore genuinely redundant: they need no constraint, no
multiplier, and no semicanonical-condition derivative.  The residual asymmetry of
the final generalized Fock is reported as a diagnostic of exactly this.

The IPEA shift biases the active *diagonal* of H0 in a particular basis and so
breaks that invariance; it is refused rather than differentiated approximately.

Amplitudes, denominators and shifts
------------------------------------
``E2 = V.T`` is not stationary in ``T`` once a level shift is applied, so the
amplitudes carry an explicit multiplier ``lambda``.  For the single-state case
the multiplier equation ``V + G(D) lambda = 0`` gives ``lambda = T`` exactly, and
the Lagrangian collapses to the shifted Hylleraas functional; in the multistate
case ``lambda_I = -u_I G(D_I)^{-1} (sum_J u_J V_J)`` and is solved for.

The derivative of ``G(D)`` is the derivative of a *matrix function*, given by the
Daleckii-Krein divided-difference formula

    lambda . dG(D) . T = Tr[Omega dD] ,
    Omega~_ij = G^[1](w_i, w_j) lambda~_i T~_j     (symmetrized)

with ``G^[1](a,b) = (G(a)-G(b))/(a-b)``.  The option validator in
``caspt2_dyall._caspt2_options`` admits only three regularization families
(real shift; real+imaginary; GAMESS ISA ``edshft``), and all three have

    G(w) = d + e/d ,  d = w + level_shift ,  e in {0, imaginary_shift^2, edshft}

so ``G^[1](a,b) = 1 - e/(d_a d_b)`` and ``Omega`` is exactly rank <= 2 (rank <= 4
in the multistate case).  No ``n_ext x n_ext`` object is ever formed and the
shift derivative is exact rather than neglected.

Reference relaxation
--------------------
The reference orbitals and CI vectors are not variational for the PT2 energy, so
both carry multipliers:

* CI -- one projected linear solve per reference root,
  ``(H_model - E_r) v_r = -(1 - c_r c_r^T) g_r``;
* orbitals -- a Z-vector against whatever actually fixes the reference orbitals:
  ``[pt2] reference=casscf`` uses the (state-averaged) CASSCF orbital-rotation
  gradient, ``reference=casci`` with RHF orbitals uses ``F^RHF_pq = 0``.

Both are expressed in the REFERENCE orbitals, not the semicanonical ones: for
``reference=casci`` the RHF occupied space is the first ``nocc`` reference
orbitals, and there is no way to say that in a basis whose active block has been
rotated.  Which rotations are parameters follows from the same fact and is set
out in :func:`_orbital_pair_sets`; briefly, the PT2 energy depends on the
reference orbitals only through the SPANS of the orbital blocks, so within-block
rotations are redundant for the energy -- but for ``reference=casci`` the
active-active ones are still constrained, because the active block straddles the
RHF occupied/virtual boundary and such a rotation changes ``D^RHF``.

Both Jacobians are closed-form.  With ``kappa_A`` an antisymmetric generator and
the rotation-derivative rule ``dM = [M, kappa]`` on every tensor index,

    dF^RHF/dt_A = [F^RHF, kappa_A] - G([D^RHF, kappa_A])
    dF^gen/dt_A = [F^gen,  kappa_A] - F^gen(dD_A, dGamma_A)

which for the RHF case reduces to the textbook CPHF ``A`` matrix.  The systems
are assembled and solved densely: this module is validation-grade, in the same
sense the CASPT2 energy path is (the external determinant space is enumerated in
full), so the orbital-pair and CI dimensions are small by construction.

XMS
---
For ``xms-caspt2``/``xmcqdpt2`` the references are pre-rotated by the eigenvectors
``R`` of the model-space state-averaged Fock.  ``R`` depends on the geometry, so
its response is included -- eliminated directly with first-order eigenvector
perturbation theory rather than carried as a constrained parameter, which turns
the whole rotation response into one extra weight matrix on the Fock operator
plus one extra term in the reference-vector gradient (:func:`_xms_elimination`).
A near-degenerate model-space Fock makes ``dR/dx`` undefined and is refused.

The *effective-Hamiltonian* eigenvectors need no response (Hellmann-Feynman for
an eigenvalue).  Which published state the derivative belongs to is pinned by
requiring this module's own reconstruction of the spectrum to reproduce the
reported ``OQP::CASPT2_ENERGIES`` (:func:`_check_reported_energies`), and the gap
to the neighbouring root is checked before differentiating, because a degenerate
root has no well-defined mixing vector.

Scope
-----
See ``_gate_variant``: single-state CASPT2/MRMP2, MCQDPT2 and XMS-CASPT2/XMCQDPT2
are differentiated; the multi-set MS-CASPT2 construction, Dyall H0 (NEVPT2),
strong contraction, a nonzero IPEA shift and imported (non-RHF, non-CASSCF)
orbitals are refused with a specific message.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from oqp.library.caspt2_dyall import (
    CASPT2Setup,
    _build_operators,
    _caspt2_setup,
    _effective_fock,
    _embed_reference,
    _freeze_core,
    _pt2_frozen_count,
    _pt2_memory,
    _semicanonicalize,
    _xms_rotation,
)
from oqp.library.casscf import _full_rdms, _nonredundant_pairs
from oqp.library.fci import _symmetric_eigh, _transform_integrals
from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial

#: Largest within-block off-diagonal of the closed+active Fock still accepted as
#: "semicanonical".  Above this the coded ``diag(eps)`` H0 and the invariant
#: block-diagonal H0 this module differentiates are different operators, so the
#: gradient would not belong to the reported energy.
SEMICANONICAL_TOL = 1.0e-8

#: Largest accepted CASSCF orbital-rotation gradient norm for a
#: ``reference=casscf`` run.  The reference constraint is ``g_orb = 0``; a run
#: that did not reach it has no well-defined orbital response.
CASSCF_STATIONARITY_TOL = 1.0e-5

#: Largest accepted RHF Fock off-diagonal for a ``reference=casci`` run whose
#: orbitals come from the SCF.  The constraint is canonicality.
RHF_CANONICAL_TOL = 1.0e-6

#: Energy agreement required between this module's reconstruction and the value
#: the energy path reported.  A mismatch means the gradient does not belong to
#: the printed energy.
ENERGY_MATCH_TOL = 1.0e-9

#: Smallest accepted gap between the differentiated effective-Hamiltonian root
#: and its neighbour.  Below this the mixing vector is not well defined.
STATE_GAP_TOL = 1.0e-6

#: Guard on the dense orbital-response system.
MAX_ORBITAL_PAIRS = 4000


class CASPT2GradientNotImplemented(NotImplementedError):
    """A PT2 option combination outside the analytic gradient's scope."""


# --------------------------------------------------------------------- helpers
def _sym8(g):
    """Symmetrize a two-particle density over the eight permutations of an ERI.

    ``grd2`` contracts one representative per symmetry orbit with a multiplicity
    factor, so a density that carries fewer than eight symmetries would be
    mis-weighted for mixed quartets.  The contraction value is unchanged by this
    projection because ``(pq|rs)`` already has all eight.
    """
    g = np.asarray(g, dtype=float)
    out = g + g.transpose(1, 0, 2, 3)
    out = out + out.transpose(0, 1, 3, 2)
    out = out + out.transpose(2, 3, 0, 1)
    return out * 0.125


def _jk_operator(dens, eri):
    """``J(dens) - 1/2 K(dens)``, the two-electron part of the closed+active Fock.

    Same contraction as :func:`oqp.library.caspt2_dyall._effective_fock` without
    the core Hamiltonian, and symmetric under exchange of its argument with the
    matrix it is contracted against -- which is what lets the same routine serve
    both ``f`` itself and ``dL/dD_sa``.
    """
    j = np.einsum("rs,pqrs->pq", dens, eri, optimize=True)
    k = np.einsum("rs,prsq->pq", dens, eri, optimize=True)
    return j - 0.5 * k


def _fock_like_densities(weight_matrix, dens, eri):
    """Effective 1-/2-PDM of ``sum_pq Y_pq f_pq(dens)`` for a fixed ``dens``.

    ``f = h + J(dens) - 1/2 K(dens)``, so with the convention
    ``L = sum d h + 1/2 sum Gamma (pq|rs)``

        d      = Y
        Gamma  = 2 Y_ab dens_cd - Y_ad dens_bc

    (the second term is the exchange contraction ``-1/2 sum Y_pq dens_rs (pr|sq)``
    relabelled into chemist order).  ``Gamma`` is returned unsymmetrized; the
    caller symmetrizes once at the end.
    """
    y = np.asarray(weight_matrix, dtype=float)
    d = np.asarray(dens, dtype=float)
    gamma = 2.0 * np.einsum("ab,cd->abcd", y, d, optimize=True)
    gamma -= np.einsum("ad,bc->abcd", y, d, optimize=True)
    return y.copy(), gamma


def _block_project(mat, ncore, nact, norb):
    """Zero everything outside the inactive/active/virtual diagonal blocks."""
    out = np.zeros_like(mat)
    for lo, hi in ((0, ncore), (ncore, ncore + nact), (ncore + nact, norb)):
        if hi > lo:
            out[lo:hi, lo:hi] = mat[lo:hi, lo:hi]
    return out


def _invariant_h0(fock, ncore, nact, norb):
    """Block-diagonal part of the closed+active Fock (the invariant H0 matrix)."""
    sym = 0.5 * (fock + fock.T)
    return _block_project(sym, ncore, nact, norb)


def _generalized_fock(d1, gamma2, h1e, eri):
    """``F_pq = sum_r d_pr h_qr + sum_rst Gamma_prst (qr|st)``.

    The same object ``casscf_kernel.F90`` builds.  With ``L`` expressed through
    the MO integrals, ``dL/dkappa_mn = 2 (F_nm - F_mn)`` and the energy-weighted
    (Pulay) density is ``(F + F^T)/2``.
    """
    return d1 @ h1e + np.einsum("prst,qrst->pq", gamma2, eri, optimize=True)


def _orbital_gradient_matrix(fock):
    """``dL/dkappa`` as an antisymmetric matrix, ``2 (F^T - F)``."""
    return 2.0 * (fock.T - fock)


def _rotate_1(mat, u):
    """``u M u^T``: MO matrix from the semicanonical basis into the reference one."""
    return u @ mat @ u.T


def _rotate_2(tensor, u):
    """Four-index analogue of :func:`_rotate_1`."""
    out = np.einsum("pa,abcd->pbcd", u, tensor, optimize=True)
    out = np.einsum("qb,pbcd->pqcd", u, out, optimize=True)
    out = np.einsum("rc,pqcd->pqrd", u, out, optimize=True)
    out = np.einsum("sd,pqrd->pqrs", u, out, optimize=True)
    return out


def _kappa_generator(norb, m, n):
    """The antisymmetric unit generator ``E_mn - E_nm``."""
    k = np.zeros((norb, norb))
    k[m, n] = 1.0
    k[n, m] = -1.0
    return k


def _rotation_derivative_1(mat, kappa):
    """``dM = [M, kappa]`` -- the first-order change of an MO matrix under
    ``C -> C exp(kappa)``."""
    return mat @ kappa - kappa @ mat


def _rotation_derivative_2(tensor, kappa):
    """Four-index rotation derivative: ``-kappa`` applied to every index."""
    out = -np.einsum("pm,mqrs->pqrs", kappa, tensor, optimize=True)
    out -= np.einsum("qm,pmrs->pqrs", kappa, tensor, optimize=True)
    out -= np.einsum("rm,pqms->pqrs", kappa, tensor, optimize=True)
    out -= np.einsum("sm,pqrm->pqrs", kappa, tensor, optimize=True)
    return out


# ------------------------------------------------------------------ transition densities
def _rdm12(vec, dets, norb, want2=True):
    g1 = make_rdm1_spatial(vec, dets, norb)
    g2 = make_rdm2_spatial(vec, dets, norb) if want2 else None
    return g1, g2


def _transition_density(a, b, dets, norb, want2=True):
    """``(d, Gamma)`` with ``sum d h + 1/2 sum Gamma (pq|rs) = <a|H|b>`` exactly.

    Only the part of a transition density that is symmetric under the ERI
    permutations is ever contracted, and

        gamma^{a+b} - gamma^a - gamma^b = gamma^{ab} + gamma^{ba}
        Gamma^{a+b} - Gamma^a - Gamma^b = Gamma^{ab} + Gamma^{ba},
        Gamma^{ba}_pqrs = Gamma^{ab}_qpsr

    so half the polarization identity is exactly ``<a|H|b>``'s density and no
    transition-RDM kernel is needed -- three calls to the existing same-vector
    engines suffice.  ``a is b`` is special-cased to one call.
    """
    a = np.asarray(a, dtype=float)
    if a is b or b is None:
        return _rdm12(a, dets, norb, want2)
    b = np.asarray(b, dtype=float)
    ga, Ga = _rdm12(a, dets, norb, want2)
    gb, Gb = _rdm12(b, dets, norb, want2)
    gs, Gs = _rdm12(a + b, dets, norb, want2)
    d1 = 0.5 * (gs - ga - gb)
    d2 = None if not want2 else 0.5 * (Gs - Ga - Gb)
    return d1, d2


# ------------------------------------------------------------------------ shifts
def _shift_parameters(options):
    """``(s, e)`` of the denominator function ``G(w) = d + e/d``, ``d = w + s``.

    ``caspt2_dyall._caspt2_options`` makes ``edshft`` mutually exclusive with the
    real and imaginary shifts, so these three families are exhaustive:

    ==================================  ===============  ============
    options                             ``s``            ``e``
    ==================================  ===============  ============
    real level shift only               ``level_shift``  ``0``
    real + imaginary                    ``level_shift``  ``sigma^2``
    GAMESS ISA (``edshft``)             ``0``            ``edshft``
    ==================================  ===============  ============
    """
    s = float(options.level_shift)
    if options.edshft:
        return 0.0, float(options.edshft)
    if options.imaginary_shift:
        return s, float(options.imaginary_shift) ** 2
    return s, 0.0


def _denominator_apply(h0op, e0, vec_ext, func):
    """Apply a scalar function of ``D = H0_QQ - E0`` in its own eigenbasis.

    ``func`` receives the bare eigenvalues ``w`` and returns the multiplier to
    apply.  ``h0op`` owns the (root-independent) eigenbasis, exactly as
    :func:`caspt2_dyall._first_order` uses it.
    """
    w0, transforms = h0op.external_eigenbasis()
    w = w0 - e0
    tilde = np.empty_like(vec_ext)
    for idx, u in transforms:
        tilde[idx] = vec_ext[idx] if u is None else u.T @ vec_ext[idx]
    tilde = tilde * func(w)
    out = np.empty_like(tilde)
    for idx, u in transforms:
        out[idx] = tilde[idx] if u is None else u @ tilde[idx]
    return out


def _denominator_factor(w, options):
    """The ``1/G(w)`` factor, byte-for-byte the arithmetic of ``_first_order``."""
    d = w + options.level_shift
    if options.edshft:
        with np.errstate(divide="ignore", invalid="ignore"):
            d = d + options.edshft / d
    if options.imaginary_shift:
        return d / (d * d + options.imaginary_shift ** 2)
    return np.divide(1.0, d, out=np.zeros_like(d),
                     where=(np.abs(d) > 1.0e-300) & np.isfinite(d))


def _shifted_inverse(w, s):
    """``1/(w + s)``, the factor of ``T' = (D + s)^{-1} T``."""
    d = w + s
    return np.divide(1.0, d, out=np.zeros_like(d),
                     where=(np.abs(d) > 1.0e-300) & np.isfinite(d))


# ------------------------------------------------------------------------ state
@dataclass
class CASPT2GradientState:
    """Everything the derivative needs, in the semicanonical MO basis.

    Two orbital bases appear.  ``coeff_ref`` is the *reference* set -- the RHF
    canonical orbitals or the (SA-)CASSCF-optimized orbitals -- and is where the
    orbital constraint lives.  ``coeff_sc`` is the semicanonical set the PT2 is
    evaluated in; ``u_sc`` relates them, ``coeff_sc = coeff_ref @ u_sc``, and is
    block diagonal.  The two differ only by rotations inside the orbital blocks,
    which the invariant H0 makes redundant.
    """
    setup: CASPT2Setup
    mol: object

    # orbitals
    coeff_ref: np.ndarray
    coeff_sc: np.ndarray
    u_sc: np.ndarray

    # pre-freeze semicanonical quantities (full nbf space)
    h1e_full: np.ndarray
    eri_full: np.ndarray
    # the SAME integrals in the reference orbitals, where the orbital
    # constraint and its response live (see _reference_orbital_data)
    h1e_ref: np.ndarray
    eri_ref: np.ndarray
    ci_coeffs: np.ndarray
    ci_energies: np.ndarray
    act_dets: list
    eps_full: np.ndarray
    dsa_full: np.ndarray

    # frozen-core fold
    nfrozen: int
    h1e: np.ndarray
    eri: np.ndarray
    eps: np.ndarray
    dsa: np.ndarray
    ncore: int
    norb: int
    enuc: float

    # determinant-space operators
    full_dets: list
    det_index: dict
    hfull: np.ndarray
    h0op: object
    external: np.ndarray
    internal: np.ndarray
    h0mat: np.ndarray
    h0_offdiag: float

    # per-reference-state PT2 data (index = position in ``setup.roots``)
    refs: list = field(default_factory=list)       # model-space vectors (XMS rotated)
    amplitudes: list = field(default_factory=list)  # T_I, full determinant length
    e0s: list = field(default_factory=list)
    e2s: list = field(default_factory=list)
    min_denoms: list = field(default_factory=list)
    vvecs: list = field(default_factory=list)       # V_I on the external space

    # multistate
    xms_rotation: np.ndarray = None
    heff: np.ndarray = None
    mixing: np.ndarray = None
    total_energies: np.ndarray = None

    @property
    def nstate(self):
        return len(self.setup.roots)

    @property
    def nbf(self):
        return self.h1e_full.shape[0]


def _internal_indices(act_dets, ncore, nact, norb, det_index):
    """Full-determinant positions of the CAS (model-space) determinants."""
    core_mask = sum(1 << i for i in range(ncore))
    act_mask = (1 << nact) - 1
    out = np.empty(len(act_dets), dtype=np.int64)
    for k, ad in enumerate(act_dets):
        a_full = core_mask | ((ad & act_mask) << ncore)
        b_full = core_mask | ((ad >> nact) << ncore)
        out[k] = det_index[a_full | (b_full << norb)]
    return out


def _gate_variant(options, setup):
    """Refuse, with a specific reason, every combination outside scope."""
    if options.contraction != "none":
        raise CASPT2GradientNotImplemented(
            "[pt2] contraction='strong' (SC-NEVPT2) forms no external "
            "determinant space, so the amplitude Lagrangian this gradient is "
            "built on does not exist for it.  Use contraction=none.")
    if options.h0 != "fock":
        raise CASPT2GradientNotImplemented(
            "the analytic PT2 gradient is implemented for the CASPT2 (Fock) "
            "zeroth-order Hamiltonian; [pt2] h0=dyall (NEVPT2) has a different "
            "zeroth-order derivative and is out of scope.  Use "
            "runtype=grad with h0=fock, or keep the numerical gradient for "
            "NEVPT2.")
    if float(options.ipea_shift or 0.0) != 0.0:
        raise CASPT2GradientNotImplemented(
            "[pt2] ipea_shift biases the ACTIVE DIAGONAL of H0, which is not "
            "invariant under rotations inside the active block.  The analytic "
            "gradient uses the block-invariant form of H0 and would need the "
            "semicanonical constraint and its multipliers to carry an IPEA "
            "shift.  Set ipea_shift=0.0 for an analytic gradient.")
    if options.variant == "ms" and options.family == "caspt2":
        raise CASPT2GradientNotImplemented(
            "MS-CASPT2 uses the MULTI-SET construction (per-state orbitals from "
            "each root's own Fock, a per-state full-Fock-matrix H0 and "
            "inter-state Loewdin-minor rotations).  Its response is not "
            "implemented; use method=xms-caspt2, which is single-set and has an "
            "analytic gradient, or method=mcqdpt2 for the single-set multistate "
            "QDPT construction.")
    if options.reference not in {"casci", "casscf"}:
        raise CASPT2GradientNotImplemented(
            f"unknown [pt2] reference={options.reference}")
    if options.reference == "casci" and setup.orbital_source != "rhf":
        raise CASPT2GradientNotImplemented(
            "[pt2] reference=casci with [cas] orbital_source=%s: imported "
            "orbitals are not a differentiable function of the geometry, so "
            "there is no orbital response to solve for.  Use "
            "orbital_source=rhf or reference=casscf."
            % setup.orbital_source)


def _build_state(mol, ref_energy=None) -> CASPT2GradientState:
    """Rebuild the exact state the energy path produced, in one pass."""
    setup = _caspt2_setup(mol, ref_energy, run_reference=False)
    options = setup.options
    _gate_variant(options, setup)

    nbf = setup.nbf
    try:
        stash = mol.data["OQP::CASPT2_REFERENCE_MO"]
    except Exception:                                    # noqa: BLE001
        stash = None
    if stash is None:
        raise RuntimeError(
            "the analytic CASPT2 gradient needs the PT2 reference orbitals, "
            "which the energy path stores as OQP::CASPT2_REFERENCE_MO.  Run the "
            "CASPT2 energy in the same process before the gradient.")
    coeff_ref = np.ascontiguousarray(
        np.asarray(stash, dtype=float).reshape((nbf, nbf)))

    roots = setup.roots
    weights = setup.weights
    coeff_sc, h1e_full, eri_full, ci_coeffs, ci_energies, act_dets, eps_full, dsa_full = \
        _semicanonicalize(setup.hcore_ao, setup.eri_ao, coeff_ref, setup.ncore,
                          setup.nact, setup.active_nelec, setup.enuc,
                          setup.settings, roots, weights)

    ovl = _ao_overlap(mol, nbf)
    u_sc = coeff_ref.T @ ovl @ coeff_sc
    _check_block_rotation(u_sc, setup.ncore, setup.nact, nbf)

    nfrozen = _pt2_frozen_count(mol, options, setup.ncore)
    h1e, eri, eps, dsa, ncore, norb, enuc = (h1e_full, eri_full, eps_full,
                                             dsa_full, setup.ncore, nbf, setup.enuc)
    if nfrozen:
        h1e, eri, eps, dsa, ncore, norb, enuc = _freeze_core(
            h1e_full, eri_full, eps_full, dsa_full, setup.ncore, nbf,
            setup.enuc, nfrozen)

    nact = setup.nact
    active_occ = np.diag(dsa)[ncore:ncore + nact]
    full_dets, det_index, hfull, h0op, external = _build_operators(
        h1e, eri, eps, ncore, nact, setup.active_nelec, norb, options.max_det,
        options.h0, active_occ=active_occ, ipea=options.ipea_shift,
        max_memory=_pt2_memory(options, setup.settings)[0])
    internal = _internal_indices(act_dets, ncore, nact, norb, det_index)

    # Invariant H0 matrix and the semicanonical check that makes it equal to the
    # diag(eps) the energy path actually used.
    fock = _effective_fock(h1e, eri, dsa)
    h0mat = _invariant_h0(fock, ncore, nact, norb)
    offdiag = float(np.max(np.abs(h0mat - np.diag(np.diag(h0mat))))) if norb else 0.0
    if offdiag > SEMICANONICAL_TOL:
        raise RuntimeError(
            "the CASPT2 orbitals are not semicanonical to %.1e (largest "
            "within-block Fock off-diagonal %.3e).  The energy path used "
            "H0 = diag(eps); the analytic gradient differentiates the "
            "block-diagonal Fock operator, which only coincides with it at a "
            "semicanonical point.  Raise the semicanonicalization quality or "
            "keep the numerical gradient." % (SEMICANONICAL_TOL, offdiag))
    if abs(float(np.max(np.abs(np.diag(h0mat) - np.asarray(eps))))) > SEMICANONICAL_TOL:
        raise RuntimeError(
            "the invariant H0 diagonal does not reproduce the orbital energies "
            "the energy path used; refusing to differentiate a different "
            "operator.")

    h1e_ref, eri_ref = _transform_integrals(setup.hcore_ao, setup.eri_ao,
                                           coeff_ref)
    state = CASPT2GradientState(
        setup=setup, mol=mol, coeff_ref=coeff_ref, coeff_sc=coeff_sc, u_sc=u_sc,
        h1e_full=h1e_full, eri_full=eri_full, h1e_ref=h1e_ref, eri_ref=eri_ref,
        ci_coeffs=ci_coeffs,
        ci_energies=np.asarray(ci_energies), act_dets=act_dets,
        eps_full=eps_full, dsa_full=dsa_full, nfrozen=nfrozen, h1e=h1e, eri=eri,
        eps=eps, dsa=dsa, ncore=ncore, norb=norb, enuc=enuc,
        full_dets=full_dets, det_index=det_index, hfull=hfull, h0op=h0op,
        external=external, internal=internal, h0mat=h0mat, h0_offdiag=offdiag)

    _solve_first_order(state)
    return state


def _ao_overlap(mol, nbf):
    from oqp.library.fci import _unpack_lower_triangle
    return _unpack_lower_triangle(np.asarray(mol.data["OQP::SM"], dtype=float), nbf)


def _check_block_rotation(u_sc, ncore, nact, nbf):
    """The semicanonicalization must be a pure within-block rotation."""
    mask = np.ones((nbf, nbf), dtype=bool)
    for lo, hi in ((0, ncore), (ncore, ncore + nact), (ncore + nact, nbf)):
        mask[lo:hi, lo:hi] = False
    leak = float(np.max(np.abs(u_sc[mask]))) if mask.any() else 0.0
    if leak > 1.0e-8:
        raise RuntimeError(
            "the semicanonical orbitals are not a within-block rotation of the "
            "reference orbitals (largest inter-block element %.3e).  The "
            "analytic gradient relies on that to make the rotation redundant."
            % leak)


def _solve_first_order(state):
    """Reference vectors (XMS rotated if requested), amplitudes and E2."""
    setup = state.setup
    options = setup.options
    roots = setup.roots
    refs = [_embed_reference(state.ci_coeffs[:, r], state.act_dets, state.ncore,
                             setup.nact, state.norb, state.det_index)
            for r in roots]

    if options.variant == "xms":
        rot = _xms_rotation(state.h1e, state.eri, state.dsa, state.ci_coeffs,
                            state.act_dets, state.ncore, setup.nact, state.norb,
                            roots, state.det_index, state.hfull)
        state.xms_rotation = np.asarray(rot, dtype=float)
        refs = [sum(state.xms_rotation[k, i] * refs[k] for k in range(len(roots)))
                for i in range(len(roots))]
    state.refs = refs

    ext = state.external
    for vec in refs:
        e0 = state.h0op.expectation(vec)
        v_ext = (state.hfull @ vec)[ext]
        factor = _denominator_apply(
            state.h0op, e0, v_ext, lambda w: _denominator_factor(w, options))
        amp = np.zeros(len(vec))
        amp[ext] = -factor
        w0, _tr = state.h0op.external_eigenbasis()
        bare = w0 - e0
        state.e0s.append(float(e0))
        state.vvecs.append(v_ext)
        state.amplitudes.append(amp)
        state.e2s.append(float(v_ext @ amp[ext]))
        state.min_denoms.append(float(np.min(np.abs(bare))) if bare.size else float("inf"))

    nstate = len(roots)
    if options.variant == "caspt2":
        state.heff = np.array([[float(refs[0] @ state.hfull @ refs[0]) + state.enuc
                                + state.e2s[0]]])
        state.mixing = np.array([[1.0]])
        state.total_energies = np.array([state.heff[0, 0]])
        return

    hmodel = np.zeros((nstate, nstate))
    for i in range(nstate):
        hi = state.hfull @ refs[i]
        for j in range(nstate):
            hmodel[i, j] = refs[j] @ hi
        hmodel[i, i] += state.enuc
    heff = hmodel.copy()
    for i in range(nstate):
        for j in range(nstate):
            heff[i, j] += 0.5 * (state.vvecs[i] @ state.amplitudes[j][state.external]
                                 + state.vvecs[j] @ state.amplitudes[i][state.external])
    heff = 0.5 * (heff + heff.T)
    energies, mixing = _symmetric_eigh(heff)
    state.heff = heff
    state.mixing = mixing
    state.total_energies = np.asarray(energies)


# ------------------------------------------------------------- Lagrangian pieces
@dataclass
class LagrangianPieces:
    """Everything the derivative of one target state needs at fixed parameters.

    ``d1``/``d2`` are the *explicit* effective densities: ``dL/dh`` and
    ``2 dL/d(pq|rs)`` at fixed CI coefficients, amplitudes and orbitals.  The
    remaining fields are the parameter derivatives the response equations
    consume.
    """
    d1: np.ndarray
    d2: np.ndarray
    ymat: np.ndarray                 # H0 weight, block projected
    multipliers: list                # lambda_I, one per reference state
    omegas: list                     # Tr Omega_I
    psi_grad: list                   # dL/dPsi_I in the full determinant space
    mixing: np.ndarray               # u, the Heff mixing vector of the target


def _lagrangian_multipliers(state, u):
    """``lambda_I = -u_I G(D_I)^{-1} sum_J u_J V_J``.

    Reduces to ``lambda = T`` for a single state, which is the statement that
    the shifted Hylleraas functional is exactly stationary in the amplitudes.
    """
    ext = state.external
    options = state.setup.options
    w_ext = np.zeros(ext.size)
    for i, v in enumerate(state.vvecs):
        w_ext = w_ext + u[i] * v
    out = []
    for i, e0 in enumerate(state.e0s):
        factor = _denominator_apply(
            state.h0op, e0, w_ext, lambda w: _denominator_factor(w, options))
        lam = np.zeros(len(state.refs[i]))
        lam[ext] = -u[i] * factor
        out.append(lam)
    return out


def _omega_density(state, lam, amp, e0):
    """``(gamma^Omega, omega)`` for one reference state.

    ``Omega`` is the Daleckii-Krein weight of ``lambda . dG(D) . T``.  With
    ``G(w) = d + e/d`` the divided difference is ``G^[1](a,b) = 1 - e/(d_a d_b)``,
    so ``Omega = sym(lambda T^T) - e sym(lambda' T'^T)`` with
    ``X' = (D + s)^{-1} X`` -- rank at most four, never an ``n_ext x n_ext``
    matrix.
    """
    s, e = _shift_parameters(state.setup.options)
    dets, norb = state.full_dets, state.norb
    g1, _ = _transition_density(lam, amp, dets, norb, want2=False)
    omega = float(lam @ amp)
    if e:
        lam_p = np.zeros_like(lam)
        amp_p = np.zeros_like(amp)
        ext = state.external
        lam_p[ext] = _denominator_apply(state.h0op, e0, lam[ext],
                                        lambda w: _shifted_inverse(w, s))
        amp_p[ext] = _denominator_apply(state.h0op, e0, amp[ext],
                                        lambda w: _shifted_inverse(w, s))
        g1p, _ = _transition_density(lam_p, amp_p, dets, norb, want2=False)
        g1 = g1 - e * g1p
        omega -= e * float(lam_p @ amp_p)
    return g1, omega


def _unrelaxed_densities(state, u) -> LagrangianPieces:
    """``dL/dh`` and ``2 dL/d(pq|rs)`` at fixed orbitals, CI vectors and amplitudes.

    The Lagrangian of the target state (mixing vector ``u``) is

        L = <Psi_u|H|Psi_u> + E_nuc + <T_u|H|Psi_u>
          + sum_I [ lambda_I . G(D_I) T_I + <lambda_I|H|Psi_I> ]

    with ``Psi_u = sum_I u_I Psi_I`` and ``T_u = sum_I u_I T_I``.  Its value is
    the reported total energy and its parameter derivatives vanish by
    construction for the amplitudes.  Only ``H`` and ``H0`` carry an explicit
    integral dependence, and ``H0 = block-diag(f)`` propagates through
    ``f = h + J(D_sa) - 1/2 K(D_sa)``.
    """
    dets, norb = state.full_dets, state.norb
    nstate = state.nstate
    refs = state.refs

    psi_u = np.zeros(len(refs[0]))
    amp_u = np.zeros(len(refs[0]))
    for i in range(nstate):
        psi_u = psi_u + u[i] * refs[i]
        amp_u = amp_u + u[i] * state.amplitudes[i]

    d1, d2 = _transition_density(psi_u, psi_u, dets, norb)
    d1 = np.array(d1, dtype=float)
    d2 = np.array(d2, dtype=float)

    t1, t2 = _transition_density(amp_u, psi_u, dets, norb)
    d1 += t1
    d2 += t2

    lams = _lagrangian_multipliers(state, u)
    omegas = []
    ymat = np.zeros((norb, norb))
    for i in range(nstate):
        l1, l2 = _transition_density(lams[i], refs[i], dets, norb)
        d1 += l1
        d2 += l2
        gomega, omega = _omega_density(state, lams[i], state.amplitudes[i],
                                       state.e0s[i])
        gref, _ = _transition_density(refs[i], refs[i], dets, norb, want2=False)
        ymat += gomega - omega * gref
        omegas.append(omega)

    ymat = _block_project(0.5 * (ymat + ymat.T), state.ncore, state.setup.nact,
                          norb)
    y1, y2 = _fock_like_densities(ymat, state.dsa, state.eri)
    d1 += y1
    d2 += y2

    psi_grad = _reference_gradients(state, u, lams, amp_u, psi_u, omegas)
    return LagrangianPieces(d1=d1, d2=d2, ymat=ymat, multipliers=lams,
                            omegas=omegas, psi_grad=psi_grad, mixing=u)


def _reference_gradients(state, u, lams, amp_u, psi_u, omegas):
    """``dL/dPsi_I`` in the full determinant space, at fixed integrals.

    ``L`` sees ``Psi_I`` through ``<Psi_u|H|Psi_u>``, ``<T_u|H|Psi_u>``,
    ``<lambda_I|H|Psi_I>`` and ``E0_I = <Psi_I|H0|Psi_I>``; the last enters the
    denominator operator, whose derivative with respect to a scalar shift of
    ``D_I`` is ``-Tr Omega_I``.
    """
    h_psi_u = state.hfull @ psi_u
    h_amp_u = state.hfull @ amp_u
    out = []
    for i in range(state.nstate):
        g = 2.0 * u[i] * h_psi_u + u[i] * h_amp_u + state.hfull @ lams[i]
        g = g - 2.0 * omegas[i] * _apply_h0(state, state.refs[i])
        out.append(g)
    return out


def _apply_h0(state, vec):
    """``H0 |vec>`` in the full determinant space.

    ``H0`` is the one-electron operator whose MO matrix is the block-diagonal
    Fock; ``h0op`` carries the same operator in the representation the energy
    path chose (an exact diagonal for the Fock zeroth order).
    """
    if state.h0op.diagonal is not None:
        return state.h0op.diagonal * vec
    return state.h0op.dense @ vec


# ------------------------------------------------------------------- XMS response
def _xms_model_fock(state):
    """``(F_model, R, lambda)`` of the XMS pre-rotation.

    Rebuilds exactly what :func:`caspt2_dyall._xms_rotation` diagonalizes: the
    state-averaged Fock operator projected on the UNROTATED reference roots.
    """
    from oqp.library.fci import _build_dense_hamiltonian, _spin_orbital_integrals
    setup = state.setup
    fock = _effective_fock(state.h1e, state.eri, state.dsa)
    fspin, gzero = _spin_orbital_integrals(fock, np.zeros_like(state.eri))
    fock_det = _build_dense_hamiltonian(state.full_dets, state.det_index, fspin,
                                        gzero, 2 * state.norb, 0.0)
    base = [_embed_reference(state.ci_coeffs[:, r], state.act_dets, state.ncore,
                             setup.nact, state.norb, state.det_index)
            for r in setup.roots]
    n = len(base)
    fmodel = np.zeros((n, n))
    for i in range(n):
        fi = fock_det @ base[i]
        for j in range(n):
            fmodel[i, j] = base[j] @ fi
    fmodel = 0.5 * (fmodel + fmodel.T)
    lam, rot = np.linalg.eigh(fmodel)
    return fmodel, rot, lam, base, fock_det


def _xms_elimination(state, pieces):
    """Response of the XMS state rotation ``R``.

    ``R`` diagonalizes the model-space state-averaged Fock ``M``.  Eliminating
    ``dR`` with first-order eigenvector perturbation theory,

        sum_ki (dL/dR_ki) dR_ki = Tr[Z_R dM] ,
        Z_R = sym(R Phi^T R^T) ,  Phi_ji = (R^T dL/dR)_ji / (lam_i - lam_j)

    so the whole rotation response is one extra weight matrix on the Fock
    operator plus one extra term in the reference-vector gradient.  Returns
    ``(y_fock, base_grad, base_refs)``: the extra Fock weight, ``dL/dPsi_k`` for
    the UNROTATED roots, and those roots' model-space vectors.
    """
    setup = state.setup
    nstate = state.nstate
    if setup.options.variant != "xms":
        return np.zeros((state.norb, state.norb)), list(pieces.psi_grad), \
            list(state.refs)

    fmodel, rot, lam, base, fock_det = _xms_model_fock(state)
    if not np.allclose(rot, state.xms_rotation, atol=1.0e-10):
        # eigh is deterministic, but a degenerate model-space Fock leaves the
        # eigenvectors arbitrary and the two calls can disagree.  That also makes
        # dR/dx undefined, so refuse rather than differentiate a random basis.
        raise RuntimeError(
            "the XMS model-space Fock rotation could not be reproduced "
            "deterministically; its eigenvectors are not well defined (near "
            "degeneracy in the state-averaged Fock), so dR/dx does not exist.")
    gap = np.min(np.abs(lam[:, None] - lam[None, :]
                        + np.eye(nstate) * 1.0e30)) if nstate > 1 else np.inf
    if gap < STATE_GAP_TOL:
        raise RuntimeError(
            "the XMS model-space Fock has two eigenvalues within %.1e (%.3e); "
            "the state rotation and therefore its response are ill-conditioned."
            % (STATE_GAP_TOL, gap))

    gr = np.zeros((nstate, nstate))
    for k in range(nstate):
        for i in range(nstate):
            gr[k, i] = float(base[k] @ pieces.psi_grad[i])
    rtg = rot.T @ gr
    phi = np.zeros((nstate, nstate))
    for i in range(nstate):
        for j in range(nstate):
            if i != j:
                phi[j, i] = rtg[j, i] / (lam[i] - lam[j])
    zr = rot @ phi.T @ rot.T
    zr = 0.5 * (zr + zr.T)

    # dL/dPsi_k for the unrotated roots: through the rotation, and through M.
    base_grad = []
    for k in range(nstate):
        g = sum(rot[k, i] * pieces.psi_grad[i] for i in range(nstate))
        g = g + 2.0 * sum(zr[k, l] * (fock_det @ base[l]) for l in range(nstate))
        base_grad.append(g)

    # dTr[Z_R M]/df: the model-space Fock is <Psi_k| sum_pq f_pq E_pq |Psi_l>.
    yfock = np.zeros((state.norb, state.norb))
    for k in range(nstate):
        for l in range(nstate):
            if zr[k, l] == 0.0:
                continue
            g1, _ = _transition_density(base[k], base[l], state.full_dets,
                                        state.norb, want2=False)
            yfock += zr[k, l] * g1
    return 0.5 * (yfock + yfock.T), base_grad, base


#: Smallest accepted inactive orbital-energy gap across the frozen/correlated
#: boundary.  Below it the frozen SPACE is not determined by the closed+active
#: Fock, so neither is its response.
FROZEN_SPLIT_GAP_TOL = 1.0e-6


def _frozen_split_weight(state, gmat_sc):
    """Response of the frozen/correlated inactive split, as a Fock weight.

    The PT2 frozen core takes the ``nfrozen`` LOWEST eigenvectors of the
    closed+active Fock restricted to the inactive span.  That space is invariant
    under a within-inactive rotation at fixed ``f`` -- which is why the inactive
    block is redundant for everything else -- but ``f`` itself moves with the
    geometry, so the split moves with it, and ``L`` is not invariant under the
    resulting rotation.  With ``0 < nfrozen < ncore`` this is a real term:
    leaving it out costs about 1e-4 in the gradient and breaks rotational
    invariance.

    It is eliminated the same way the XMS rotation is.  In the semicanonical
    basis ``f`` is already diagonal on the inactive block, so first-order
    eigenvector perturbation theory gives the induced rotation directly,
    ``Theta_pq = df_pq / (lam_q - lam_p)``, and

        dL = sum_{p>q} (dL/dt_pq) Theta_pq = sum_pq W_pq df_pq ,
        W_pq = (1/2) (dL/dt_pq) / (lam_q - lam_p)

    so the whole response is one symmetric weight on ``f``, which then
    propagates through ``f = h + J(D_sa) - 1/2 K(D_sa)`` exactly like the H0
    weight does.  ``W`` is symmetric because ``dL/dt`` and ``lam_q - lam_p`` are
    both antisymmetric.

    ``gmat_sc`` is ``dL/dt`` in the SEMICANONICAL basis, taken before the orbital
    multipliers and without the CI response: a within-inactive rotation leaves
    the inactive span, hence the core dressing, hence the active Hamiltonian and
    its eigenvectors untouched, so the CI contributes nothing here and there is
    no circularity.

    Within-frozen and within-correlated pairs come out with a vanishing weight on
    their own -- ``L`` really is invariant under those -- so they need no special
    case.
    """
    setup = state.setup
    ncore = setup.ncore
    nbf = state.nbf
    w = np.zeros((nbf, nbf))
    if state.nfrozen <= 0 or state.nfrozen >= ncore or ncore < 2:
        return w
    lam = np.asarray(state.eps_full[:ncore], dtype=float)
    for p in range(ncore):
        for q in range(ncore):
            if p == q:
                continue
            gap = lam[q] - lam[p]
            if abs(gap) < FROZEN_SPLIT_GAP_TOL:
                if abs(float(gmat_sc[p, q])) > 1.0e-10:
                    raise RuntimeError(
                        "inactive orbitals %d and %d are degenerate to %.3e in "
                        "the closed+active Fock, and the PT2 frozen core splits "
                        "them: the frozen space is not determined, so neither is "
                        "its response.  Change [pt2] frozen so the split does "
                        "not fall inside a degenerate shell."
                        % (p, q, abs(gap)))
                continue
            w[p, q] = 0.5 * float(gmat_sc[p, q]) / gap
    return 0.5 * (w + w.T)


# -------------------------------------------------------------------- CI response
def _active_one_electron_operator(state, mat_act):
    """Dense active-space one-electron operator in the CAS determinant basis."""
    from oqp.library.fci import _build_dense_hamiltonian, _spin_orbital_integrals
    nact = state.setup.nact
    hspin, gspin = _spin_orbital_integrals(
        np.ascontiguousarray(mat_act, dtype=float),
        np.zeros((nact, nact, nact, nact)))
    index = {d: i for i, d in enumerate(state.act_dets)}
    return _build_dense_hamiltonian(state.act_dets, index, hspin, gspin,
                                    2 * nact, 0.0)


def _ci_gradients(state, pieces, yfock_total, base_grad, yfock_full=None):
    """``dL/dc_r`` in the active CI space, for each reference root.

    Two channels reach the CI vectors: the model-space projection of
    ``dL/dPsi_k``, and the state-averaged density ``D_sa`` that defines the Fock
    operator inside ``H0`` (and, for XMS, inside the model-space rotation, and
    for a split frozen core inside the inactive-block eigenproblem).  The second
    is ``dL/dD_sa = J(Y) - 1/2 K(Y)`` contracted with ``dD_sa/dc``, which for a
    weight-averaged CAS density is ``2 w_r <I| O |c_r>``.

    ``yfock_total`` lives in the frozen-core-reduced space (H0 and the XMS
    rotation are both built there); ``yfock_full`` is the frozen-split weight,
    which lives in the full space because the inactive-block eigenproblem does.
    The active orbitals are the same either way, so the two active blocks add.
    """
    setup = state.setup
    nact = setup.nact
    internal = state.internal
    smat = _jk_operator(yfock_total, state.eri)
    act = slice(state.ncore, state.ncore + nact)
    smat_act = 0.5 * (smat + smat.T)[act, act]
    if yfock_full is not None and np.any(yfock_full):
        full = _jk_operator(yfock_full, state.eri_full)
        act_full = slice(setup.ncore, setup.ncore + nact)
        smat_act = smat_act + 0.5 * (full + full.T)[act_full, act_full]
    sop = _active_one_electron_operator(state, smat_act)
    out = []
    for k in range(state.nstate):
        g = np.asarray(base_grad[k])[internal].copy()
        g = g + 2.0 * float(setup.weights[k]) * (sop @ state.ci_coeffs[:, setup.roots[k]])
        out.append(g)
    return out


def _ci_multipliers(state, ci_grads):
    """Solve ``(H_model - E_r) y_r = -(1 - c_r c_r^T) g_r`` with ``y_r . c_r = 0``."""
    hmodel = state.hfull[np.ix_(state.internal, state.internal)]
    hmodel = 0.5 * (hmodel + hmodel.T)
    ndet = hmodel.shape[0]
    out = []
    for k, g in enumerate(ci_grads):
        c = state.ci_coeffs[:, state.setup.roots[k]]
        e = float(c @ hmodel @ c)
        rhs = -(g - c * float(c @ g))
        amat = hmodel - e * np.eye(ndet)
        # Project the singular direction out of both sides and pin it with a
        # rank-one term, so the solve is nonsingular and returns the unique
        # solution orthogonal to c.
        proj = np.eye(ndet) - np.outer(c, c)
        amat = proj @ amat @ proj + np.outer(c, c)
        y = np.linalg.solve(amat, proj @ rhs)
        y = y - c * float(c @ y)
        out.append(y)
    return out


def _embed_ci(state, vec_act):
    """Lift an active CI vector into the full determinant space."""
    out = np.zeros(len(state.full_dets))
    out[state.internal] = vec_act
    return out


# ------------------------------------------------------------ frozen-core unfolding
def _frozen_core_density(nfrozen, nbf):
    d = np.zeros((nbf, nbf))
    for i in range(nfrozen):
        d[i, i] = 2.0
    return d


def _separable_gamma(dens):
    return (np.einsum("pq,rs->pqrs", dens, dens, optimize=True)
            - 0.5 * np.einsum("ps,rq->pqrs", dens, dens, optimize=True))


def _unfold_frozen(d_red, g_red, nfrozen, nbf, eri_full):
    """Lift the reduced-space effective densities back to the full MO space.

    ``_freeze_core`` is an exact algebraic fold: the frozen orbitals become an
    extra closed shell whose mean field dresses ``h`` and whose energy goes into
    ``E_nuc``.  Undoing it is therefore the same closed-shell algebra --
    ``d += D^fc``, ``Gamma += Gamma^sep(D^fc)`` for the constant, plus the
    ``Y``-against-``D^fc`` Fock term for the dressing of ``h``.
    """
    if nfrozen == 0:
        return np.array(d_red, dtype=float), np.array(g_red, dtype=float)
    f = int(nfrozen)
    d_full = np.zeros((nbf, nbf))
    d_full[f:, f:] = d_red
    g_full = np.zeros((nbf,) * 4)
    g_full[f:, f:, f:, f:] = g_red
    dfc = _frozen_core_density(f, nbf)
    _y1, y2 = _fock_like_densities(d_full, dfc, eri_full)
    g_full += y2
    g_full += _separable_gamma(dfc)
    d_full += dfc
    return d_full, g_full


# ------------------------------------------------------------- orbital constraints
def _reference_orbital_data(state):
    """``(label, F, D, Gamma)`` of the orbital constraint, REFERENCE basis.

    The constraint lives where the reference wavefunction is defined, not where
    the perturbation is evaluated.  That distinction is load-bearing for
    ``reference=casci``: the RHF *occupied space* is the first ``nocc``
    REFERENCE orbitals, and the CASPT2 active block straddles that boundary, so
    a rotation inside the active block changes which orbitals are occupied and
    therefore changes ``F^RHF``.  Expressed in the semicanonical basis the
    condition would be neither diagonal nor invariant, and the active-active
    rotation could not be treated as redundant.

    ``reference=casci`` constrains the RHF Fock built from the RHF density;
    ``reference=casscf`` the (state-averaged) CASSCF generalized Fock built from
    the CASSCF RDMs, both in the reference orbitals.
    """
    setup = state.setup
    nbf = state.nbf
    if setup.options.reference == "casci":
        nocc = int(state.mol.data["nelec_A"])
        dens = np.zeros((nbf, nbf))
        for i in range(nocc):
            dens[i, i] = 2.0
        fock = state.h1e_ref + _jk_operator(dens, state.eri_ref)
        return "rhf", 0.5 * (fock + fock.T), dens, None

    gammas = np.zeros((setup.nact, setup.nact))
    bigg = np.zeros((setup.nact,) * 4)
    for k, r in enumerate(setup.roots):
        w = float(setup.weights[k])
        gammas += w * make_rdm1_spatial(state.ci_coeffs[:, r], state.act_dets,
                                        setup.nact)
        bigg += w * make_rdm2_spatial(state.ci_coeffs[:, r], state.act_dets,
                                      setup.nact)
    dens_sc, gam_sc = _full_rdms(gammas, bigg, setup.ncore, setup.nact, nbf)
    dens = _rotate_1(dens_sc, state.u_sc)
    gam = _rotate_2(gam_sc, state.u_sc)
    fock = _generalized_fock(dens, gam, state.h1e_ref, state.eri_ref)
    return "casscf", fock, dens, gam


def _orbital_pair_sets(state):
    """Orbital rotations that carry a constraint, in the reference basis.

    ``reference=casscf`` -- the CASSCF solution fixes only the non-redundant
    (inter-block) rotations; rotations inside the inactive, active and virtual
    blocks are redundant for the CASSCF energy AND for the PT2 energy (the block
    SPANS are all the PT2 depends on -- including the frozen core, whose orbitals
    are the lowest eigenvectors of the closed+active Fock restricted to the
    inactive span, which a within-inactive rotation leaves unchanged).  They also
    leave ``g_orb`` on the inter-block pairs at zero, so they are a true gauge
    freedom and are excluded from both sets.

    ``reference=casci`` -- the ACTIVE-ACTIVE rotations are additionally
    constrained, and they are not optional.  The RHF occupied space is the first
    ``nocc`` reference orbitals, and the CASPT2 active block straddles that
    boundary, so an active-active rotation moves an orbital across it: ``D^RHF``
    changes, ``F^RHF`` does not transform covariantly, and the constraint on the
    inter-block pairs is broken too.  RHF canonicality inside the active block is
    what fixes them.

    Within-inactive and within-virtual rotations stay out for ``casci`` as well.
    Those blocks lie wholly inside the RHF occupied and virtual spaces
    respectively, so such a rotation leaves ``D^RHF`` untouched and every
    constraint transforms covariantly -- the corresponding Jacobian columns are
    identically zero (which would make the system singular for a degenerate pair,
    e.g. the 2p shell of Li), and the direction is pure gauge.
    """
    setup = state.setup
    nbf = state.nbf
    pairs = _nonredundant_pairs(setup.ncore, setup.nact, nbf)
    if setup.options.reference == "casci":
        act = range(setup.ncore, setup.ncore + setup.nact)
        pairs = pairs + [(t, u) for t in act for u in act if t > u]
    return pairs


def _orbital_jacobian(state, kind, fock, dens, gam, pairs):
    """``dR_B/dt_A`` for every constraint ``B`` and orbital parameter ``A``.

    The rotation-derivative rule ``dM = [M, kappa]`` on each tensor index gives,
    for a quantity built from the integrals with a density that is FIXED in the
    MO basis,

        d/dt [ X(rotated integrals, fixed rho) ] = [X, kappa] - X(d rho)

    with ``d rho = [rho, kappa]`` the rotation the density would have undergone
    had it rotated too.  ``D^RHF`` is fixed in the MO basis precisely because the
    occupied space is defined by orbital INDEX, so it follows the orbitals; the
    same holds for the CASSCF RDMs at fixed CI vectors.  Applied to the RHF Fock
    this reproduces the textbook CPHF ``A`` matrix; applied to the CASSCF
    generalized Fock it is the CASSCF orbital Hessian.
    """
    nbf = state.nbf
    npar = len(pairs)
    if npar > MAX_ORBITAL_PAIRS:
        raise CASPT2GradientNotImplemented(
            "the analytic CASPT2 gradient assembles a dense %d x %d orbital "
            "response system, above the %d guard.  This module is "
            "validation-grade; reduce the basis." % (npar, npar, MAX_ORBITAL_PAIRS))

    jac = np.zeros((npar, npar))
    for a, (m, n) in enumerate(pairs):
        kap = _kappa_generator(nbf, m, n)
        if kind == "rhf":
            dfock = (_rotation_derivative_1(fock, kap)
                     - _jk_operator(_rotation_derivative_1(dens, kap),
                                    state.eri_ref))
            dref = dfock
        else:
            ddens = _rotation_derivative_1(dens, kap)
            dgam = _rotation_derivative_2(gam, kap)
            dfock = (_rotation_derivative_1(fock, kap)
                     - _generalized_fock(ddens, dgam, state.h1e_ref,
                                         state.eri_ref))
            dref = 2.0 * (dfock.T - dfock)
        for b, (p, q) in enumerate(pairs):
            jac[b, a] = dref[p, q]
    return jac


def _constraint_densities(state, kind, dens, gam, pairs, zvec):
    """Effective densities of ``sum_B z_B R_B``, in the reference basis."""
    nbf = state.nbf
    d1 = np.zeros((nbf, nbf))
    d2 = np.zeros((nbf,) * 4)

    if kind == "rhf":
        zmat = np.zeros((nbf, nbf))
        for b, (p, q) in enumerate(pairs):
            zmat[p, q] += 0.5 * zvec[b]
            zmat[q, p] += 0.5 * zvec[b]
        y1, y2 = _fock_like_densities(zmat, dens, state.eri_ref)
        d1 += y1
        d2 += y2
    else:
        zeta = np.zeros((nbf, nbf))
        for b, (p, q) in enumerate(pairs):
            zeta[q, p] += 2.0 * zvec[b]
            zeta[p, q] -= 2.0 * zvec[b]
        d1 += zeta.T @ dens
        d2 += 2.0 * np.einsum("ab,arst->brst", zeta, gam, optimize=True)
    return d1, d2


def _full_rdms_linear(dgamma, dbigg, dens, ncore, nact, nbf):
    """First variation of :func:`casscf._full_rdms` about the current RDMs.

    ``_full_rdms`` is ``D = 2 P_inact + gamma`` and ``Gamma = D (x) D - 1/2 ...``
    with the all-active block replaced by the CI two-particle RDM, so its
    linearization is the product rule on the separable part plus the direct
    replacement inside the active block.
    """
    act = list(range(ncore, ncore + nact))
    d1 = np.zeros((nbf, nbf))
    d1[np.ix_(act, act)] = dgamma
    g1 = (np.einsum("pq,rs->pqrs", d1, dens, optimize=True)
          + np.einsum("pq,rs->pqrs", dens, d1, optimize=True)
          - 0.5 * np.einsum("ps,rq->pqrs", d1, dens, optimize=True)
          - 0.5 * np.einsum("ps,rq->pqrs", dens, d1, optimize=True))
    g1[np.ix_(act, act, act, act)] = dbigg
    return d1, g1


def _ci_rdm_derivatives(state, k):
    """``(dgamma_I, dGamma_I)`` -- the CI derivative of root ``k``'s active RDMs.

    ``d gamma_tu / d c_I = <I|E_tu|c> + <c|E_tu|I>``, i.e. twice the symmetrized
    transition density between the unit CI vector and ``c``; the same holds for
    the two-particle RDM.  Evaluated with the existing same-vector engines
    through the polarization identity, in the small ACTIVE determinant space.
    """
    setup = state.setup
    nact = setup.nact
    c = state.ci_coeffs[:, setup.roots[k]]
    ndet = len(state.act_dets)
    dg = np.zeros((ndet, nact, nact))
    dG = np.zeros((ndet, nact, nact, nact, nact))
    unit = np.zeros(ndet)
    for i in range(ndet):
        unit[:] = 0.0
        unit[i] = 1.0
        t1, t2 = _transition_density(unit, c, state.act_dets, nact)
        dg[i] = 2.0 * t1
        dG[i] = 2.0 * t2
    return dg, dG


def _constraint_weights(state, kind, pairs):
    """``(a1_B, a2_B)`` such that ``R_B = sum a1 D + sum a2 Gamma``, SEMICANONICAL.

    Only the CI-dependent part matters here.  The RHF condition has none -- its
    density is the SCF one, independent of the CASCI vectors -- so
    ``reference=casci`` needs no coupling at all.  The CASSCF orbital gradient
    does, through the state-averaged RDMs; its weights are formed in the
    reference basis (where the constraint lives) and rotated into the
    semicanonical basis (where the CI vectors live).
    """
    nbf = state.nbf
    u = state.u_sc
    out = []
    for (p, q) in pairs:
        if kind == "rhf":
            out.append(None)
            continue
        zeta = np.zeros((nbf, nbf))
        zeta[q, p] = 2.0
        zeta[p, q] = -2.0
        a1 = zeta @ state.h1e_ref
        a2 = np.einsum("ab,brst->arst", zeta, state.eri_ref, optimize=True)
        # weights are dual to the densities, so they carry the INVERSE rotation
        out.append((_rotate_1(a1, u.T), _rotate_2(a2, u.T)))
    return out


def _constraint_ci_coupling(state, kind, dens_sc, pairs):
    """``dR_B/dc_r`` for every constraint that sees the CI vectors."""
    setup = state.setup
    nbf = state.nbf
    weights = _constraint_weights(state, kind, pairs)
    npar = len(pairs)
    ndet = len(state.act_dets)
    have = np.zeros(npar, dtype=bool)
    out = [[None] * state.nstate for _ in range(npar)]
    if not any(w is not None for w in weights):
        return out, have
    derivs = [_ci_rdm_derivatives(state, k) for k in range(state.nstate)]
    for b, w in enumerate(weights):
        if w is None:
            continue
        a1, a2 = w
        have[b] = True
        for k in range(state.nstate):
            dg, dG = derivs[k]
            wk = float(setup.weights[k])
            vec = np.zeros(ndet)
            for i in range(ndet):
                dd, gg = _full_rdms_linear(dg[i], dG[i], dens_sc, setup.ncore,
                                           setup.nact, nbf)
                val = float(np.einsum("pq,pq->", a1, dd, optimize=True))
                if a2 is not None:
                    val += float(np.einsum("pqrs,pqrs->", a2, gg, optimize=True))
                vec[i] = wk * val
            out[b][k] = vec
    return out, have


@dataclass
class RelaxedDensities:
    """MO-basis effective densities of one target state, plus diagnostics.

    ``stationarity`` is the quantity that must vanish: the residual of
    ``dL_total/dkappa = 0`` over the CONSTRAINED orbital pairs.  ``asymmetry``
    (``max |F - F^T|`` over all pairs) is informational and is *expected* to be
    nonzero on the redundant within-block pairs -- there the multiplier terms
    have a gradient the Lagrangian itself does not, and the corresponding
    rotation is pure gauge (the semicanonical frame is held fixed, ``dU/dx = 0``),
    so it contributes nothing to the nuclear derivative.
    """
    d1: np.ndarray
    d2: np.ndarray
    fock: np.ndarray
    xmat: np.ndarray
    asymmetry: float
    stationarity: float
    orbital_residual: float
    zorb: np.ndarray
    ci_multipliers: list
    energy: float


def _ci_response_densities(state, ymuls):
    """``sum_r TD(y_r, c_r)`` in the reduced determinant space."""
    norb = state.norb
    d1 = np.zeros((norb, norb))
    d2 = np.zeros((norb,) * 4)
    for k, y in enumerate(ymuls):
        if not np.any(y):
            continue
        a = _embed_ci(state, y)
        b = _embed_ci(state, state.ci_coeffs[:, state.setup.roots[k]])
        t1, t2 = _transition_density(a, b, state.full_dets, norb)
        d1 += t1
        d2 += t2
    return d1, d2


def _orbital_gradient_vector(state, d1_ref, d2_ref, pairs):
    """``dL/dt_A`` for every orbital parameter, from the effective densities.

    Everything is in the REFERENCE basis, which is where the orbital constraint
    and therefore the response live.
    """
    d1s = 0.5 * (d1_ref + d1_ref.T)
    d2s = _sym8(d2_ref)
    fock = _generalized_fock(d1s, d2s, state.h1e_ref, state.eri_ref)
    gmat = _orbital_gradient_matrix(fock)
    return np.array([gmat[p, q] for (p, q) in pairs]), fock


#: Condition number above which the orbital-response system is refused.  A
#: near-singular system means the reference conditions do not actually determine
#: the orbitals -- two orbitals degenerate under the constraint operator, say --
#: and the multipliers would be arbitrary.
MAX_RESPONSE_CONDITION = 1.0e10


def _solve_orbital_response(amat, rhs):
    """Solve the orbital multiplier equation, refusing a singular system.

    The redundant rotations are already excluded from the pair set, so a
    singular system here is a genuine degeneracy of the reference conditions
    (two orbitals the constraint cannot tell apart) rather than a known gauge
    direction, and its multipliers would be arbitrary.
    """
    if amat.size == 0:
        return np.zeros(0)
    cond = float(np.linalg.cond(amat))
    if not np.isfinite(cond) or cond > MAX_RESPONSE_CONDITION:
        raise RuntimeError(
            "the CASPT2 orbital-response system is singular (condition number "
            "%.3e).  Two reference orbitals are degenerate under the reference "
            "condition, so the orbital multipliers -- and the gradient -- are "
            "not determined.  This is usually a symmetry-degenerate pair that "
            "should be inside one orbital block, not split across two." % cond)
    return np.linalg.solve(amat, rhs)


def _check_orbital_constraints(state, kind, fock_c, pairs):
    """The constraint must actually hold at the point being differentiated."""
    if not pairs:
        return 0.0
    if kind == "rhf":
        res = max(abs(float(fock_c[p, q])) for (p, q) in pairs)
        if res > RHF_CANONICAL_TOL:
            raise RuntimeError(
                "the CASCI reference orbitals do not diagonalize the RHF Fock "
                "(largest off-diagonal element %.3e > %.1e).  The analytic "
                "gradient's orbital constraint is exactly that condition; "
                "tighten [scf] conv." % (res, RHF_CANONICAL_TOL))
    else:
        res = max(abs(2.0 * (float(fock_c[q, p]) - float(fock_c[p, q])))
                  for (p, q) in pairs)
        if res > CASSCF_STATIONARITY_TOL:
            raise RuntimeError(
                "the CASSCF reference is not stationary (largest "
                "orbital-rotation gradient %.3e > %.1e).  The analytic CASPT2 "
                "gradient's orbital constraint is g_orb = 0; tighten "
                "[casscf] gradient_norm_tol." % (res, CASSCF_STATIONARITY_TOL))
    return res


def relaxed_densities(state, target_index) -> RelaxedDensities:
    """Solve every response equation and return the fully relaxed densities.

    Returns MO-basis densities in the REFERENCE orbitals, so the caller
    transforms them to the AO basis with ``state.coeff_ref``.
    """
    nbf = state.nbf
    u = np.asarray(state.mixing[:, target_index], dtype=float)

    pieces = _unrelaxed_densities(state, u)
    yfock_xms, base_grad, _base = _xms_elimination(state, pieces)

    d1_red = pieces.d1.copy()
    d2_red = pieces.d2.copy()
    if np.any(yfock_xms):
        x1, x2 = _fock_like_densities(yfock_xms, state.dsa, state.eri)
        d1_red += x1
        d2_red += x2
    yfock_total = pieces.ymat + yfock_xms

    # The frozen/correlated inactive split responds through the inactive-block
    # Fock eigenproblem.  Its driving gradient is taken in the SEMICANONICAL
    # basis, before the CI and orbital multipliers: a within-inactive rotation
    # leaves the CI vectors untouched, so nothing here is circular.
    d1_pre, d2_pre = _unfold_frozen(d1_red, d2_red, state.nfrozen, nbf,
                                    state.eri_full)
    fock_pre = _generalized_fock(0.5 * (d1_pre + d1_pre.T), _sym8(d2_pre),
                                 state.h1e_full, state.eri_full)
    wfrozen = _frozen_split_weight(state, _orbital_gradient_matrix(fock_pre))

    ci_grad0 = _ci_gradients(state, pieces, yfock_total, base_grad,
                             yfock_full=wfrozen)
    ymul0 = _ci_multipliers(state, ci_grad0)

    pairs = _orbital_pair_sets(state)
    kind, fock_c, dens_c, gam_c = _reference_orbital_data(state)
    residual = _check_orbital_constraints(state, kind, fock_c, pairs)

    def total_ref(ymuls, with_constant=True):
        """Reduced-space densities plus the CI response, unfolded and rotated."""
        c1, c2 = _ci_response_densities(state, ymuls)
        base1 = d1_red + c1 if with_constant else c1
        base2 = d2_red + c2 if with_constant else c2
        f1, f2 = _unfold_frozen(base1, base2, state.nfrozen, nbf,
                                state.eri_full)
        if with_constant and np.any(wfrozen):
            w1, w2 = _fock_like_densities(wfrozen, state.dsa_full,
                                          state.eri_full)
            f1 = f1 + w1
            f2 = f2 + w2
        if state.nfrozen and not with_constant:
            # _unfold_frozen adds the constant frozen-shell density, which is
            # not part of a linear response column; subtract it.
            z1, z2 = _unfold_frozen(np.zeros_like(base1), np.zeros_like(base2),
                                    state.nfrozen, nbf, state.eri_full)
            f1 = f1 - z1
            f2 = f2 - z2
        return _rotate_1(f1, state.u_sc), _rotate_2(f2, state.u_sc)

    d1_ref0, d2_ref0 = total_ref(ymul0)
    b0, _f0 = _orbital_gradient_vector(state, d1_ref0, d2_ref0, pairs)

    npar = len(pairs)
    if npar:
        jac = _orbital_jacobian(state, kind, fock_c, dens_c, gam_c, pairs)
        coupling, have = _constraint_ci_coupling(state, kind, state.dsa_full,
                                                 pairs)
        mmat = np.zeros((npar, npar))
        for b in range(npar):
            if not have[b]:
                continue
            ymul_b = _ci_multipliers(state, coupling[b])
            f1, f2 = total_ref(ymul_b, with_constant=False)
            mmat[b, :], _fb = _orbital_gradient_vector(state, f1, f2, pairs)
        zorb = _solve_orbital_response((jac + mmat).T, -b0)
        if have.any():
            # _ci_multipliers is linear in its right-hand side, so solving once
            # with the assembled gradient is the same as accumulating the
            # per-column solutions -- and keeps one code path.
            ymul = _ci_multipliers(
                state, [ci_grad0[k] + sum(zorb[b] * coupling[b][k]
                                          for b in range(npar) if have[b])
                        for k in range(state.nstate)])
        else:
            ymul = ymul0
        d1_ref, d2_ref = total_ref(ymul)
        z1, z2 = _constraint_densities(state, kind, dens_c, gam_c, pairs, zorb)
        d1_ref = d1_ref + z1
        d2_ref = d2_ref + z2
    else:
        zorb = np.zeros(0)
        ymul = ymul0
        d1_ref, d2_ref = d1_ref0, d2_ref0

    d1s = 0.5 * (d1_ref + d1_ref.T)
    d2s = _sym8(d2_ref)
    fock = _generalized_fock(d1s, d2s, state.h1e_ref, state.eri_ref)
    asym = float(np.max(np.abs(fock - fock.T))) if nbf else 0.0
    gmat = _orbital_gradient_matrix(fock)
    stat = max((abs(float(gmat[p, q])) for (p, q) in pairs), default=0.0)
    xmat = 0.5 * (fock + fock.T)
    return RelaxedDensities(d1=d1s, d2=d2s, fock=fock, xmat=xmat,
                            asymmetry=asym, stationarity=stat,
                            orbital_residual=residual,
                            zorb=zorb, ci_multipliers=ymul,
                            energy=float(state.total_energies[target_index]))


# --------------------------------------------------------------------- AO layer
#: Relative threshold below which a factorization eigenvalue is dropped.
FACTOR_DROP_REL = 1.0e-14


def _pack_upper(mat):
    """Packed upper triangle in the layout the gradient kernels read.

    ``mathlib::pack_matrix`` uses LAPACK's column-major packed-upper format,
    which for a symmetric matrix is the same sequence
    :func:`oqp.library.fci._unpack_lower_triangle` reads back.
    """
    m = 0.5 * (np.asarray(mat, dtype=float) + np.asarray(mat, dtype=float).T)
    n = m.shape[0]
    return np.ascontiguousarray(
        np.concatenate([m[:c + 1, c] for c in range(n)]), dtype=np.float64)


def _factorize_two_particle(gamma_mo, coeff):
    """``Gamma^AO = sum_k lam_k A^k (x) A^k`` with ``A^k = C V^k C^T``.

    ``gamma_mo`` is already eight-fold symmetric, so as a matrix over the
    composite index ``(pq)`` it is symmetric and its eigendecomposition is an
    EXACT rewriting -- a change of storage, not a fit.  The antisymmetric half
    of the composite space is annihilated by the symmetry, so about half the
    eigenvalues are numerically zero and are dropped.
    """
    n = gamma_mo.shape[0]
    m2 = np.asarray(gamma_mo, dtype=float).reshape(n * n, n * n)
    m2 = 0.5 * (m2 + m2.T)
    evals, evecs = np.linalg.eigh(m2)
    big = float(np.max(np.abs(evals))) if evals.size else 0.0
    thr = FACTOR_DROP_REL * max(big, 1.0)
    keep = np.flatnonzero(np.abs(evals) > thr)
    lam = np.ascontiguousarray(evals[keep], dtype=np.float64)
    avec = np.empty((n, n, keep.size), dtype=np.float64)
    for idx, k in enumerate(keep):
        vmat = evecs[:, k].reshape(n, n)
        vmat = 0.5 * (vmat + vmat.T)
        avec[:, :, idx] = coeff @ vmat @ coeff.T
    return lam, np.ascontiguousarray(avec)


def _call_native_gradient(mol, nbf, dm_packed, wm_packed, lam, avec):
    """Contract the relaxed AO densities with the derivative integrals."""
    from oqp.library.fci import _lib_backend
    backend = _lib_backend()
    if backend is None:
        raise RuntimeError(
            "the analytic CASPT2 gradient needs liboqp; the native backend is "
            "unavailable.")
    lib, ffi = backend
    if not hasattr(lib, "caspt2_gradient"):
        raise RuntimeError(
            "liboqp does not provide the caspt2_gradient entry point; rebuild "
            "OpenQP against a source tree that contains "
            "source/modules/caspt2_gradient.F90.")
    info = np.zeros(2, dtype=np.float64)
    dm = np.ascontiguousarray(dm_packed, dtype=np.float64)
    wm = np.ascontiguousarray(wm_packed, dtype=np.float64)
    lm = np.ascontiguousarray(lam, dtype=np.float64)
    av = np.ascontiguousarray(avec, dtype=np.float64)
    status = int(lib.caspt2_gradient(
        mol.data._data, int(nbf),
        ffi.cast("double *", dm.ctypes.data),
        ffi.cast("double *", wm.ctypes.data),
        int(lm.size),
        ffi.cast("double *", av.ctypes.data),
        ffi.cast("double *", lm.ctypes.data),
        ffi.cast("double *", info.ctypes.data)))
    if status != 0:
        raise RuntimeError(
            "the CASPT2 gradient kernel failed with status %d "
            "(-21 non-Hartree-Fock Hamiltonian, -22 allocation, -25 size "
            "mismatch)." % status)
    return info


def _check_reported_energies(state):
    """Refuse to differentiate a state the energy path did not report."""
    try:
        published = np.asarray(state.mol.data["OQP::CASPT2_ENERGIES"], dtype=float)
    except Exception:                                    # noqa: BLE001
        published = None
    if published is None or published.size == 0:
        raise RuntimeError(
            "no CASPT2 energies are stored on the molecule; run the CASPT2 "
            "energy before the gradient.")
    mine = np.asarray(state.total_energies, dtype=float)
    if mine.size != published.size:
        raise RuntimeError(
            "the analytic gradient reconstructed %d state(s) but the energy "
            "path published %d; refusing to differentiate a different "
            "calculation." % (mine.size, published.size))
    worst = float(np.max(np.abs(mine - published)))
    if worst > ENERGY_MATCH_TOL:
        raise RuntimeError(
            "the analytic gradient's own reconstruction of the CASPT2 energy "
            "differs from the reported one by %.3e Eh (limit %.1e).  The "
            "gradient would not belong to the printed energy; this usually "
            "means the reported energy came from a different engine "
            "([pt2] engine=direct/fortran with a numerically different "
            "reduction) or from a non-reproducible reference."
            % (worst, ENERGY_MATCH_TOL))
    return worst


def _check_state_gap(state, index):
    if state.nstate < 2:
        return float("inf")
    energies = np.asarray(state.total_energies, dtype=float)
    others = np.delete(energies, index)
    gap = float(np.min(np.abs(others - energies[index])))
    if gap < STATE_GAP_TOL:
        raise RuntimeError(
            "the requested effective-Hamiltonian root is degenerate with a "
            "neighbour to %.3e Eh (limit %.1e): the mixing vector, and with it "
            "the state the gradient belongs to, is not well defined."
            % (gap, STATE_GAP_TOL))
    return gap


def caspt2_analytic_gradient(mol, grad_list=(0,), ref_energy=None):
    """Analytic CASPT2 / XMS-CASPT2 nuclear gradient.

    Returns an array of shape ``(nstate, natom, 3)`` in Hartree/Bohr, indexed by
    published PT2 state; rows that were not requested stay zero, matching the
    TDDFT convention the driver expects.
    """
    from oqp.library.caspt2_dyall import _log

    state = _build_state(mol, ref_energy)
    drift = _check_reported_energies(state)

    natom = int(mol.data["natom"])
    nstate = state.nstate
    out = np.zeros((max(nstate, 1), natom, 3))
    requested = [int(s) for s in np.atleast_1d(np.asarray(grad_list, dtype=int))]
    for s in requested:
        if s < 0 or s >= nstate:
            raise ValueError(
                "CASPT2 gradient state %d is out of range: this calculation "
                "published %d state(s), indices 0..%d."
                % (s, nstate, nstate - 1))

    nbf = state.nbf
    coeff = state.coeff_ref
    for s in requested:
        gap = _check_state_gap(state, s)
        rel = relaxed_densities(state, s)
        dao = coeff @ rel.d1 @ coeff.T
        xao = coeff @ rel.xmat @ coeff.T
        lam, avec = _factorize_two_particle(rel.d2, coeff)
        # The overlap-derivative kernel ADDS its argument and the gradient term
        # is -sum X^AO S^x, so the negated Lagrangian goes in -- the same sign
        # convention grd1::eijden uses for RHF.
        info = _call_native_gradient(mol, nbf, _pack_upper(dao),
                                     _pack_upper(-xao), lam, avec)
        out[s] = np.asarray(mol.get_grad(), dtype=float).reshape(natom, 3)
        _log(mol, "")
        _log(mol, "   ==============================================")
        _log(mol, "   PyOQP: analytic CASPT2 nuclear gradient")
        _log(mol, "   ==============================================")
        _log(mol, f"   {'differentiated PT2 state':<36}{s:>20d}")
        _log(mol, f"   {'energy of that state':<36}{rel.energy:>20.10f}")
        _log(mol, f"   {'energy reconstruction drift':<36}{drift:>20.3e}")
        _log(mol, f"   {'reference-constraint residual':<36}"
                  f"{rel.orbital_residual:>20.3e}")
        _log(mol, f"   {'Lagrangian stationarity residual':<36}"
                  f"{rel.stationarity:>20.3e}")
        _log(mol, f"   {'redundant-rotation asymmetry':<36}"
                  f"{rel.asymmetry:>20.3e}")
        _log(mol, f"   {'relaxed 2-PDM factorization rank':<36}"
                  f"{int(lam.size):>20d}")
        _log(mol, f"   {'relaxed AO density trace':<36}{float(info[1]):>20.6f}")
        _log(mol, f"   {'orbital-multiplier norm |z|':<36}"
                  f"{float(np.linalg.norm(rel.zorb)):>20.3e}")
        if np.isfinite(gap):
            _log(mol, f"   {'effective-Hamiltonian state gap':<36}{gap:>20.3e}")
        _log(mol, f"   {'semicanonical residual':<36}"
                  f"{state.h0_offdiag:>20.3e}")
        _log(mol, "   ==============================================")
    return out
