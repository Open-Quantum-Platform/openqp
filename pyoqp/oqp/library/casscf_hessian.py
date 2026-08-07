"""Analytic MCSCF orbital-rotation Hessian for the native CASSCF.

Where this now RUNS
-------------------
The assembly below runs inside liboqp: ``casscf_anhess.F90`` performs the same
sequence -- SA RDMs, ``casscf_hess_bmat``, the ``0.5 (B + B^T)``
symmetrization, the inactive-Fock fold, the active-space Hamiltonian and its
dense spectrum, then the per-root ``wmat`` / ``amp`` / ``relax`` chain -- with
the excitation stack cached for the whole run.  ``casscf_energy`` therefore
crosses the boundary once for ``[casscf] hessian = analytic`` too, where this
function crossed it ~10 times per macroiteration plus one per RDM.

The Python here is unchanged and stays as the numerical pin (both arms are
compared in the same process by ``tests/test_casscf_energy.py``) and as the
fallback for every declined run.  In particular the two REFUSALS are still
raised from here: a root degeneracy with non-zero orbital coupling, and an
active space past ``_MAX_STACK_BYTES``.  The driver returns a distinguished
status for each rather than inventing a message or, worse, a number.

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

from oqp.library.fci import _as_f64c, _as_i64c

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


def _z_matrix_general(D, G, t, T):
    """Reference form of :func:`_z_matrix`, assuming no symmetry of ``t``/``T``.

    Kept as the numerical pin for the symmetry argument documented on
    :func:`_z_matrix`; not used on the production path."""
    z = (D @ t - t @ D).T
    z += 0.5 * (
        np.einsum("nqrs,mqrs->mn", G, T, optimize=True)
        + np.einsum("qnrs,qmrs->mn", G, T, optimize=True)
        + np.einsum("qrns,qrms->mn", G, T, optimize=True)
        + np.einsum("qrsn,qrsm->mn", G, T, optimize=True)
    )
    return z


def _z_matrix(D, G, t, T):
    """Directional-derivative intermediate: ``dE_fix along K = sum_mn K_mn Z_mn``.

    ``Z(h, g)`` contracted over a pair reproduces the production gradient
    ``2 (F_qp - F_pq)``; evaluated on derivative integrals it yields the
    fixed-CI Hessian columns.

    The four four-index contractions of :func:`_z_matrix_general` are all the
    same number here, so only one is evaluated.  Both operands carry the
    pairwise index symmetry that makes them coincide:

      * ``G[p,q,r,s] = G[q,p,s,r] = G[r,s,p,q]`` -- the spin-summed chemist
        2-RDM (note ``G`` is NOT symmetric under ``q<->p`` alone; the exchange
        term breaks that, and the argument does not need it);
      * ``T`` inherits the full eight-fold symmetry of the real MO ERIs,
        because the one-index derivative of :func:`_one_index_derivative_g`
        applies the same generator to all four slots.

    Then ``G[a,n,c,d] T[a,m,c,d] = G[n,a,d,c] T[m,a,d,c]`` relabels term 2 onto
    term 1, and the ``[r,s,p,q]`` pair swap does the same for terms 3 and 4.
    Verified over 3220 live calls of a CASSCF/cc-pVDZ CAS(6,6) run: the four
    agree to 5.3e-14 relative.  ``_z_matrix_general`` remains as that pin (see
    ``tests/test_casscf_zmatrix.py``); do not use this form on tensors lacking
    those symmetries."""
    z = (D @ t - t @ D).T
    z += 2.0 * np.einsum("nqrs,mqrs->mn", G, T, optimize=True)
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
def _lib_excitation_stack(nact, dets):
    """Excitation-matrix stack through the Fortran engine, or None."""
    from oqp.library.fci import _lib_backend

    backend = _lib_backend()
    if backend is None or 2 * int(nact) > 62:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_excitation_stack"):
        return None
    det_arr = _as_i64c(dets)
    ndet = int(det_arr.size)
    # np.empty, not np.zeros: casscf_excitation_stack clears the whole array
    # itself, in parallel.  This buffer is nact**2 * ndet**2 doubles and the
    # clear is the dominant cost of the kernel, so paying for it twice -- once
    # single-threaded in C here, once in Fortran -- was most of its runtime.
    # A nonzero `info` means nothing was written; the caller discards it.
    stack = np.empty((int(nact), int(nact), ndet, ndet), dtype=np.float64)
    info = lib.casscf_excitation_stack(
        int(nact), ndet,
        ffi.cast("int64_t *", det_arr.ctypes.data),
        ffi.cast("double *", stack.ctypes.data))
    if int(info) != 0:
        return None
    return stack


def _excitation_matrices_python(nact, nalpha, nbeta):
    """NumPy/Python reference builder, kept as the pin for the Fortran engine."""
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


@lru_cache(maxsize=8)
def _excitation_matrices(nact, nalpha, nbeta):
    """Spin-summed excitation matrices ``(E_tu)_{row,col} = <row| E_tu |col>``.

    Same determinant ordering and bit convention (alpha bits low, beta bits
    high) as ``fci._determinants`` / ``rdm.make_rdm*``.  Built by the Fortran
    engine (``casscf_exc_stack.F90``) with
    :func:`_excitation_matrices_python` as the fallback and the numerical pin.
    Cached per active space; the returned stack is read-only."""
    dets = tuple(_determinants(int(nact), (int(nalpha), int(nbeta))))
    stack = _lib_excitation_stack(int(nact), dets)
    if stack is None:
        return _excitation_matrices_python(nact, nalpha, nbeta)
    stack.setflags(write=False)
    return dets, stack


def _active_operator_matrix(f, g, stack):
    """Dense determinant-basis matrix of ``sum f_tu E_tu + 1/2 sum g (E E - d E)``.

    Reference form, kept as the numerical pin for :func:`_active_hamiltonian`.
    The middle contraction costs ``nact^2 ndet^3``, which is what makes this the
    most expensive single call in the module."""
    ham = np.tensordot(f, stack, axes=([0, 1], [0, 1]))
    inter = np.tensordot(g, stack, axes=([2, 3], [0, 1]))  # (t, u, x, b)
    ham += 0.5 * np.tensordot(stack, inter, axes=([0, 1, 3], [0, 1, 2]))
    ham -= 0.5 * np.tensordot(np.einsum("tuuw->tw", g), stack, axes=([0, 1], [0, 1]))
    return ham


def _active_hamiltonian(f, g, stack, dets, nact):
    """Determinant-basis matrix of the active operator, through liboqp.

    ``sum f_tu E_tu + 1/2 sum g_tuvw (E_tu E_vw - d_uv E_tw)`` in the
    determinant basis *is* the CASCI Hamiltonian of the folded integrals
    ``(f, g)``, so the engine that already builds it -- ``fci_dense_hamiltonian``
    via :func:`fci._build_dense_hamiltonian` -- produces exactly this matrix.
    Routing here replaces the ``nact^2 ndet^3`` tensordot of
    :func:`_active_operator_matrix` with a validated Fortran kernel instead of a
    second implementation of the same object.  Agreement is 2.8e-14 absolute on
    a CAS(6,6) (see ``tests/test_casscf_active_hamiltonian.py``).

    Falls back to the NumPy reference when the native path is unavailable."""
    try:
        from oqp.library.fci import _build_dense_hamiltonian, _spin_orbital_integrals
    except Exception:
        return _active_operator_matrix(f, g, stack)
    if 2 * int(nact) > 62:            # determinant encoding must fit int64
        return _active_operator_matrix(f, g, stack)
    try:
        hspin, gspin = _spin_orbital_integrals(np.asarray(f), np.asarray(g))
        det_list = list(dets)
        det_index = {d: i for i, d in enumerate(det_list)}
        ham = _build_dense_hamiltonian(det_list, det_index, hspin, gspin,
                                       2 * int(nact), 0.0)
    except Exception:
        return _active_operator_matrix(f, g, stack)
    if ham is None:
        return _active_operator_matrix(f, g, stack)
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


def _hess_bmat_backend():
    """liboqp (lib, ffi) for the Fortran fixed-CI Hessian engine, or None."""
    from oqp.library.fci import _lib_backend

    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_hess_bmat"):
        return None
    return lib, ffi


def _lib_hess_bmat(nbf, ncore, nact, pairs, D, G, h1e, eri):
    """Fixed-CI Hessian columns and folded derivative integrals, in Fortran.

    Reproduces the per-pair Python loop -- ``_z_matrix`` over the one-index
    derivative tensor plus ``_fold_active`` -- but assembles both from the eight
    sparse slabs of the derivative, so the per-pair cost is O(nbf^4) instead of
    O(nbf^5) and no nbf^4 temporary is formed.  Returns ``(B, f_der, g_der)``,
    or ``None`` when the engine is unavailable."""
    backend = _hess_bmat_backend()
    if backend is None:
        return None
    lib, ffi = backend
    npar = len(pairs)
    pr = np.ascontiguousarray(np.asarray(pairs, dtype=np.int32).reshape(npar, 2))
    dd = _as_f64c(D)
    gg = _as_f64c(G)
    hh = _as_f64c(h1e)
    vv = _as_f64c(eri)
    B = np.zeros((npar, npar), dtype=np.float64)
    f_der = np.zeros((npar, nact, nact), dtype=np.float64)
    g_der = np.zeros((npar,) + (nact,) * 4, dtype=np.float64)
    lib.casscf_hess_bmat(
        int(nbf), int(ncore), int(nact), int(npar),
        ffi.cast("int32_t *", pr.ctypes.data),
        ffi.cast("double *", dd.ctypes.data),
        ffi.cast("double *", gg.ctypes.data),
        ffi.cast("double *", hh.ctypes.data),
        ffi.cast("double *", vv.ctypes.data),
        ffi.cast("double *", B.ctypes.data),
        ffi.cast("double *", f_der.ctypes.data),
        ffi.cast("double *", g_der.ctypes.data))
    return B, f_der, g_der


def _hess_amp_backend():
    """liboqp (lib, ffi) for the Fortran CI-relaxation amplitude engine, or None."""
    from oqp.library.fci import _lib_backend

    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_hess_amp"):
        return None
    return lib, ffi


def _lib_hess_relax(npar, ndet, iavg, ovl, weights, eps, e_i, amp, hess,
                    degeneracy_tol):
    """CI-relaxation accumulation through the Fortran engine.

    Returns ``True`` when the engine ran (``hess`` updated in place), ``False``
    when unavailable so the caller runs the Python loop.  Raises on a genuine
    root degeneracy, matching the Python path's message."""
    backend = _hess_amp_backend()
    if backend is None:
        return False
    lib, ffi = backend
    if not hasattr(lib, "casscf_hess_relax"):
        return False
    ov = _as_f64c(ovl)
    wt = _as_f64c(weights)
    ep = _as_f64c(eps)
    am = _as_f64c(amp)
    info = int(lib.casscf_hess_relax(
        int(npar), int(ndet), int(ov.shape[0]), int(iavg),
        ffi.cast("double *", ov.ctypes.data),
        ffi.cast("double *", wt.ctypes.data),
        ffi.cast("double *", ep.ctypes.data),
        float(e_i), float(degeneracy_tol), float(_COUPLING_NOISE),
        ffi.cast("double *", am.ctypes.data),
        ffi.cast("double *", hess.ctypes.data)))
    if info != 0:
        j = info - 1
        denom = float(e_i) - float(eps[j])
        raise ValueError(
            "analytic CASSCF Hessian: root degeneracy with non-zero orbital "
            f"coupling (|E_I - E_J| = {abs(denom):.2e}); the objective is not "
            "smooth here and no orbital Hessian is defined")
    return True


def _lib_hess_wmat(stack, civec):
    """``W_tua = (E_tu|c>)_a`` through the Fortran engine, or None.

    Replaces ``np.tensordot(stack, c, axes=([3], [0]))``; the C-order stack is
    already the right matrix for a single DGEMV, so this is one call with no
    repacking."""
    backend = _hess_amp_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "casscf_hess_wmat"):
        return None
    nact = int(stack.shape[0])
    ndet = int(stack.shape[2])
    st = _as_f64c(stack)
    cv = _as_f64c(civec)
    wmat = np.zeros((nact, nact, ndet), dtype=np.float64)
    lib.casscf_hess_wmat(
        nact, ndet,
        ffi.cast("double *", st.ctypes.data),
        ffi.cast("double *", cv.ctypes.data),
        ffi.cast("double *", wmat.ctypes.data))
    return wmat


def _lib_hess_amp(stack, f_der, g_der, wmat, vecs):
    """Projected derivative-operator amplitudes through the Fortran engine.

    Reproduces ``[_apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs
    for k in range(npar)]``, but batches the pair index into the GEMMs instead
    of running one einsum per pair.  Returns the ``(npar, ndet)`` array, or
    ``None`` when the engine is unavailable."""
    backend = _hess_amp_backend()
    if backend is None:
        return None
    lib, ffi = backend
    nact = int(stack.shape[0])
    ndet = int(stack.shape[2])
    npar = int(f_der.shape[0])
    st = _as_f64c(stack)
    fd = _as_f64c(f_der)
    gd = _as_f64c(g_der)
    wm = _as_f64c(wmat)
    vc = _as_f64c(vecs)
    amp = np.zeros((npar, ndet), dtype=np.float64)
    lib.casscf_hess_amp(
        nact, ndet, npar,
        ffi.cast("double *", st.ctypes.data),
        ffi.cast("double *", fd.ctypes.data),
        ffi.cast("double *", gd.ctypes.data),
        ffi.cast("double *", wm.ctypes.data),
        ffi.cast("double *", vc.ctypes.data),
        ffi.cast("double *", amp.ctypes.data))
    return amp


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
    built = _lib_hess_bmat(nbf, ncore, nact, pairs, D, G, h1e, eri)
    if built is None:
        B = np.empty((npar, npar))
        f_der = np.empty((npar, nact, nact))
        g_der = np.empty((npar, nact, nact, nact, nact))
        for l, (p, q) in enumerate(pairs):
            t_l = _one_index_derivative_h(h1e, p, q)
            T_l = _one_index_derivative_g(eri, p, q)
            B[:, l] = _pair_differences(_z_matrix(D, G, t_l, T_l), pairs)
            f_der[l], g_der[l] = _fold_active(t_l, T_l, ncore, nact)
    else:
        B, f_der, g_der = built
    hess = 0.5 * (B + B.T)

    # --- part 2: CI relaxation over the complete active spectrum
    f0, g0 = _fold_active(h1e, eri, ncore, nact)
    hact = _active_hamiltonian(f0, g0, stack, dets, nact)
    eps, vecs = _symmetric_eigh(hact)
    ovl = (ci[:, roots].T @ vecs) ** 2  # (n_avg, ndet) averaged-root overlaps

    for i_avg, (w_i, r) in enumerate(zip(weights, roots)):
        c = np.array(ci[:, r], dtype=float, copy=True)
        nrm = float(np.linalg.norm(c))
        if nrm <= 0.0:
            raise ValueError("CI vector for an averaged root has zero norm")
        c /= nrm
        e_i = float(c @ (hact @ c))
        wmat = _lib_hess_wmat(stack, c)
        if wmat is None:
            wmat = np.tensordot(stack, c, axes=([3], [0]))  # (nact, nact, ndet)
        amp = _lib_hess_amp(stack, f_der, g_der, wmat, vecs)
        if amp is None:
            amp = np.empty((npar, ndet))
            for k in range(npar):
                amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

        if not _lib_hess_relax(npar, ndet, i_avg, ovl, weights, eps, e_i,
                               amp, hess, degeneracy_tol):
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
