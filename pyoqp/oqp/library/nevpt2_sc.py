"""RDM-based strongly contracted NEVPT2 (SC-NEVPT2) for OpenQP.

This is the *internally contracted* flavour of NEVPT2 (Angeli-Cimiraglia-Malrieu
strong contraction, Dyall H0).  It complements the uncontracted Dyall-H0 path in
:mod:`oqp.library.caspt2_dyall` (``[pt2] h0=dyall``):

* The uncontracted path expands the full external determinant space and is the
  validated all-valence / OpenMolcas-matching engine; for CAS(4,4) it sits
  ~0.4 mEh below PySCF/ORCA SC-NEVPT2 because the decontracted first-order space
  spans strictly more than the eight strongly-contracted perturbers.
* This module reproduces the strong contraction exactly: the first-order space is
  partitioned into the eight Dyall subspaces Sr(-1), Si(+1), Sijrs(0), Sijr(+1),
  Srsi(-1), Srs(-2), Sij(+2), Sir(0), with **one contracted perturber per external
  orbital-index tuple**, and the energy of each is

      E_l = - sum_k  norm_k / (diff_k + h_k/norm_k),

  with ``norm_k = <Psi_l^k|Psi_l^k>``, ``h_k = <Psi_l^k|(H0-E0)|Psi_l^k>`` and
  ``diff_k`` the orbital-energy differences.  The Koopmans intermediates
  (a16/a17/a19/a22/a23/a25/a3/k27/a7/a9/a12/a13 and the eri-folded 4-pdm terms
  f3ca/f3ac) are the RDM contractions up to the active 3-RDM plus the active
  4-RDM, built here directly from the determinant CI vector (small active spaces).

The formulae follow PySCF ``pyscf/mrpt/nevpt2.py`` (Guo, Sun) and reproduce its
SC-NEVPT2 to < 1 uEh on CAS(2,2) and CAS(4,4) references.

Scope: validation-grade.  The active 4-RDM is built explicitly, so this is
limited to small active spaces (``[pt2] max_active`` guards ``nact``); the
strong contraction makes the *external* (core/virtual) space cheap.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

from oqp.library.fci import _determinants
from oqp.library.rdm import _bit_count

NUMERICAL_ZERO = 1.0e-14


def _ein(subscripts, *operands):
    """``np.einsum`` with contraction-order optimization forced on.

    NumPy's default einsum evaluates the whole index expression in one naive
    nested loop.  The Koopmans intermediates below are 8- and 9-index
    expressions, so that loop is where they live: measured on ``norb = 8``,

        f3ca ``kbij,rpqjkiac->pqrabc``  0.268 s (1.0 GFlop/s) -> 0.034 s (7.8x)
        a16  ``kbia,rpqcki->pqrabc``    0.034 s (1.0 GFlop/s) -> 0.0007 s (49x)
        Sr   ``ipqr,pqrabc,iabc->i``    0.010 s (0.6 GFlop/s) -> 0.0004 s (25x)

    ``optimize=True`` instead picks a pairwise contraction order and hands each
    step to BLAS.  It is the same sum evaluated in a different (pairwise, better
    conditioned) order -- agreement measured at 4e-13 absolute on O(1) random
    operands, and the pinned example energies are unchanged to all 10 printed
    digits.
    """
    return np.einsum(subscripts, *operands, optimize=True)


# --------------------------------------------------------------------------- RDMs
def _Epq_matrix(p, q, det_list, det_index, norb):
    """Build E_pq = sum_sigma a_{p,sigma}^+ a_{q,sigma} as a sparse matrix.

    E_pq is a fixed linear operator on the determinant basis, so its nonzero
    structure depends only on ``(p, q, det_list, norb)`` -- never on the vector
    it is applied to.  Building it once and reusing it turns the RDM assembly
    from a Python walk over every determinant per application into a sparse
    mat-vec; ``make_rdms`` applies these operators O(norb^6) times.

    Determinants are OpenQP-encoded: det = a | (b << norb) with occupied spatial
    orbitals as set bits (alpha in the low ``norb`` bits, beta in the high)."""
    rows = []
    cols = []
    data = []
    for k, det in enumerate(det_list):
        for shift in (0, norb):
            qbit = 1 << (q + shift)
            if not (det & qbit):
                continue
            mask_below_q = qbit - 1
            sgn = -1.0 if (_bit_count(det & mask_below_q) & 1) else 1.0
            det1 = det & ~qbit
            pbit = 1 << (p + shift)
            if det1 & pbit:
                continue
            mask_below_p = pbit - 1
            if _bit_count(det1 & mask_below_p) & 1:
                sgn = -sgn
            idx = det_index.get(det1 | pbit)
            if idx is not None:
                rows.append(idx)
                cols.append(k)
                data.append(sgn)

    ndet = len(det_list)
    if not data:
        return sparse.csr_matrix((ndet, ndet))
    return sparse.csr_matrix((data, (rows, cols)), shape=(ndet, ndet))


def _build_Epq_operators(det_list, det_index, norb):
    """All norb^2 single-excitation operators for one determinant basis."""
    return {
        (p, q): _Epq_matrix(p, q, det_list, det_index, norb)
        for p in range(norb)
        for q in range(norb)
    }



def _lib_make_rdms(ci, norb, det_list, upto):
    """dm1..dm_upto through the Fortran engine, or None when unavailable."""
    try:
        from oqp.library.fci import _lib_backend
    except Exception:
        return None
    backend = _lib_backend()
    if backend is None or norb > 31:
        return None
    lib, ffi = backend
    if not hasattr(lib, "nevpt2_make_rdms"):
        return None
    dets = np.ascontiguousarray(np.asarray(det_list, dtype=np.int64))
    civec = np.ascontiguousarray(np.asarray(ci, dtype=np.float64).reshape(-1))
    out = []
    for n in range(1, 5):
        out.append(np.zeros((norb,) * (2 * n), dtype=np.float64)
                   if n <= upto else np.zeros(1, dtype=np.float64))
    info = lib.nevpt2_make_rdms(
        int(norb), int(dets.size), ffi.cast("int64_t *", dets.ctypes.data),
        ffi.cast("double *", civec.ctypes.data), int(upto),
        *[ffi.cast("double *", a.ctypes.data) for a in out])
    if int(info) != 0:
        return None
    return tuple(out[:upto])


def make_rdms(ci, norb, active_nelec, upto=4):
    """Spin-free density matrices dm1..dm_upto in the PySCF make_dm123/1234
    convention (``dm_n[p,q,r,s,...] = <0| E_pq E_rs ... |0>``), built directly
    from the active-space determinant CI vector.

    Validated to machine precision against ``pyscf.fci.rdm.make_dm1234``."""
    det_list = _determinants(norb, active_nelec)
    det_index = {d: i for i, d in enumerate(det_list)}
    ci = np.asarray(ci, dtype=float).reshape(-1)

    # Fortran engine first; the NumPy assembly below is the same algorithm and
    # stays as the numerical pin / fallback.
    lib_result = _lib_make_rdms(ci, norb, det_list, upto)
    if lib_result is not None:
        return lib_result

    # Built once for this determinant basis and reused for every application
    # below; the dm4 loop alone applies them norb^6 times.
    Eop = _build_Epq_operators(det_list, det_index, norb)

    E1 = {}
    for p in range(norb):
        for q in range(norb):
            E1[(p, q)] = Eop[(p, q)] @ ci

    dm1 = np.zeros((norb,) * 2)
    for p in range(norb):
        for q in range(norb):
            dm1[p, q] = ci @ E1[(p, q)]
    if upto == 1:
        return (dm1,)

    dm2 = np.zeros((norb,) * 4)
    for p in range(norb):
        for q in range(norb):
            bra = E1[(q, p)]
            for r in range(norb):
                for s in range(norb):
                    dm2[p, q, r, s] = bra @ E1[(r, s)]
    if upto == 2:
        return dm1, dm2

    E2 = {}
    for t in range(norb):
        for u in range(norb):
            base = E1[(t, u)]
            for r in range(norb):
                for s in range(norb):
                    E2[(r, s, t, u)] = Eop[(r, s)] @ base
    dm3 = np.zeros((norb,) * 6)
    for p in range(norb):
        for q in range(norb):
            bra = E1[(q, p)]
            for r in range(norb):
                for s in range(norb):
                    for t in range(norb):
                        for u in range(norb):
                            dm3[p, q, r, s, t, u] = bra @ E2[(r, s, t, u)]
    if upto == 3:
        return dm1, dm2, dm3

    # dm4[p,q,r,s,t,u,v,w] = E1[(q,p)] . (E_rs E2[(t,u,v,w)]) is norb^8 inner
    # products over vectors of length ndet.  Written one at a time that is
    # ~1.7e7 separate dots; stacking the norb^2 bras into one matrix and the
    # norb^2 kets into another turns each (t,u,v,w) block into a single
    # (norb^2 x ndet) @ (ndet x norb^2) GEMM -- same flops, but through the
    # threaded BLAS instead of the Python interpreter.
    npair = norb * norb
    ndet = len(det_list)

    # bras[a, :] = E1[(q, p)] with a = p*norb + q, so the GEMM result is
    # already indexed [p*norb+q, r*norb+s].
    bras = np.empty((npair, ndet))
    for p in range(norb):
        for q in range(norb):
            bras[p * norb + q, :] = E1[(q, p)]

    dm4 = np.empty((norb,) * 8)
    kets = np.empty((ndet, npair))
    for t in range(norb):
        for u in range(norb):
            for v in range(norb):
                for w in range(norb):
                    base2 = E2[(t, u, v, w)]
                    for r in range(norb):
                        for s in range(norb):
                            kets[:, r * norb + s] = Eop[(r, s)] @ base2
                    block = bras @ kets
                    dm4[:, :, :, :, t, u, v, w] = block.reshape(
                        norb, norb, norb, norb)
    return dm1, dm2, dm3, dm4


# ---------------------------------------------------- Koopmans intermediates
# Ported verbatim from pyscf/mrpt/nevpt2.py.  ``h2e`` is the PHYSICIST-ordered
# active two-electron tensor (= chemist (pq|rs) transposed (0,2,1,3)).
def _f3ca_f3ac(h2e, dm4):
    """The two eri-folded 4-pdm intermediates that PySCF computes via the
    compiled ``NEVPTkern_cedf_aedf`` / ``NEVPTkern_aedf_ecdf`` kernels,
    reconstructed here from the genuine spin-free 4-RDM (machine precision)."""
    dm4ca = dm4.transpose(0, 1, 4, 5, 6, 7, 2, 3)
    f3ca = _ein('kbij,rpqjkiac->pqrabc', h2e, dm4ca).transpose(
        np.argsort((1, 4, 0, 2, 5, 3)))
    f3ac = _ein('ijka,rpqbjcik->pqrabc', h2e, dm4).transpose(
        np.argsort((1, 2, 0, 4, 3, 5)))
    return f3ca, f3ac


def _a16(h1e, h2e, dm3, f3ca, f3ac):
    a16 = -_ein('ib,rpqiac->pqrabc', h1e, dm3)
    a16 += _ein('ia,rpqbic->pqrabc', h1e, dm3)
    a16 -= _ein('ci,rpqbai->pqrabc', h1e, dm3)
    a16 -= f3ca.transpose(1, 4, 0, 2, 5, 3)
    a16 -= _ein('kbia,rpqcki->pqrabc', h2e, dm3)
    a16 -= _ein('kbaj,rpqjkc->pqrabc', h2e, dm3)
    a16 += _ein('cbij,rpqjai->pqrabc', h2e, dm3)
    fdm2 = _ein('kbij,rpajki->prab', h2e, dm3)
    for i in range(h1e.shape[0]):
        a16[:, i, :, :, :, i] += fdm2
    a16 += f3ac.transpose(1, 2, 0, 4, 3, 5)
    a16 -= f3ca.transpose(1, 2, 0, 4, 3, 5)
    a16 += _ein('jbij,rpqiac->pqrabc', h2e, dm3)
    a16 -= _ein('cjka,rpqbjk->pqrabc', h2e, dm3)
    a16 += _ein('jcij,rpqbai->pqrabc', h2e, dm3)
    return a16


def _a22(h1e, h2e, dm2, dm3, f3ca, f3ac):
    a22 = -_ein('pb,kipjac->ijkabc', h1e, dm3)
    a22 -= _ein('pa,kibjpc->ijkabc', h1e, dm3)
    a22 += _ein('cp,kibjap->ijkabc', h1e, dm3)
    a22 += _ein('cqra,kibjqr->ijkabc', h2e, dm3)
    a22 -= _ein('qcpq,kibjap->ijkabc', h2e, dm3)
    a22 -= f3ac.transpose(1, 5, 0, 2, 4, 3)
    fdm2 = _ein('pqrb,kiqcpr->ikbc', h2e, dm3)
    for i in range(h1e.shape[0]):
        a22[:, i, :, i, :, :] -= fdm2
    a22 -= _ein('pqab,kiqjpc->ijkabc', h2e, dm3)
    a22 += _ein('pcrb,kiajpr->ijkabc', h2e, dm3)
    a22 += _ein('cqrb,kiqjar->ijkabc', h2e, dm3)
    a22 -= f3ac.transpose(1, 3, 0, 4, 2, 5)
    a22 += f3ca.transpose(1, 3, 0, 4, 2, 5)
    a22 += 2.0 * _ein('jb,kiac->ijkabc', h1e, dm2)
    a22 += 2.0 * _ein('pjrb,kiprac->ijkabc', h2e, dm3)
    fdm2 = _ein('pa,kipc->ikac', h1e, dm2)
    fdm2 -= _ein('cp,kiap->ikac', h1e, dm2)
    fdm2 -= _ein('cqra,kiqr->ikac', h2e, dm2)
    fdm2 += _ein('qcpq,kiap->ikac', h2e, dm2)
    fdm2 += _ein('pqra,kiqcpr->ikac', h2e, dm3)
    fdm2 -= _ein('rcpq,kiaqrp->ikac', h2e, dm3)
    for i in range(h1e.shape[0]):
        a22[:, i, :, :, i, :] += fdm2 * 2
    return a22


def _a17(h1e, h2e, dm2, dm3):
    h1e = h1e - _ein('mjjn->mn', h2e)
    return -_ein('pi,cabi->abcp', h1e, dm2) \
        - _ein('kpij,cabjki->abcp', h2e, dm3)


def _a19(h1e, h2e, dm1, dm2):
    h1e = h1e - _ein('mjjn->mn', h2e)
    return -_ein('pi,ai->ap', h1e, dm1) \
        - _ein('kpij,ajki->ap', h2e, dm2)


def _a23(h1e, h2e, dm1, dm2, dm3):
    return -_ein('ip,caib->abcp', h1e, dm2) \
        - _ein('pijk,cajbik->abcp', h2e, dm3) \
        + 2.0 * _ein('bp,ca->abcp', h1e, dm1) \
        + 2.0 * _ein('pibk,caik->abcp', h2e, dm2)


def _a25(h1e, h2e, dm1, dm2):
    return -_ein('pi,ai->ap', h1e, dm1) \
        - _ein('pijk,jaik->ap', h2e, dm2) \
        + 2.0 * _ein('ap->pa', h1e) \
        + 2.0 * _ein('piaj,ij->ap', h2e, dm1)


def _hdm1(dm1):
    return 2.0 * np.eye(dm1.shape[0]) - dm1.transpose(1, 0)


def _hdm2(dm1, dm2):
    delta = np.eye(dm2.shape[0])
    dm2 = _ein('ikjl->ijkl', dm2) - _ein('jk,il->ijkl', delta, dm1)
    return _ein('klij->ijkl', dm2) \
        + _ein('il,kj->ijkl', delta, dm1) \
        + _ein('jk,li->ijkl', delta, dm1) \
        - 2.0 * _ein('ik,lj->ijkl', delta, dm1) \
        - 2.0 * _ein('jl,ki->ijkl', delta, dm1) \
        - 2.0 * _ein('il,jk->ijkl', delta, delta) \
        + 4.0 * _ein('ik,jl->ijkl', delta, delta)


def _hdm3(dm1, dm2, dm3, hdm1, hdm2):
    delta = np.eye(dm3.shape[0])
    return - _ein('pb,qrac->pqrabc', delta, hdm2) \
        - _ein('br,pqac->pqrabc', delta, hdm2) \
        + _ein('bq,prac->pqrabc', delta, hdm2) * 2.0 \
        + _ein('ap,bqcr->pqrabc', delta, dm2) * 2.0 \
        - _ein('ap,cr,bq->pqrabc', delta, delta, dm1) * 4.0 \
        + _ein('cr,bqap->pqrabc', delta, dm2) * 2.0 \
        - _ein('bqapcr->pqrabc', dm3) \
        + _ein('ar,pc,bq->pqrabc', delta, delta, dm1) * 2.0 \
        - _ein('ar,bqcp->pqrabc', delta, dm2)


def _a3(h1e, h2e, dm1, dm2, hdm1):
    delta = np.eye(dm2.shape[0])
    return _ein('ia,ip->pa', h1e, hdm1) \
        + 2.0 * _ein('ijka,pj,ik->pa', h2e, delta, dm1) \
        - _ein('ijka,jpik->pa', h2e, dm2)


def _k27(h1e, h2e, dm1, dm2):
    return -_ein('ai,pi->pa', h1e, dm1) \
        - _ein('iajk,pkij->pa', h2e, dm2) \
        + _ein('iaji,pj->pa', h2e, dm1)


def _a7(h1e, h2e, dm1, dm2, dm3):
    delta = np.eye(dm2.shape[0])
    rm2 = _ein('iljk->ijkl', dm2) - _ein('ik,jl->ijkl', dm1, delta)
    rm3 = _ein('injmkl->ijklmn', dm3) \
        - _ein('jn,imkl->ijklmn', delta, dm2) \
        - _ein('km,ijln->ijklmn', delta, rm2) \
        - _ein('kn,ijml->ijklmn', delta, rm2)
    a7 = -_ein('bi,pqia->pqab', h1e, rm2) \
        - _ein('ai,pqbi->pqab', h1e, rm2) \
        - _ein('kbij,pqkija->pqab', h2e, rm3) \
        - _ein('kaij,pqkibj->pqab', h2e, rm3) \
        - _ein('baij,pqij->pqab', h2e, rm2)
    return rm2, a7


def _a9(h1e, h2e, hdm1, hdm2, hdm3):
    a9 = _ein('ib,pqai->pqab', h1e, hdm2)
    a9 += _ein('ijib,pqaj->pqab', h2e, hdm2) * 2.0
    a9 -= _ein('ijjb,pqai->pqab', h2e, hdm2)
    a9 -= _ein('ijkb,pkqaij->pqab', h2e, hdm3)
    a9 += _ein('ia,pqib->pqab', h1e, hdm2)
    a9 -= _ein('ijja,pqib->pqab', h2e, hdm2)
    a9 -= _ein('ijba,pqji->pqab', h2e, hdm2)
    a9 += _ein('ijia,pqjb->pqab', h2e, hdm2) * 2.0
    a9 -= _ein('ijka,pqkjbi->pqab', h2e, hdm3)
    return a9


def _a12(h1e, h2e, dm1, dm2, dm3):
    return _ein('ia,qpib->pqab', h1e, dm2) \
        - _ein('bi,qpai->pqab', h1e, dm2) \
        + _ein('ijka,qpjbik->pqab', h2e, dm3) \
        - _ein('kbij,qpajki->pqab', h2e, dm3) \
        - _ein('bjka,qpjk->pqab', h2e, dm2) \
        + _ein('jbij,qpai->pqab', h2e, dm2)


def _a13(h1e, h2e, dm1, dm2, dm3):
    delta = np.eye(dm3.shape[0])
    a13 = -_ein('ia,qbip->pqab', h1e, dm2)
    a13 += _ein('pa,qb->pqab', h1e, dm1) * 2.0
    a13 += _ein('bi,qiap->pqab', h1e, dm2)
    a13 -= _ein('pa,bi,qi->pqab', delta, h1e, dm1) * 2.0
    a13 -= _ein('ijka,qbjpik->pqab', h2e, dm3)
    a13 += _ein('kbij,qjapki->pqab', h2e, dm3)
    a13 += _ein('blma,qmlp->pqab', h2e, dm2)
    a13 += _ein('kpma,qbkm->pqab', h2e, dm2) * 2.0
    a13 -= _ein('bpma,qm->pqab', h2e, dm1) * 2.0
    a13 -= _ein('lbkl,qkap->pqab', h2e, dm2)
    a13 -= _ein('ap,mbkl,qlmk->pqab', delta, h2e, dm2) * 2.0
    a13 += _ein('ap,lbkl,qk->pqab', delta, h2e, dm1) * 2.0
    return a13


def _norm_to_energy(norm, h, diff):
    idx = abs(norm) > NUMERICAL_ZERO
    ener = -(norm[idx] / (diff[idx] + h[idx] / norm[idx])).sum()
    return norm.sum(), ener


# --------------------------------------------------------------------- subspaces
def _Sr(dm1, dm2, dm3, dm4, h1e, h2e, h1e_v, h2e_v, e_virt, f3=None):
    f3ca, f3ac = _f3ca_f3ac(h2e, dm4) if f3 is None else f3
    a16 = _a16(h1e, h2e, dm3, f3ca, f3ac)
    a17 = _a17(h1e, h2e, dm2, dm3)
    a19 = _a19(h1e, h2e, dm1, dm2)
    ener = _ein('ipqr,pqrabc,iabc->i', h2e_v, a16, h2e_v) \
        + _ein('ipqr,pqra,ia->i', h2e_v, a17, h1e_v) * 2.0 \
        + _ein('ip,pa,ia->i', h1e_v, a19, h1e_v)
    norm = _ein('ipqr,rpqbac,iabc->i', h2e_v, dm3, h2e_v) \
        + _ein('ipqr,rpqa,ia->i', h2e_v, dm2, h1e_v) * 2.0 \
        + _ein('ip,pa,ia->i', h1e_v, dm1, h1e_v)
    return _norm_to_energy(norm, ener, e_virt)


def _Si(dm1, dm2, dm3, dm4, h1e, h2e, h1e_v, h2e_v, e_core, f3=None):
    f3ca, f3ac = _f3ca_f3ac(h2e, dm4) if f3 is None else f3
    a22 = _a22(h1e, h2e, dm2, dm3, f3ca, f3ac)
    a23 = _a23(h1e, h2e, dm1, dm2, dm3)
    a25 = _a25(h1e, h2e, dm1, dm2)
    ncas = dm1.shape[0]
    delta = np.eye(ncas)
    dm3_h = _ein('abef,cd->abcdef', dm2, delta) * 2 - dm3.transpose(0, 1, 3, 2, 4, 5)
    dm2_h = _ein('ab,cd->abcd', dm1, delta) * 2 - dm2.transpose(0, 1, 3, 2)
    dm1_h = 2 * delta - dm1.transpose(1, 0)
    ener = _ein('qpir,pqrabc,baic->i', h2e_v, a22, h2e_v) \
        + _ein('qpir,pqra,ai->i', h2e_v, a23, h1e_v) * 2.0 \
        + _ein('pi,pa,ai->i', h1e_v, a25, h1e_v)
    norm = _ein('qpir,rpqbac,baic->i', h2e_v, dm3_h, h2e_v) \
        + _ein('qpir,rpqa,ai->i', h2e_v, dm2_h, h1e_v) * 2.0 \
        + _ein('pi,pa,ai->i', h1e_v, dm1_h, h1e_v)
    return _norm_to_energy(norm, ener, -e_core)


def _Sijrs(e_core, e_virt, g_cvcv):
    ncore = len(e_core)
    nvirt = len(e_virt)
    if ncore == 0 or nvirt == 0:
        return 0.0, 0.0
    eia = e_core[:, None] - e_virt[None, :]
    norm = 0.0
    e = 0.0
    for i in range(ncore):
        gi = g_cvcv[i].transpose(1, 0, 2)  # (j,a,b)
        djba = (eia.reshape(-1, 1) + eia[i].reshape(1, -1)).ravel()
        t2i = (gi.ravel() / djba).reshape(ncore, nvirt, nvirt)
        theta = gi * 2 - gi.transpose(0, 2, 1)
        norm += _ein('jab,jab', gi, theta)
        e += _ein('jab,jab', t2i, theta)
    return norm, e


def _Sijr(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0, 0.0
    hdm1 = _hdm1(dm1)
    a3 = _a3(h1e, h2e, dm1, dm2, hdm1)
    ncore = len(e_core)
    ci_diag = np.diag_indices(ncore)
    ci_triu = np.triu_indices(ncore)
    norm = 2.0 * _ein('rpji,raji,pa->rji', h2e_v, h2e_v, hdm1) \
        - 1.0 * _ein('rpji,raij,pa->rji', h2e_v, h2e_v, hdm1)
    norm += norm.transpose(0, 2, 1)
    norm[:, ci_diag[0], ci_diag[1]] *= 0.5
    h = 2.0 * _ein('rpji,raji,pa->rji', h2e_v, h2e_v, a3) \
        - 1.0 * _ein('rpji,raij,pa->rji', h2e_v, h2e_v, a3)
    h += h.transpose(0, 2, 1)
    h[:, ci_diag[0], ci_diag[1]] *= 0.5
    diff = e_virt[:, None, None] - e_core[None, :, None] - e_core[None, None, :]
    return _norm_to_energy(norm[:, ci_triu[0], ci_triu[1]],
                           h[:, ci_triu[0], ci_triu[1]],
                           diff[:, ci_triu[0], ci_triu[1]])


def _Srsi(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0, 0.0
    k27 = _k27(h1e, h2e, dm1, dm2)
    nvirt = len(e_virt)
    vi_diag = np.diag_indices(nvirt)
    vi_triu = np.triu_indices(nvirt)
    norm = 2.0 * _ein('rsip,rsia,pa->rsi', h2e_v, h2e_v, dm1) \
        - 1.0 * _ein('rsip,sria,pa->rsi', h2e_v, h2e_v, dm1)
    norm += norm.transpose(1, 0, 2)
    norm[vi_diag] *= 0.5
    h = 2.0 * _ein('rsip,rsia,pa->rsi', h2e_v, h2e_v, k27) \
        - 1.0 * _ein('rsip,sria,pa->rsi', h2e_v, h2e_v, k27)
    h += h.transpose(1, 0, 2)
    h[vi_diag] *= 0.5
    diff = e_virt[:, None, None] + e_virt[None, :, None] - e_core[None, None, :]
    return _norm_to_energy(norm[vi_triu], h[vi_triu], diff[vi_triu])


def _Srs(dm1, dm2, dm3, h1e, h2e, h2e_v, e_virt):
    if len(e_virt) == 0:
        return 0.0, 0.0
    rm2, a7 = _a7(h1e, h2e, dm1, dm2, dm3)
    norm = 0.5 * _ein('rsqp,rsba,pqba->rs', h2e_v, h2e_v, rm2)
    h = 0.5 * _ein('rsqp,rsba,pqab->rs', h2e_v, h2e_v, a7)
    diff = e_virt[:, None] + e_virt[None, :]
    return _norm_to_energy(norm, h, diff)


def _Sij(dm1, dm2, dm3, h1e, h2e, h2e_v, e_core):
    if len(e_core) == 0:
        return 0.0, 0.0
    hdm1 = _hdm1(dm1)
    hdm2 = _hdm2(dm1, dm2)
    hdm3 = _hdm3(dm1, dm2, dm3, hdm1, hdm2)
    a9 = _a9(h1e, h2e, hdm1, hdm2, hdm3)
    norm = 0.5 * _ein('qpij,baij,pqab->ij', h2e_v, h2e_v, hdm2)
    h = 0.5 * _ein('qpij,baij,pqab->ij', h2e_v, h2e_v, a9)
    diff = e_core[:, None] + e_core[None, :]
    return _norm_to_energy(norm, h, -diff)


def _Sir(dm1, dm2, dm3, h1e, h2e, h1e_v, h2e_v1, h2e_v2, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0, 0.0
    norm = _ein('rpiq,raib,qpab->ir', h2e_v1, h2e_v1, dm2) * 2.0 \
        - _ein('rpiq,rabi,qpab->ir', h2e_v1, h2e_v2, dm2) \
        - _ein('rpqi,raib,qpab->ir', h2e_v2, h2e_v1, dm2) \
        + _ein('raqi,rabi,qb->ir', h2e_v2, h2e_v2, dm1) * 2.0 \
        - _ein('rpqi,rabi,qbap->ir', h2e_v2, h2e_v2, dm2) \
        + _ein('rpqi,raai,qp->ir', h2e_v2, h2e_v2, dm1) \
        + _ein('rpiq,ri,qp->ir', h2e_v1, h1e_v, dm1) * 4.0 \
        - _ein('rpqi,ri,qp->ir', h2e_v2, h1e_v, dm1) * 2.0 \
        + _ein('ri,ri->ir', h1e_v, h1e_v) * 2.0
    a12 = _a12(h1e, h2e, dm1, dm2, dm3)
    a13 = _a13(h1e, h2e, dm1, dm2, dm3)
    h = _ein('rpiq,raib,pqab->ir', h2e_v1, h2e_v1, a12) * 2.0 \
        - _ein('rpiq,rabi,pqab->ir', h2e_v1, h2e_v2, a12) \
        - _ein('rpqi,raib,pqab->ir', h2e_v2, h2e_v1, a12) \
        + _ein('rpqi,rabi,pqab->ir', h2e_v2, h2e_v2, a13)
    diff = e_core[:, None] - e_virt[None, :]
    return _norm_to_energy(norm, h, -diff)


# ---------------------------------------------------------------- block builder
def _blocks(h1e_mo, eri_mo, ncore, nact, eps):
    """Build all SC-NEVPT2 integral blocks from the full-MO (semicanonical)
    bare one-electron matrix ``h1e_mo`` and chemist two-electron tensor
    ``eri_mo`` (= the objects ``_transform_integrals`` returns)."""
    nbf = h1e_mo.shape[0]
    nocc = ncore + nact
    C = slice(0, ncore)
    A = slice(ncore, nocc)
    V = slice(nocc, nbf)

    vhf = (2.0 * _ein('pqii->pq', eri_mo[:, :, C, C])
           - _ein('piiq->pq', eri_mo[:, C, C, :]))
    heff = h1e_mo + vhf

    def phys(block):
        return block.transpose(0, 2, 1, 3)

    h1e = heff[A, A]
    h2e = eri_mo[A, A, A, A].transpose(0, 2, 1, 3)

    blocks = {'h1e': h1e, 'h2e': h2e}
    v = phys(eri_mo[V, A, A, A])
    blocks['Sr'] = (v, heff[V, A] - _ein('mbbn->mn', v))
    v = phys(eri_mo[A, C, A, A])
    blocks['Si'] = (v, heff[A, C])
    blocks['Sijrs'] = eri_mo[C, V, C, V]
    blocks['Sijr'] = phys(eri_mo[V, C, A, C])
    blocks['Srsi'] = phys(eri_mo[V, C, V, A])
    blocks['Srs'] = phys(eri_mo[V, A, V, A])
    blocks['Sij'] = phys(eri_mo[A, C, A, C])
    v1 = phys(eri_mo[V, C, A, A])
    v2 = phys(eri_mo[V, A, A, C])
    blocks['Sir'] = (v1, v2, heff[V, C])
    blocks['e_core'] = np.asarray(eps[:ncore], dtype=float)
    blocks['e_virt'] = np.asarray(eps[nocc:], dtype=float)
    return blocks


def sc_nevpt2_energy(h1e_mo, eri_mo, eps, ncore, nact, active_nelec, ci_vector):
    """Strongly contracted NEVPT2 correlation energy.

    Parameters
    ----------
    h1e_mo, eri_mo : full-MO bare one-electron matrix and chemist (pq|rs) tensor
        in the **semicanonical** orbital basis (``_transform_integrals`` output).
    eps : semicanonical orbital energies (diag generalized Fock).
    ncore, nact : inactive (doubly occupied) and active orbital counts.
    active_nelec : (na, nb) active electrons.
    ci_vector : the reference-root active-space CI vector over OpenQP determinants.

    Returns
    -------
    (e2, components) : total SC-NEVPT2 correlation and a dict of the eight
        per-subspace energies.
    """
    dm1, dm2, dm3, dm4 = make_rdms(ci_vector, nact, active_nelec, upto=4)
    B = _blocks(h1e_mo, eri_mo, ncore, nact, eps)
    h1e, h2e = B['h1e'], B['h2e']
    ec, ev = B['e_core'], B['e_virt']

    comp = {}
    # Sr and Si consume the SAME pair of eri-folded 4-pdm intermediates -- they
    # depend only on (h2e, dm4).  Building them once instead of once per subspace
    # removes the single most expensive contraction in the module (measured 1.9 s
    # of a 11.0 s CAS(8,8) block before the einsum-order fix).
    f3 = _f3ca_f3ac(h2e, dm4)
    v, h1v = B['Sr']
    comp['Sr'] = _Sr(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ev, f3=f3)[1]
    v, h1v = B['Si']
    comp['Si'] = _Si(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ec, f3=f3)[1]
    comp['Sijrs'] = _Sijrs(ec, ev, B['Sijrs'])[1]
    comp['Sijr'] = _Sijr(dm1, dm2, h1e, h2e, B['Sijr'], ec, ev)[1]
    comp['Srsi'] = _Srsi(dm1, dm2, h1e, h2e, B['Srsi'], ec, ev)[1]
    comp['Srs'] = _Srs(dm1, dm2, dm3, h1e, h2e, B['Srs'], ev)[1]
    comp['Sij'] = _Sij(dm1, dm2, dm3, h1e, h2e, B['Sij'], ec)[1]
    v1, v2, h1v = B['Sir']
    comp['Sir'] = _Sir(dm1, dm2, dm3, h1e, h2e, h1v, v1, v2, ec, ev)[1]
    e2 = float(sum(comp.values()))
    return e2, comp
