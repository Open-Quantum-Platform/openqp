"""Internally contracted CASPT2 (IC-CASPT2) for OpenQP.

Reuses the strongly-contracted NEVPT2 8-subspace first-order partitioning
(Sr, Si, Sijrs, Sijr, Srsi, Srs, Sij, Sir) from :mod:`oqp.library.nevpt2_sc`
with Fock's zeroth-order Hamiltonian (CASPT2 H0 = diag(eps)).

In SC-NEVPT2 each subspace perturber l contributes

    E_l = - norm_k / (diff_k + h_k / norm_k)

where ``h_k = <Psi_l^k|(H0 - E0)|Psi_l^k>`` is the (nonzero) fluctuation
potential of Dyall's H0 inside the contracted perturber.  For Fock H0 the
one-body operator is already diagonal in the MO basis so

    h_k = 0   (exactly),

and the expression reduces to the diagonal-denominator

    E_l = - norm_k / diff_k

with ``diff_k`` the orbital-energy denominator for that subspace.
The *norm* expressions (norms of the contracted perturbers) are identical
to SC-NEVPT2 because they depend only on integrals and RDMs, not on H0.

No external determinant space is formed -- the perturbers are contracted
combinations of active RDMs and integrals, so the cost scales as n_act^6
(3-RDM) instead of ndet x n_core^2 x n_virt^2.

The module reuses all :mod:`oqp.library.nevpt2_sc` building blocks:
Koopmans intermediates (a16, a22, a7, a9, a12, a13, f3ca_f3ac, ...),
subspace functions (_Sr, _Si, ..., _Sij, _Sir), and the integral-block
builder (_blocks).  Only the energy evaluation changes: h_k = 0.

Limitations
-----------
- Single-state only (MS/XMS would need off-diagonal coupling, future work).
- Fock H0 only (Dyall H0 -> SC-NEVPT2, which already exists).
- Closed-shell RHF reference only (same as the determinant-based CASPT2).
"""
from __future__ import annotations

import numpy as np

from oqp.library.nevpt2_sc import (
    _blocks,
    _ein,
)

NUMERICAL_ZERO = 1.0e-14


def ic_caspt2_energy(h1e_mo, eri_mo, eps, ncore, nact, active_nelec, ci_vector,
                     max_memory=None):
    """Internally contracted CASPT2 correlation energy (Fock H0).

    Parameters
    ----------
    h1e_mo, eri_mo : full-MO bare one-electron matrix and chemist (pq|rs) tensor
        in the **semicanonical** orbital basis.
    eps : semicanonical orbital energies (diag generalized Fock).
    ncore, nact : inactive (doubly occupied) and active orbital counts.
    active_nelec : (na, nb) active electrons.
    ci_vector : the reference-root active-space CI vector over OpenQP determinants.

    Returns
    -------
    (e2, components) : total IC-CASPT2 correlation energy and a dict of eight
        per-subspace energies.
    """
    from oqp.library.nevpt2_sc import make_rdms
    dm1, dm2, dm3, dm4 = make_rdms(ci_vector, nact, active_nelec, upto=4)

    B = _blocks(h1e_mo, eri_mo, ncore, nact, eps)
    h1e, h2e = B['h1e'], B['h2e']
    ec, ev = B['e_core'], B['e_virt']

    comp = {}
    from oqp.library.nevpt2_sc import _f3ca_f3ac
    f3 = _f3ca_f3ac(h2e, dm4)   # shared by Sr and Si

    v, h1v = B['Sr']
    comp['Sr'] = _ic_Sr(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ev, f3=f3)
    v, h1v = B['Si']
    comp['Si'] = _ic_Si(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ec, f3=f3)
    comp['Sijrs'] = _ic_Sijrs(ec, ev, B['Sijrs'])
    comp['Sijr'] = _ic_Sijr(dm1, dm2, h1e, h2e, B['Sijr'], ec, ev)
    comp['Srsi'] = _ic_Srsi(dm1, dm2, h1e, h2e, B['Srsi'], ec, ev)
    comp['Srs'] = _ic_Srs(dm1, dm2, dm3, h1e, h2e, B['Srs'], ev)
    comp['Sij'] = _ic_Sij(dm1, dm2, dm3, h1e, h2e, B['Sij'], ec)
    v1, v2, h1v = B['Sir']
    comp['Sir'] = _ic_Sir(dm1, dm2, dm3, h1e, h2e, h1v, v1, v2, ec, ev)

    e2 = float(sum(comp.values()))
    return e2, comp


# ------------------------------------------------------------------------- IC variants
# Each subspace function mirrors the SC-NEVPT2 equivalent but sets h_k = 0,
# so E_l = -norm_l / diff_l (diagonal denominators).

def _ic_Sr(dm1, dm2, dm3, dm4, h1e, h2e, h1e_v, h2e_v, e_virt, f3=None):
    from oqp.library.nevpt2_sc import _f3ca_f3ac, _a16, _a17, _a19
    f3ca, f3ac = _f3ca_f3ac(h2e, dm4) if f3 is None else f3
    a16 = _a16(h1e, h2e, dm3, f3ca, f3ac)
    a17 = _a17(h1e, h2e, dm2, dm3)
    a19 = _a19(h1e, h2e, dm1, dm2)
    norm = _ein('ipqr,rpqbac,iabc->i', h2e_v, dm3, h2e_v) \
        + _ein('ipqr,rpqa,ia->i', h2e_v, dm2, h1e_v) * 2.0 \
        + _ein('ip,pa,ia->i', h1e_v, dm1, h1e_v)
    return _ic_norm_to_energy(norm, e_virt)


def _ic_Si(dm1, dm2, dm3, dm4, h1e, h2e, h1e_v, h2e_v, e_core, f3=None):
    from oqp.library.nevpt2_sc import _f3ca_f3ac, _a22, _a23, _a25
    f3ca, f3ac = _f3ca_f3ac(h2e, dm4) if f3 is None else f3
    a22 = _a22(h1e, h2e, dm2, dm3, f3ca, f3ac)
    a23 = _a23(h1e, h2e, dm1, dm2, dm3)
    a25 = _a25(h1e, h2e, dm1, dm2)
    ncas = dm1.shape[0]
    delta = np.eye(ncas)
    dm3_h = _ein('abef,cd->abcdef', dm2, delta) * 2 - dm3.transpose(0, 1, 3, 2, 4, 5)
    dm2_h = _ein('ab,cd->abcd', dm1, delta) * 2 - dm2.transpose(0, 1, 3, 2)
    dm1_h = 2 * delta - dm1.transpose(1, 0)
    norm = _ein('qpir,rpqbac,baic->i', h2e_v, dm3_h, h2e_v) \
        + _ein('qpir,rpqa,ai->i', h2e_v, dm2_h, h1e_v) * 2.0 \
        + _ein('pi,pa,ai->i', h1e_v, dm1_h, h1e_v)
    return _ic_norm_to_energy(norm, -e_core)


def _ic_Sijrs(e_core, e_virt, g_cvcv):
    ncore = len(e_core)
    nvirt = len(e_virt)
    if ncore == 0 or nvirt == 0:
        return 0.0
    eia = e_core[:, None] - e_virt[None, :]
    e = 0.0
    for i in range(ncore):
        gi = g_cvcv[i].transpose(1, 0, 2)
        djba = (eia.reshape(-1, 1) + eia[i].reshape(1, -1)).ravel()
        t2i = (gi.ravel() / djba).reshape(ncore, nvirt, nvirt)
        theta = gi * 2 - gi.transpose(0, 2, 1)
        e += _ein('jab,jab', t2i, theta)
    return e


def _ic_Sijr(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0
    from oqp.library.nevpt2_sc import _hdm1, _a3
    hdm1 = _hdm1(dm1)
    a3 = _a3(h1e, h2e, dm1, dm2, hdm1)
    ncore = len(e_core)
    ci_diag = np.diag_indices(ncore)
    ci_triu = np.triu_indices(ncore)
    norm = 2.0 * _ein('rpji,raji,pa->rji', h2e_v, h2e_v, hdm1) \
        - 1.0 * _ein('rpji,raij,pa->rji', h2e_v, h2e_v, hdm1)
    norm += norm.transpose(0, 2, 1)
    norm[:, ci_diag[0], ci_diag[1]] *= 0.5
    diff = e_virt[:, None, None] - e_core[None, :, None] - e_core[None, None, :]
    return _ic_norm_to_energy(norm[:, ci_triu[0], ci_triu[1]],
                              diff[:, ci_triu[0], ci_triu[1]])


def _ic_Srsi(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0
    from oqp.library.nevpt2_sc import _k27
    k27 = _k27(h1e, h2e, dm1, dm2)
    nvirt = len(e_virt)
    vi_diag = np.diag_indices(nvirt)
    vi_triu = np.triu_indices(nvirt)
    norm = 2.0 * _ein('rsip,rsia,pa->rsi', h2e_v, h2e_v, dm1) \
        - 1.0 * _ein('rsip,sria,pa->rsi', h2e_v, h2e_v, dm1)
    norm += norm.transpose(1, 0, 2)
    norm[vi_diag] *= 0.5
    diff = e_virt[:, None, None] + e_virt[None, :, None] - e_core[None, None, :]
    return _ic_norm_to_energy(norm[vi_triu], diff[vi_triu])


def _ic_Srs(dm1, dm2, dm3, h1e, h2e, h2e_v, e_virt):
    if len(e_virt) == 0:
        return 0.0
    from oqp.library.nevpt2_sc import _a7
    rm2, a7 = _a7(h1e, h2e, dm1, dm2, dm3)
    norm = 0.5 * _ein('rsqp,rsba,pqba->rs', h2e_v, h2e_v, rm2)
    diff = e_virt[:, None] + e_virt[None, :]
    return _ic_norm_to_energy(norm, diff)


def _ic_Sij(dm1, dm2, dm3, h1e, h2e, h2e_v, e_core):
    if len(e_core) == 0:
        return 0.0
    from oqp.library.nevpt2_sc import _hdm1, _hdm2, _hdm3, _a9
    hdm1 = _hdm1(dm1)
    hdm2 = _hdm2(dm1, dm2)
    hdm3 = _hdm3(dm1, dm2, dm3, hdm1, hdm2)
    a9 = _a9(h1e, h2e, hdm1, hdm2, hdm3)
    norm = 0.5 * _ein('qpij,baij,pqab->ij', h2e_v, h2e_v, hdm2)
    diff = e_core[:, None] + e_core[None, :]
    return _ic_norm_to_energy(norm, -diff)


def _ic_Sir(dm1, dm2, dm3, h1e, h2e, h1e_v, h2e_v1, h2e_v2, e_core, e_virt):
    if len(e_core) == 0 or len(e_virt) == 0:
        return 0.0
    from oqp.library.nevpt2_sc import _a12, _a13
    norm = _ein('rpiq,raib,qpab->ir', h2e_v1, h2e_v1, dm2) * 2.0 \
        - _ein('rpiq,rabi,qpab->ir', h2e_v1, h2e_v2, dm2) \
        - _ein('rpqi,raib,qpab->ir', h2e_v2, h2e_v1, dm2) \
        + _ein('raqi,rabi,qb->ir', h2e_v2, h2e_v2, dm1) * 2.0 \
        - _ein('rpqi,rabi,qbap->ir', h2e_v2, h2e_v2, dm2) \
        + _ein('rpqi,raai,qp->ir', h2e_v2, h2e_v2, dm1) \
        + _ein('rpiq,ri,qp->ir', h2e_v1, h1e_v, dm1) * 4.0 \
        - _ein('rpqi,ri,qp->ir', h2e_v2, h1e_v, dm1) * 2.0 \
        + _ein('ri,ri->ir', h1e_v, h1e_v) * 2.0
    diff = e_core[:, None] - e_virt[None, :]
    return _ic_norm_to_energy(norm, -diff)


def _ic_norm_to_energy(norm, diff):
    """E = -sum_k norm_k / diff_k  (h_k = 0 for Fock H0)."""
    if np.isscalar(norm):
        return 0.0 if abs(diff) < NUMERICAL_ZERO else -float(norm) / diff
    idx = np.abs(diff) > NUMERICAL_ZERO
    if not np.any(idx):
        return 0.0
    return -float(np.sum(norm[idx] / diff[idx]))
