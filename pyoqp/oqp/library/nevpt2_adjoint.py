"""Exact adjoints (reverse-mode derivatives) of the SC-NEVPT2 energy functional.

This module differentiates :mod:`oqp.library.nevpt2_sc` analytically.  Nothing
here is a finite difference: every routine below is the closed-form derivative
of the corresponding forward routine, written as another contraction of the
same tensors.  For a multilinear step

    O = sum_{contracted} A B ...          (``np.einsum(sub, A, B, ...)``)

the derivative of a scalar ``E`` with respect to ``A`` follows from the chain
rule with no approximation,

    dE/dA = sum  (dE/dO)  B ...           (the same sum, re-routed),

which is what :func:`_ein_bar` forms.  The two non-multilinear steps -- the
strongly contracted denominator of :func:`oqp.library.nevpt2_sc._norm_to_energy`
and the ``Sijrs`` amplitude division -- are differentiated explicitly.

What the adjoints are
---------------------
:func:`sc_nevpt2_energy_adjoints` returns, for the SC-NEVPT2 correlation energy
``E2(h1e_mo, eri_mo, eps, ci)`` evaluated in the semicanonical basis,

    hbar[p,q]     = dE2 / d h1e_mo[p,q]        (bare one-electron MO matrix)
    gbar[p,q,r,s] = dE2 / d eri_mo[p,q,r,s]    (chemist (pq|rs))
    epsbar[p]     = dE2 / d eps[p]             (semicanonical orbital energies)
    cibar[K]      = dE2 / d ci[K]              (active-space CI coefficients)

``hbar``/``gbar`` here are the *explicit* integral dependence only: the further
dependence of ``eps`` on the integrals runs through the generalized Fock and is
resolved by the caller (:mod:`oqp.library.nevpt2_gradient`), which owns the
orbital-response step and therefore also owns the Fock chain rule.

Index conventions are inherited unchanged from :mod:`oqp.library.nevpt2_sc`:
``h2e`` is the PHYSICIST-ordered active two-electron tensor, ``eri_mo`` is
chemist ``(pq|rs)``, and the density matrices follow the PySCF
``make_dm1234`` convention ``dm_n[p,q,r,s,...] = <0| E_pq E_rs ... |0>``.
"""
from __future__ import annotations

import numpy as np

from oqp.library.nevpt2_sc import (
    NUMERICAL_ZERO,
    _blocks,
    _ein,
    _f3ca_f3ac,
    _hdm1,
    _hdm2,
    _hdm3,
    make_rdms,
)


# --------------------------------------------------------------------- einsum
def _split_subscripts(subscripts):
    """``'ab,bc->ac'`` -> ``(['ab', 'bc'], 'ac')``; implicit output is rejected.

    Every call site in :mod:`oqp.library.nevpt2_sc` writes the output
    explicitly, so refusing the implicit form keeps the adjoint from silently
    guessing an index order.
    """
    if "->" not in subscripts:
        raise ValueError("nevpt2 adjoints require an explicit '->' in einsum "
                         f"subscripts, got {subscripts!r}")
    lhs, out = subscripts.split("->")
    return [s.strip() for s in lhs.split(",")], out.strip()


def _ein_bar(subscripts, out_bar, operands, which, shape):
    """Exact ``d(scalar)/d(operands[which])`` for ``np.einsum(subscripts, *operands)``.

    ``out_bar`` is ``d(scalar)/d(output)``.  The result is the same sum with the
    output and the differentiated operand exchanged, which is the chain rule for
    a multilinear expression -- no approximation is involved.

    Two structural cases the plain exchange does not cover are handled here:

    * an index that appears ONLY in the differentiated operand was summed away,
      so the adjoint must broadcast back over it (an explicit ones factor);
    * an index REPEATED inside the differentiated operand (``'mjjn->mn'``,
      ``'qcpq,...'``) is a diagonal extraction, whose adjoint is a scatter onto
      that diagonal rather than a contraction.
    """
    subs, out_sub = _split_subscripts(subscripts)
    target = subs[which]
    # Deduplicate the target's own indices: 'mjjn' -> 'mjn'.  The forward step
    # read a diagonal, so the adjoint writes one.
    uniq = []
    for ch in target:
        if ch not in uniq:
            uniq.append(ch)
    uniq_sub = "".join(uniq)

    others = [(subs[i], operands[i]) for i in range(len(operands)) if i != which]
    have = set(out_sub)
    for s, _a in others:
        have.update(s)

    args = [out_bar]
    lhs = [out_sub]
    for s, a in others:
        lhs.append(s)
        args.append(a)
    # Indices private to the differentiated operand were contracted away in the
    # forward pass; the adjoint is constant along them.
    missing = [ch for ch in uniq_sub if ch not in have]
    if missing:
        dims = []
        for ch in missing:
            dims.append(shape[target.index(ch)])
        lhs.append("".join(missing))
        args.append(np.ones(dims, dtype=float))

    contracted = _ein(",".join(lhs) + "->" + uniq_sub, *args)
    if len(uniq_sub) == len(target):
        return contracted

    # Scatter the deduplicated result back onto the repeated-index diagonal.
    # The map (unique indices) -> (full index tuple) is injective, so a plain
    # fancy-index assignment is exact and needs no np.add.at accumulation.
    grids = np.indices([shape[target.index(ch)] for ch in uniq_sub])
    index = tuple(grids[uniq_sub.index(ch)] for ch in target)
    full = np.zeros(shape, dtype=float)
    full[index] = contracted
    return full


class _Bar:
    """Accumulator for the adjoints of the SC-NEVPT2 inputs.

    Keeps one buffer per differentiated tensor so the eight Dyall subspaces can
    add their contributions independently, exactly as the forward code adds
    their energies.
    """

    def __init__(self, nact, nbf):
        self.h1e = np.zeros((nact, nact))            # core-dressed active 1e
        self.h2e = np.zeros((nact,) * 4)             # physicist active 2e
        self.dm1 = np.zeros((nact,) * 2)
        self.dm2 = np.zeros((nact,) * 4)
        self.dm3 = np.zeros((nact,) * 6)
        self.dm4 = np.zeros((nact,) * 8)
        self.f3ca = np.zeros((nact,) * 6)
        self.f3ac = np.zeros((nact,) * 6)
        self.blocks = {}                             # per-subspace integral bars
        self.e_core = np.zeros(nbf)                  # sized by the caller's slice
        self.e_virt = np.zeros(nbf)

    def add_block(self, key, value):
        cur = self.blocks.get(key)
        if cur is None:
            self.blocks[key] = value
        else:
            self.blocks[key] = cur + value


# ------------------------------------------------------- denominator adjoints
def _norm_to_energy_bar(norm, h, diff, ener_bar):
    """Adjoint of ``nevpt2_sc._norm_to_energy``'s energy return.

    Forward (over the retained entries ``|norm| > NUMERICAL_ZERO``):

        d     = diff + h / norm
        ener  = - sum_k norm_k / d_k

    so, differentiating the quotient exactly,

        d ener / d norm = -1/d - h / (norm d^2)
        d ener / d h    =  1/d^2
        d ener / d diff =  norm / d^2
    """
    norm = np.asarray(norm, dtype=float)
    h = np.asarray(h, dtype=float)
    diff = np.asarray(diff, dtype=float)
    keep = np.abs(norm) > NUMERICAL_ZERO
    nb = np.zeros_like(norm)
    hb = np.zeros_like(norm)
    db = np.zeros_like(norm)
    if not np.any(keep):
        return nb, hb, db
    n_k = norm[keep]
    h_k = h[keep]
    d_k = diff[keep] + h_k / n_k
    inv = 1.0 / d_k
    inv2 = inv * inv
    nb[keep] = ener_bar * (-inv - h_k * inv2 / n_k)
    hb[keep] = ener_bar * inv2
    db[keep] = ener_bar * n_k * inv2
    return nb, hb, db


# ------------------------------------------- Koopmans intermediates (adjoints)
def _f3ca_f3ac_bar(h2e, dm4, ca_bar, ac_bar, bar):
    """Adjoint of ``nevpt2_sc._f3ca_f3ac``.

    Forward:
        dm4ca = dm4.transpose(0,1,4,5,6,7,2,3)
        f3ca  = ein('kbij,rpqjkiac->pqrabc', h2e, dm4ca).transpose(argsort(1,4,0,2,5,3))
        f3ac  = ein('ijka,rpqbjcik->pqrabc', h2e, dm4  ).transpose(argsort(1,2,0,4,3,5))
    """
    n = h2e.shape[0]
    dm4ca = dm4.transpose(0, 1, 4, 5, 6, 7, 2, 3)

    # undo the trailing transposes: X.transpose(argsort(p)) has adjoint
    # Xbar = out_bar.transpose(p)
    ca_core = ca_bar.transpose((1, 4, 0, 2, 5, 3))
    ac_core = ac_bar.transpose((1, 2, 0, 4, 3, 5))

    ops = (h2e, dm4ca)
    bar.h2e += _ein_bar('kbij,rpqjkiac->pqrabc', ca_core, ops, 0, h2e.shape)
    dm4ca_bar = _ein_bar('kbij,rpqjkiac->pqrabc', ca_core, ops, 1, dm4ca.shape)
    # inverse of dm4.transpose(0,1,4,5,6,7,2,3)
    bar.dm4 += dm4ca_bar.transpose(np.argsort((0, 1, 4, 5, 6, 7, 2, 3)))

    ops = (h2e, dm4)
    bar.h2e += _ein_bar('ijka,rpqbjcik->pqrabc', ac_core, ops, 0, h2e.shape)
    bar.dm4 += _ein_bar('ijka,rpqbjcik->pqrabc', ac_core, ops, 1, dm4.shape)


def _a16_bar(h1e, h2e, dm3, f3ca, f3ac, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a16`` (Sr subspace)."""
    n = h1e.shape[0]
    ob = out_bar

    def add(sub, ops, sign):
        for k, op in enumerate(ops):
            contrib = sign * _ein_bar(sub, ob, ops, k, op.shape)
            if op is h1e:
                bar.h1e += contrib
            elif op is h2e:
                bar.h2e += contrib
            elif op is dm3:
                bar.dm3 += contrib
            else:  # pragma: no cover - defensive
                raise AssertionError("unexpected a16 operand")

    add('ib,rpqiac->pqrabc', (h1e, dm3), -1.0)
    add('ia,rpqbic->pqrabc', (h1e, dm3), +1.0)
    add('ci,rpqbai->pqrabc', (h1e, dm3), -1.0)
    bar.f3ca += -ob.transpose(np.argsort((1, 4, 0, 2, 5, 3)))
    add('kbia,rpqcki->pqrabc', (h2e, dm3), -1.0)
    add('kbaj,rpqjkc->pqrabc', (h2e, dm3), -1.0)
    add('cbij,rpqjai->pqrabc', (h2e, dm3), +1.0)

    # fdm2 = ein('kbij,rpajki->prab', h2e, dm3);  a16[:,i,:,:,:,i] += fdm2
    fdm2_bar = np.zeros((n, n, n, n))
    for i in range(n):
        fdm2_bar += ob[:, i, :, :, :, i]
    ops = (h2e, dm3)
    bar.h2e += _ein_bar('kbij,rpajki->prab', fdm2_bar, ops, 0, h2e.shape)
    bar.dm3 += _ein_bar('kbij,rpajki->prab', fdm2_bar, ops, 1, dm3.shape)

    bar.f3ac += ob.transpose(np.argsort((1, 2, 0, 4, 3, 5)))
    bar.f3ca += -ob.transpose(np.argsort((1, 2, 0, 4, 3, 5)))
    add('jbij,rpqiac->pqrabc', (h2e, dm3), +1.0)
    add('cjka,rpqbjk->pqrabc', (h2e, dm3), -1.0)
    add('jcij,rpqbai->pqrabc', (h2e, dm3), +1.0)


def _a22_bar(h1e, h2e, dm2, dm3, f3ca, f3ac, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a22`` (Si subspace)."""
    n = h1e.shape[0]
    ob = out_bar

    def add(sub, ops, sign):
        for k, op in enumerate(ops):
            contrib = sign * _ein_bar(sub, ob, ops, k, op.shape)
            if op is h1e:
                bar.h1e += contrib
            elif op is h2e:
                bar.h2e += contrib
            elif op is dm2:
                bar.dm2 += contrib
            elif op is dm3:
                bar.dm3 += contrib
            else:  # pragma: no cover - defensive
                raise AssertionError("unexpected a22 operand")

    add('pb,kipjac->ijkabc', (h1e, dm3), -1.0)
    add('pa,kibjpc->ijkabc', (h1e, dm3), -1.0)
    add('cp,kibjap->ijkabc', (h1e, dm3), +1.0)
    add('cqra,kibjqr->ijkabc', (h2e, dm3), +1.0)
    add('qcpq,kibjap->ijkabc', (h2e, dm3), -1.0)
    bar.f3ac += -ob.transpose(np.argsort((1, 5, 0, 2, 4, 3)))

    # fdm2 = ein('pqrb,kiqcpr->ikbc', h2e, dm3);  a22[:,i,:,i,:,:] -= fdm2
    fdm2_bar = np.zeros((n, n, n, n))
    for i in range(n):
        fdm2_bar -= ob[:, i, :, i, :, :]
    ops = (h2e, dm3)
    bar.h2e += _ein_bar('pqrb,kiqcpr->ikbc', fdm2_bar, ops, 0, h2e.shape)
    bar.dm3 += _ein_bar('pqrb,kiqcpr->ikbc', fdm2_bar, ops, 1, dm3.shape)

    add('pqab,kiqjpc->ijkabc', (h2e, dm3), -1.0)
    add('pcrb,kiajpr->ijkabc', (h2e, dm3), +1.0)
    add('cqrb,kiqjar->ijkabc', (h2e, dm3), +1.0)
    bar.f3ac += -ob.transpose(np.argsort((1, 3, 0, 4, 2, 5)))
    bar.f3ca += ob.transpose(np.argsort((1, 3, 0, 4, 2, 5)))
    add('jb,kiac->ijkabc', (h1e, dm2), +2.0)
    add('pjrb,kiprac->ijkabc', (h2e, dm3), +2.0)

    # the second folded block: a22[:,i,:,:,i,:] += 2 * fdm2(ikac)
    f2_bar = np.zeros((n, n, n, n))
    for i in range(n):
        f2_bar += 2.0 * ob[:, i, :, :, i, :]
    for sub, ops, sign in (
            ('pa,kipc->ikac', (h1e, dm2), +1.0),
            ('cp,kiap->ikac', (h1e, dm2), -1.0),
            ('cqra,kiqr->ikac', (h2e, dm2), -1.0),
            ('qcpq,kiap->ikac', (h2e, dm2), +1.0),
            ('pqra,kiqcpr->ikac', (h2e, dm3), +1.0),
            ('rcpq,kiaqrp->ikac', (h2e, dm3), -1.0)):
        for k, op in enumerate(ops):
            contrib = sign * _ein_bar(sub, f2_bar, ops, k, op.shape)
            if op is h1e:
                bar.h1e += contrib
            elif op is h2e:
                bar.h2e += contrib
            elif op is dm2:
                bar.dm2 += contrib
            else:
                bar.dm3 += contrib


def _simple_bar(sub, ops, sign, out_bar, bar, names):
    """Add ``sign * d(scalar)/d(op)`` for every operand of one einsum term."""
    for k, op in enumerate(ops):
        name = names[k]
        if name is None:                    # constant (Kronecker delta)
            continue
        contrib = sign * _ein_bar(sub, out_bar, ops, k, op.shape)
        setattr(bar, name, getattr(bar, name) + contrib)


def _a17_bar(h1e, h2e, dm2, dm3, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a17``.

    Forward dresses the one-electron matrix first, ``h = h1e - ein('mjjn->mn', h2e)``,
    so the dressed matrix's adjoint has to be pushed back through that step.
    """
    n = h1e.shape[0]
    hd = h1e - _ein('mjjn->mn', h2e)
    hd_bar = np.zeros_like(h1e)

    ops = (hd, dm2)
    hd_bar += -_ein_bar('pi,cabi->abcp', out_bar, ops, 0, hd.shape)
    bar.dm2 += -_ein_bar('pi,cabi->abcp', out_bar, ops, 1, dm2.shape)
    ops = (h2e, dm3)
    bar.h2e += -_ein_bar('kpij,cabjki->abcp', out_bar, ops, 0, h2e.shape)
    bar.dm3 += -_ein_bar('kpij,cabjki->abcp', out_bar, ops, 1, dm3.shape)

    bar.h1e += hd_bar
    bar.h2e += -_ein_bar('mjjn->mn', hd_bar, (h2e,), 0, h2e.shape)


def _a19_bar(h1e, h2e, dm1, dm2, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a19`` (same one-electron dressing as a17)."""
    hd = h1e - _ein('mjjn->mn', h2e)
    hd_bar = np.zeros_like(h1e)

    ops = (hd, dm1)
    hd_bar += -_ein_bar('pi,ai->ap', out_bar, ops, 0, hd.shape)
    bar.dm1 += -_ein_bar('pi,ai->ap', out_bar, ops, 1, dm1.shape)
    ops = (h2e, dm2)
    bar.h2e += -_ein_bar('kpij,ajki->ap', out_bar, ops, 0, h2e.shape)
    bar.dm2 += -_ein_bar('kpij,ajki->ap', out_bar, ops, 1, dm2.shape)

    bar.h1e += hd_bar
    bar.h2e += -_ein_bar('mjjn->mn', hd_bar, (h2e,), 0, h2e.shape)


def _a23_bar(h1e, h2e, dm1, dm2, dm3, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a23``."""
    _simple_bar('ip,caib->abcp', (h1e, dm2), -1.0, out_bar, bar, ('h1e', 'dm2'))
    _simple_bar('pijk,cajbik->abcp', (h2e, dm3), -1.0, out_bar, bar, ('h2e', 'dm3'))
    _simple_bar('bp,ca->abcp', (h1e, dm1), +2.0, out_bar, bar, ('h1e', 'dm1'))
    _simple_bar('pibk,caik->abcp', (h2e, dm2), +2.0, out_bar, bar, ('h2e', 'dm2'))


def _a25_bar(h1e, h2e, dm1, dm2, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a25``."""
    _simple_bar('pi,ai->ap', (h1e, dm1), -1.0, out_bar, bar, ('h1e', 'dm1'))
    _simple_bar('pijk,jaik->ap', (h2e, dm2), -1.0, out_bar, bar, ('h2e', 'dm2'))
    bar.h1e += 2.0 * _ein_bar('ap->pa', out_bar, (h1e,), 0, h1e.shape)
    _simple_bar('piaj,ij->ap', (h2e, dm1), +2.0, out_bar, bar, ('h2e', 'dm1'))


def _a3_bar(h1e, h2e, dm1, dm2, hdm1, out_bar, bar, hdm1_bar):
    """Adjoint of ``nevpt2_sc._a3`` (Sijr subspace)."""
    n = dm2.shape[0]
    delta = np.eye(n)
    ops = (h1e, hdm1)
    bar.h1e += _ein_bar('ia,ip->pa', out_bar, ops, 0, h1e.shape)
    hdm1_bar += _ein_bar('ia,ip->pa', out_bar, ops, 1, hdm1.shape)
    ops = (h2e, delta, dm1)
    bar.h2e += 2.0 * _ein_bar('ijka,pj,ik->pa', out_bar, ops, 0, h2e.shape)
    bar.dm1 += 2.0 * _ein_bar('ijka,pj,ik->pa', out_bar, ops, 2, dm1.shape)
    _simple_bar('ijka,jpik->pa', (h2e, dm2), -1.0, out_bar, bar, ('h2e', 'dm2'))


def _k27_bar(h1e, h2e, dm1, dm2, out_bar, bar):
    """Adjoint of ``nevpt2_sc._k27`` (Srsi subspace)."""
    _simple_bar('ai,pi->pa', (h1e, dm1), -1.0, out_bar, bar, ('h1e', 'dm1'))
    _simple_bar('iajk,pkij->pa', (h2e, dm2), -1.0, out_bar, bar, ('h2e', 'dm2'))
    _simple_bar('iaji,pj->pa', (h2e, dm1), +1.0, out_bar, bar, ('h2e', 'dm1'))


def _a7_bar(h1e, h2e, dm1, dm2, dm3, rm2, rm2_bar, a7_bar, bar):
    """Adjoint of ``nevpt2_sc._a7`` (Srs subspace), including its rm2/rm3 build."""
    n = dm2.shape[0]
    delta = np.eye(n)
    rm3 = (_ein('injmkl->ijklmn', dm3)
           - _ein('jn,imkl->ijklmn', delta, dm2)
           - _ein('km,ijln->ijklmn', delta, rm2)
           - _ein('kn,ijml->ijklmn', delta, rm2))

    rm2_b = np.array(rm2_bar, dtype=float, copy=True)
    rm3_b = np.zeros_like(rm3)

    ob = a7_bar
    ops = (h1e, rm2)
    bar.h1e += -_ein_bar('bi,pqia->pqab', ob, ops, 0, h1e.shape)
    rm2_b += -_ein_bar('bi,pqia->pqab', ob, ops, 1, rm2.shape)
    ops = (h1e, rm2)
    bar.h1e += -_ein_bar('ai,pqbi->pqab', ob, ops, 0, h1e.shape)
    rm2_b += -_ein_bar('ai,pqbi->pqab', ob, ops, 1, rm2.shape)
    ops = (h2e, rm3)
    bar.h2e += -_ein_bar('kbij,pqkija->pqab', ob, ops, 0, h2e.shape)
    rm3_b += -_ein_bar('kbij,pqkija->pqab', ob, ops, 1, rm3.shape)
    ops = (h2e, rm3)
    bar.h2e += -_ein_bar('kaij,pqkibj->pqab', ob, ops, 0, h2e.shape)
    rm3_b += -_ein_bar('kaij,pqkibj->pqab', ob, ops, 1, rm3.shape)
    ops = (h2e, rm2)
    bar.h2e += -_ein_bar('baij,pqij->pqab', ob, ops, 0, h2e.shape)
    rm2_b += -_ein_bar('baij,pqij->pqab', ob, ops, 1, rm2.shape)

    # rm3 -> dm3, dm2, rm2
    bar.dm3 += _ein_bar('injmkl->ijklmn', rm3_b, (dm3,), 0, dm3.shape)
    ops = (delta, dm2)
    bar.dm2 += -_ein_bar('jn,imkl->ijklmn', rm3_b, ops, 1, dm2.shape)
    ops = (delta, rm2)
    rm2_b += -_ein_bar('km,ijln->ijklmn', rm3_b, ops, 1, rm2.shape)
    rm2_b += -_ein_bar('kn,ijml->ijklmn', rm3_b, ops, 1, rm2.shape)

    # rm2 = ein('iljk->ijkl', dm2) - ein('ik,jl->ijkl', dm1, delta)
    bar.dm2 += _ein_bar('iljk->ijkl', rm2_b, (dm2,), 0, dm2.shape)
    ops = (dm1, delta)
    bar.dm1 += -_ein_bar('ik,jl->ijkl', rm2_b, ops, 0, dm1.shape)


def _a9_bar(h1e, h2e, hdm2, hdm3, out_bar, hdm2_bar, hdm3_bar, bar):
    """Adjoint of ``nevpt2_sc._a9`` (Sij subspace)."""
    terms = (
        ('ib,pqai->pqab', (h1e, hdm2), +1.0, 'h1e', 2),
        ('ijib,pqaj->pqab', (h2e, hdm2), +2.0, 'h2e', 2),
        ('ijjb,pqai->pqab', (h2e, hdm2), -1.0, 'h2e', 2),
        ('ijkb,pkqaij->pqab', (h2e, hdm3), -1.0, 'h2e', 3),
        ('ia,pqib->pqab', (h1e, hdm2), +1.0, 'h1e', 2),
        ('ijja,pqib->pqab', (h2e, hdm2), -1.0, 'h2e', 2),
        ('ijba,pqji->pqab', (h2e, hdm2), -1.0, 'h2e', 2),
        ('ijia,pqjb->pqab', (h2e, hdm2), +2.0, 'h2e', 2),
        ('ijka,pqkjbi->pqab', (h2e, hdm3), -1.0, 'h2e', 3),
    )
    for sub, ops, sign, name, hkind in terms:
        contrib = sign * _ein_bar(sub, out_bar, ops, 0, ops[0].shape)
        setattr(bar, name, getattr(bar, name) + contrib)
        contrib = sign * _ein_bar(sub, out_bar, ops, 1, ops[1].shape)
        if hkind == 2:
            hdm2_bar += contrib
        else:
            hdm3_bar += contrib


def _a12_bar(h1e, h2e, dm2, dm3, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a12`` (Sir subspace)."""
    _simple_bar('ia,qpib->pqab', (h1e, dm2), +1.0, out_bar, bar, ('h1e', 'dm2'))
    _simple_bar('bi,qpai->pqab', (h1e, dm2), -1.0, out_bar, bar, ('h1e', 'dm2'))
    _simple_bar('ijka,qpjbik->pqab', (h2e, dm3), +1.0, out_bar, bar, ('h2e', 'dm3'))
    _simple_bar('kbij,qpajki->pqab', (h2e, dm3), -1.0, out_bar, bar, ('h2e', 'dm3'))
    _simple_bar('bjka,qpjk->pqab', (h2e, dm2), -1.0, out_bar, bar, ('h2e', 'dm2'))
    _simple_bar('jbij,qpai->pqab', (h2e, dm2), +1.0, out_bar, bar, ('h2e', 'dm2'))


def _a13_bar(h1e, h2e, dm1, dm2, dm3, out_bar, bar):
    """Adjoint of ``nevpt2_sc._a13`` (Sir subspace)."""
    n = dm3.shape[0]
    delta = np.eye(n)
    _simple_bar('ia,qbip->pqab', (h1e, dm2), -1.0, out_bar, bar, ('h1e', 'dm2'))
    _simple_bar('pa,qb->pqab', (h1e, dm1), +2.0, out_bar, bar, ('h1e', 'dm1'))
    _simple_bar('bi,qiap->pqab', (h1e, dm2), +1.0, out_bar, bar, ('h1e', 'dm2'))
    _simple_bar('pa,bi,qi->pqab', (delta, h1e, dm1), -2.0, out_bar, bar,
                (None, 'h1e', 'dm1'))
    _simple_bar('ijka,qbjpik->pqab', (h2e, dm3), -1.0, out_bar, bar, ('h2e', 'dm3'))
    _simple_bar('kbij,qjapki->pqab', (h2e, dm3), +1.0, out_bar, bar, ('h2e', 'dm3'))
    _simple_bar('blma,qmlp->pqab', (h2e, dm2), +1.0, out_bar, bar, ('h2e', 'dm2'))
    _simple_bar('kpma,qbkm->pqab', (h2e, dm2), +2.0, out_bar, bar, ('h2e', 'dm2'))
    _simple_bar('bpma,qm->pqab', (h2e, dm1), -2.0, out_bar, bar, ('h2e', 'dm1'))
    _simple_bar('lbkl,qkap->pqab', (h2e, dm2), -1.0, out_bar, bar, ('h2e', 'dm2'))
    _simple_bar('ap,mbkl,qlmk->pqab', (delta, h2e, dm2), -2.0, out_bar, bar,
                (None, 'h2e', 'dm2'))
    _simple_bar('ap,lbkl,qk->pqab', (delta, h2e, dm1), +2.0, out_bar, bar,
                (None, 'h2e', 'dm1'))


# ------------------------------------------------------------ hole-RDM chains
def _hdm1_bar(hdm1_bar, bar):
    """``hdm1 = 2 I - dm1^T``."""
    bar.dm1 += -hdm1_bar.T


def _hdm2_bar(dm1, dm2, hdm2_bar, bar):
    """Adjoint of ``nevpt2_sc._hdm2``.

    Forward:
        t   = ein('ikjl->ijkl', dm2) - ein('jk,il->ijkl', delta, dm1)
        out = ein('klij->ijkl', t) + ein('il,kj->ijkl', delta, dm1)
                                   + ein('jk,li->ijkl', delta, dm1)
                                   - 2 ein('ik,lj->ijkl', delta, dm1)
                                   - 2 ein('jl,ki->ijkl', delta, dm1)
                                   - 2 ein('il,jk->ijkl', delta, delta)
                                   + 4 ein('ik,jl->ijkl', delta, delta)
    """
    n = dm1.shape[0]
    delta = np.eye(n)
    ob = hdm2_bar
    t_bar = _ein_bar('klij->ijkl', ob, (np.zeros((n,) * 4),), 0, (n,) * 4)
    for sub, ops in (('il,kj->ijkl', (delta, dm1)),
                     ('jk,li->ijkl', (delta, dm1))):
        bar.dm1 += _ein_bar(sub, ob, ops, 1, dm1.shape)
    for sub, ops in (('ik,lj->ijkl', (delta, dm1)),
                     ('jl,ki->ijkl', (delta, dm1))):
        bar.dm1 += -2.0 * _ein_bar(sub, ob, ops, 1, dm1.shape)
    # t -> dm2, dm1
    bar.dm2 += _ein_bar('ikjl->ijkl', t_bar, (dm2,), 0, dm2.shape)
    bar.dm1 += -_ein_bar('jk,il->ijkl', t_bar, (delta, dm1), 1, dm1.shape)


def _hdm3_bar(dm1, dm2, dm3, hdm2, hdm3_bar, bar, hdm2_bar):
    """Adjoint of ``nevpt2_sc._hdm3``."""
    n = dm3.shape[0]
    delta = np.eye(n)
    ob = hdm3_bar
    terms = (
        ('pb,qrac->pqrabc', (delta, hdm2), -1.0, 'hdm2'),
        ('br,pqac->pqrabc', (delta, hdm2), -1.0, 'hdm2'),
        ('bq,prac->pqrabc', (delta, hdm2), +2.0, 'hdm2'),
        ('ap,bqcr->pqrabc', (delta, dm2), +2.0, 'dm2'),
        ('ap,cr,bq->pqrabc', (delta, delta, dm1), -4.0, 'dm1'),
        ('cr,bqap->pqrabc', (delta, dm2), +2.0, 'dm2'),
        ('bqapcr->pqrabc', (dm3,), -1.0, 'dm3'),
        ('ar,pc,bq->pqrabc', (delta, delta, dm1), +2.0, 'dm1'),
        ('ar,bqcp->pqrabc', (delta, dm2), -1.0, 'dm2'),
    )
    for sub, ops, sign, target in terms:
        which = len(ops) - 1          # the non-delta operand is always last
        contrib = sign * _ein_bar(sub, ob, ops, which, ops[which].shape)
        if target == 'hdm2':
            hdm2_bar += contrib
        else:
            setattr(bar, target, getattr(bar, target) + contrib)


# ---------------------------------------------------------------- subspaces
def _Sr_bar(dm1, dm2, dm3, h1e, h2e, h1e_v, h2e_v, e_virt, f3, bar):
    """Adjoint of ``nevpt2_sc._Sr``: returns the (h2e_v, h1e_v, e_virt) bars."""
    f3ca, f3ac = f3
    from oqp.library.nevpt2_sc import _a16, _a17, _a19
    a16 = _a16(h1e, h2e, dm3, f3ca, f3ac)
    a17 = _a17(h1e, h2e, dm2, dm3)
    a19 = _a19(h1e, h2e, dm1, dm2)

    ener = (_ein('ipqr,pqrabc,iabc->i', h2e_v, a16, h2e_v)
            + _ein('ipqr,pqra,ia->i', h2e_v, a17, h1e_v) * 2.0
            + _ein('ip,pa,ia->i', h1e_v, a19, h1e_v))
    norm = (_ein('ipqr,rpqbac,iabc->i', h2e_v, dm3, h2e_v)
            + _ein('ipqr,rpqa,ia->i', h2e_v, dm2, h1e_v) * 2.0
            + _ein('ip,pa,ia->i', h1e_v, dm1, h1e_v))
    nb, hb, db = _norm_to_energy_bar(norm, ener, e_virt, 1.0)

    v_bar = np.zeros_like(h2e_v)
    h1v_bar = np.zeros_like(h1e_v)
    a16_bar = np.zeros_like(a16)
    a17_bar = np.zeros_like(a17)
    a19_bar = np.zeros_like(a19)

    ops = (h2e_v, a16, h2e_v)
    v_bar += _ein_bar('ipqr,pqrabc,iabc->i', hb, ops, 0, h2e_v.shape)
    a16_bar += _ein_bar('ipqr,pqrabc,iabc->i', hb, ops, 1, a16.shape)
    v_bar += _ein_bar('ipqr,pqrabc,iabc->i', hb, ops, 2, h2e_v.shape)
    ops = (h2e_v, a17, h1e_v)
    v_bar += 2.0 * _ein_bar('ipqr,pqra,ia->i', hb, ops, 0, h2e_v.shape)
    a17_bar += 2.0 * _ein_bar('ipqr,pqra,ia->i', hb, ops, 1, a17.shape)
    h1v_bar += 2.0 * _ein_bar('ipqr,pqra,ia->i', hb, ops, 2, h1e_v.shape)
    ops = (h1e_v, a19, h1e_v)
    h1v_bar += _ein_bar('ip,pa,ia->i', hb, ops, 0, h1e_v.shape)
    a19_bar += _ein_bar('ip,pa,ia->i', hb, ops, 1, a19.shape)
    h1v_bar += _ein_bar('ip,pa,ia->i', hb, ops, 2, h1e_v.shape)

    ops = (h2e_v, dm3, h2e_v)
    v_bar += _ein_bar('ipqr,rpqbac,iabc->i', nb, ops, 0, h2e_v.shape)
    bar.dm3 += _ein_bar('ipqr,rpqbac,iabc->i', nb, ops, 1, dm3.shape)
    v_bar += _ein_bar('ipqr,rpqbac,iabc->i', nb, ops, 2, h2e_v.shape)
    ops = (h2e_v, dm2, h1e_v)
    v_bar += 2.0 * _ein_bar('ipqr,rpqa,ia->i', nb, ops, 0, h2e_v.shape)
    bar.dm2 += 2.0 * _ein_bar('ipqr,rpqa,ia->i', nb, ops, 1, dm2.shape)
    h1v_bar += 2.0 * _ein_bar('ipqr,rpqa,ia->i', nb, ops, 2, h1e_v.shape)
    ops = (h1e_v, dm1, h1e_v)
    h1v_bar += _ein_bar('ip,pa,ia->i', nb, ops, 0, h1e_v.shape)
    bar.dm1 += _ein_bar('ip,pa,ia->i', nb, ops, 1, dm1.shape)
    h1v_bar += _ein_bar('ip,pa,ia->i', nb, ops, 2, h1e_v.shape)

    _a16_bar(h1e, h2e, dm3, f3ca, f3ac, a16_bar, bar)
    _a17_bar(h1e, h2e, dm2, dm3, a17_bar, bar)
    _a19_bar(h1e, h2e, dm1, dm2, a19_bar, bar)
    return v_bar, h1v_bar, db


def _Si_bar(dm1, dm2, dm3, h1e, h2e, h1e_v, h2e_v, e_core, f3, bar):
    """Adjoint of ``nevpt2_sc._Si``: returns the (h2e_v, h1e_v, e_core) bars."""
    f3ca, f3ac = f3
    from oqp.library.nevpt2_sc import _a22, _a23, _a25
    a22 = _a22(h1e, h2e, dm2, dm3, f3ca, f3ac)
    a23 = _a23(h1e, h2e, dm1, dm2, dm3)
    a25 = _a25(h1e, h2e, dm1, dm2)
    ncas = dm1.shape[0]
    delta = np.eye(ncas)
    dm3_h = _ein('abef,cd->abcdef', dm2, delta) * 2 - dm3.transpose(0, 1, 3, 2, 4, 5)
    dm2_h = _ein('ab,cd->abcd', dm1, delta) * 2 - dm2.transpose(0, 1, 3, 2)
    dm1_h = 2 * delta - dm1.transpose(1, 0)

    ener = (_ein('qpir,pqrabc,baic->i', h2e_v, a22, h2e_v)
            + _ein('qpir,pqra,ai->i', h2e_v, a23, h1e_v) * 2.0
            + _ein('pi,pa,ai->i', h1e_v, a25, h1e_v))
    norm = (_ein('qpir,rpqbac,baic->i', h2e_v, dm3_h, h2e_v)
            + _ein('qpir,rpqa,ai->i', h2e_v, dm2_h, h1e_v) * 2.0
            + _ein('pi,pa,ai->i', h1e_v, dm1_h, h1e_v))
    nb, hb, db = _norm_to_energy_bar(norm, ener, -np.asarray(e_core), 1.0)

    v_bar = np.zeros_like(h2e_v)
    h1v_bar = np.zeros_like(h1e_v)
    a22_bar = np.zeros_like(a22)
    a23_bar = np.zeros_like(a23)
    a25_bar = np.zeros_like(a25)
    dm3h_bar = np.zeros_like(dm3_h)
    dm2h_bar = np.zeros_like(dm2_h)
    dm1h_bar = np.zeros_like(dm1_h)

    ops = (h2e_v, a22, h2e_v)
    v_bar += _ein_bar('qpir,pqrabc,baic->i', hb, ops, 0, h2e_v.shape)
    a22_bar += _ein_bar('qpir,pqrabc,baic->i', hb, ops, 1, a22.shape)
    v_bar += _ein_bar('qpir,pqrabc,baic->i', hb, ops, 2, h2e_v.shape)
    ops = (h2e_v, a23, h1e_v)
    v_bar += 2.0 * _ein_bar('qpir,pqra,ai->i', hb, ops, 0, h2e_v.shape)
    a23_bar += 2.0 * _ein_bar('qpir,pqra,ai->i', hb, ops, 1, a23.shape)
    h1v_bar += 2.0 * _ein_bar('qpir,pqra,ai->i', hb, ops, 2, h1e_v.shape)
    ops = (h1e_v, a25, h1e_v)
    h1v_bar += _ein_bar('pi,pa,ai->i', hb, ops, 0, h1e_v.shape)
    a25_bar += _ein_bar('pi,pa,ai->i', hb, ops, 1, a25.shape)
    h1v_bar += _ein_bar('pi,pa,ai->i', hb, ops, 2, h1e_v.shape)

    ops = (h2e_v, dm3_h, h2e_v)
    v_bar += _ein_bar('qpir,rpqbac,baic->i', nb, ops, 0, h2e_v.shape)
    dm3h_bar += _ein_bar('qpir,rpqbac,baic->i', nb, ops, 1, dm3_h.shape)
    v_bar += _ein_bar('qpir,rpqbac,baic->i', nb, ops, 2, h2e_v.shape)
    ops = (h2e_v, dm2_h, h1e_v)
    v_bar += 2.0 * _ein_bar('qpir,rpqa,ai->i', nb, ops, 0, h2e_v.shape)
    dm2h_bar += 2.0 * _ein_bar('qpir,rpqa,ai->i', nb, ops, 1, dm2_h.shape)
    h1v_bar += 2.0 * _ein_bar('qpir,rpqa,ai->i', nb, ops, 2, h1e_v.shape)
    ops = (h1e_v, dm1_h, h1e_v)
    h1v_bar += _ein_bar('pi,pa,ai->i', nb, ops, 0, h1e_v.shape)
    dm1h_bar += _ein_bar('pi,pa,ai->i', nb, ops, 1, dm1_h.shape)
    h1v_bar += _ein_bar('pi,pa,ai->i', nb, ops, 2, h1e_v.shape)

    # hole densities back to the particle densities
    bar.dm2 += 2.0 * _ein_bar('abef,cd->abcdef', dm3h_bar, (dm2, delta), 0, dm2.shape)
    bar.dm3 += -dm3h_bar.transpose(np.argsort((0, 1, 3, 2, 4, 5)))
    bar.dm1 += 2.0 * _ein_bar('ab,cd->abcd', dm2h_bar, (dm1, delta), 0, dm1.shape)
    bar.dm2 += -dm2h_bar.transpose(np.argsort((0, 1, 3, 2)))
    bar.dm1 += -dm1h_bar.T

    _a22_bar(h1e, h2e, dm2, dm3, f3ca, f3ac, a22_bar, bar)
    _a23_bar(h1e, h2e, dm1, dm2, dm3, a23_bar, bar)
    _a25_bar(h1e, h2e, dm1, dm2, a25_bar, bar)
    return v_bar, h1v_bar, -db


def _Sijrs_bar(e_core, e_virt, g_cvcv):
    """Adjoint of ``nevpt2_sc._Sijrs`` (the MP2-like fully-external subspace).

    Forward, per occupied ``i``:
        gi    = g_cvcv[i].transpose(1,0,2)                   # (j,a,b)
        djba  = (eia[:,None] + eia[i][None,:]).ravel()       # e_i+e_j-e_a-e_b
        t2i   = gi / djba
        theta = 2 gi - gi.transpose(0,2,1)
        e    += sum(t2i * theta)
    """
    ncore = len(e_core)
    nvirt = len(e_virt)
    g_bar = np.zeros_like(g_cvcv)
    ec_bar = np.zeros(ncore)
    ev_bar = np.zeros(nvirt)
    if ncore == 0 or nvirt == 0:
        return g_bar, ec_bar, ev_bar
    eia = e_core[:, None] - e_virt[None, :]
    for i in range(ncore):
        gi = g_cvcv[i].transpose(1, 0, 2)
        d = (eia.reshape(-1, 1) + eia[i].reshape(1, -1)).ravel().reshape(
            ncore, nvirt, nvirt)
        t2i = gi / d
        theta = gi * 2 - gi.transpose(0, 2, 1)
        # e = sum(t2i * theta)
        t2_bar = theta
        theta_bar = t2i
        gi_bar = t2_bar / d
        gi_bar += 2.0 * theta_bar - theta_bar.transpose(0, 2, 1)
        d_bar = -t2_bar * gi / (d * d)
        # d[j,a,b] = e_i + e_j - e_a - e_b
        ec_bar[i] += float(d_bar.sum())
        ec_bar += d_bar.sum(axis=(1, 2))
        ev_bar -= d_bar.sum(axis=(0, 2))
        ev_bar -= d_bar.sum(axis=(0, 1))
        g_bar[i] += gi_bar.transpose(1, 0, 2)
    return g_bar, ec_bar, ev_bar


def _Sijr_bar(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt, bar):
    """Adjoint of ``nevpt2_sc._Sijr``."""
    from oqp.library.nevpt2_sc import _a3
    ncore = len(e_core)
    v_bar = np.zeros_like(h2e_v)
    ec_bar = np.zeros(ncore)
    ev_bar = np.zeros(len(e_virt))
    if ncore == 0 or len(e_virt) == 0:
        return v_bar, ec_bar, ev_bar
    hdm1 = _hdm1(dm1)
    a3 = _a3(h1e, h2e, dm1, dm2, hdm1)
    ci_diag = np.diag_indices(ncore)
    ci_triu = np.triu_indices(ncore)

    def _pack(mat_from_pa):
        x = (2.0 * _ein('rpji,raji,pa->rji', h2e_v, h2e_v, mat_from_pa)
             - 1.0 * _ein('rpji,raij,pa->rji', h2e_v, h2e_v, mat_from_pa))
        x = x + x.transpose(0, 2, 1)
        x[:, ci_diag[0], ci_diag[1]] *= 0.5
        return x

    norm_full = _pack(hdm1)
    h_full = _pack(a3)
    diff = e_virt[:, None, None] - e_core[None, :, None] - e_core[None, None, :]
    sel = (slice(None), ci_triu[0], ci_triu[1])
    nb_p, hb_p, db_p = _norm_to_energy_bar(
        norm_full[sel], h_full[sel], diff[sel], 1.0)

    def _unpack(packed_bar):
        full = np.zeros_like(norm_full)
        full[sel] = packed_bar
        full[:, ci_diag[0], ci_diag[1]] *= 0.5
        full = full + full.transpose(0, 2, 1)
        return full

    nb_full = _unpack(nb_p)
    hb_full = _unpack(hb_p)
    db_full = np.zeros_like(diff)
    db_full[sel] = db_p
    ev_bar += db_full.sum(axis=(1, 2))
    ec_bar -= db_full.sum(axis=(0, 2))
    ec_bar -= db_full.sum(axis=(0, 1))

    hdm1_bar = np.zeros_like(hdm1)
    a3_bar = np.zeros_like(a3)
    for full_bar, mat, mat_bar in ((nb_full, hdm1, hdm1_bar),
                                   (hb_full, a3, a3_bar)):
        ops = (h2e_v, h2e_v, mat)
        v_bar += 2.0 * _ein_bar('rpji,raji,pa->rji', full_bar, ops, 0, h2e_v.shape)
        v_bar += 2.0 * _ein_bar('rpji,raji,pa->rji', full_bar, ops, 1, h2e_v.shape)
        mat_bar += 2.0 * _ein_bar('rpji,raji,pa->rji', full_bar, ops, 2, mat.shape)
        v_bar -= _ein_bar('rpji,raij,pa->rji', full_bar, ops, 0, h2e_v.shape)
        v_bar -= _ein_bar('rpji,raij,pa->rji', full_bar, ops, 1, h2e_v.shape)
        mat_bar -= _ein_bar('rpji,raij,pa->rji', full_bar, ops, 2, mat.shape)

    _a3_bar(h1e, h2e, dm1, dm2, hdm1, a3_bar, bar, hdm1_bar)
    _hdm1_bar(hdm1_bar, bar)
    return v_bar, ec_bar, ev_bar


def _Srsi_bar(dm1, dm2, h1e, h2e, h2e_v, e_core, e_virt, bar):
    """Adjoint of ``nevpt2_sc._Srsi``."""
    from oqp.library.nevpt2_sc import _k27
    nvirt = len(e_virt)
    v_bar = np.zeros_like(h2e_v)
    ec_bar = np.zeros(len(e_core))
    ev_bar = np.zeros(nvirt)
    if len(e_core) == 0 or nvirt == 0:
        return v_bar, ec_bar, ev_bar
    k27 = _k27(h1e, h2e, dm1, dm2)
    vi_diag = np.diag_indices(nvirt)
    vi_triu = np.triu_indices(nvirt)

    def _pack(mat):
        x = (2.0 * _ein('rsip,rsia,pa->rsi', h2e_v, h2e_v, mat)
             - 1.0 * _ein('rsip,sria,pa->rsi', h2e_v, h2e_v, mat))
        x = x + x.transpose(1, 0, 2)
        x[vi_diag] *= 0.5
        return x

    norm_full = _pack(dm1)
    h_full = _pack(k27)
    diff = e_virt[:, None, None] + e_virt[None, :, None] - e_core[None, None, :]
    nb_p, hb_p, db_p = _norm_to_energy_bar(
        norm_full[vi_triu], h_full[vi_triu], diff[vi_triu], 1.0)

    def _unpack(packed_bar):
        full = np.zeros_like(norm_full)
        full[vi_triu] = packed_bar
        full[vi_diag] *= 0.5
        full = full + full.transpose(1, 0, 2)
        return full

    nb_full = _unpack(nb_p)
    hb_full = _unpack(hb_p)
    db_full = np.zeros_like(diff)
    db_full[vi_triu] = db_p
    ev_bar += db_full.sum(axis=(1, 2))
    ev_bar += db_full.sum(axis=(0, 2))
    ec_bar -= db_full.sum(axis=(0, 1))

    k27_bar = np.zeros_like(k27)
    for full_bar, mat, mat_bar in ((nb_full, dm1, None), (hb_full, k27, k27_bar)):
        ops = (h2e_v, h2e_v, mat)
        v_bar += 2.0 * _ein_bar('rsip,rsia,pa->rsi', full_bar, ops, 0, h2e_v.shape)
        v_bar += 2.0 * _ein_bar('rsip,rsia,pa->rsi', full_bar, ops, 1, h2e_v.shape)
        contrib = 2.0 * _ein_bar('rsip,rsia,pa->rsi', full_bar, ops, 2, mat.shape)
        v_bar -= _ein_bar('rsip,sria,pa->rsi', full_bar, ops, 0, h2e_v.shape)
        v_bar -= _ein_bar('rsip,sria,pa->rsi', full_bar, ops, 1, h2e_v.shape)
        contrib -= _ein_bar('rsip,sria,pa->rsi', full_bar, ops, 2, mat.shape)
        if mat_bar is None:
            bar.dm1 += contrib
        else:
            mat_bar += contrib

    _k27_bar(h1e, h2e, dm1, dm2, k27_bar, bar)
    return v_bar, ec_bar, ev_bar


def _Srs_bar(dm1, dm2, dm3, h1e, h2e, h2e_v, e_virt, bar):
    """Adjoint of ``nevpt2_sc._Srs``."""
    from oqp.library.nevpt2_sc import _a7
    v_bar = np.zeros_like(h2e_v)
    ev_bar = np.zeros(len(e_virt))
    if len(e_virt) == 0:
        return v_bar, ev_bar
    rm2, a7 = _a7(h1e, h2e, dm1, dm2, dm3)
    norm = 0.5 * _ein('rsqp,rsba,pqba->rs', h2e_v, h2e_v, rm2)
    h = 0.5 * _ein('rsqp,rsba,pqab->rs', h2e_v, h2e_v, a7)
    diff = e_virt[:, None] + e_virt[None, :]
    nb, hb, db = _norm_to_energy_bar(norm, h, diff, 1.0)
    ev_bar += db.sum(axis=1) + db.sum(axis=0)

    rm2_bar = np.zeros_like(rm2)
    a7_bar = np.zeros_like(a7)
    ops = (h2e_v, h2e_v, rm2)
    v_bar += 0.5 * _ein_bar('rsqp,rsba,pqba->rs', nb, ops, 0, h2e_v.shape)
    v_bar += 0.5 * _ein_bar('rsqp,rsba,pqba->rs', nb, ops, 1, h2e_v.shape)
    rm2_bar += 0.5 * _ein_bar('rsqp,rsba,pqba->rs', nb, ops, 2, rm2.shape)
    ops = (h2e_v, h2e_v, a7)
    v_bar += 0.5 * _ein_bar('rsqp,rsba,pqab->rs', hb, ops, 0, h2e_v.shape)
    v_bar += 0.5 * _ein_bar('rsqp,rsba,pqab->rs', hb, ops, 1, h2e_v.shape)
    a7_bar += 0.5 * _ein_bar('rsqp,rsba,pqab->rs', hb, ops, 2, a7.shape)

    _a7_bar(h1e, h2e, dm1, dm2, dm3, rm2, rm2_bar, a7_bar, bar)
    return v_bar, ev_bar


def _Sij_bar(dm1, dm2, dm3, h1e, h2e, h2e_v, e_core, bar):
    """Adjoint of ``nevpt2_sc._Sij``."""
    from oqp.library.nevpt2_sc import _a9
    v_bar = np.zeros_like(h2e_v)
    ec_bar = np.zeros(len(e_core))
    if len(e_core) == 0:
        return v_bar, ec_bar
    hdm1 = _hdm1(dm1)
    hdm2 = _hdm2(dm1, dm2)
    hdm3 = _hdm3(dm1, dm2, dm3, hdm1, hdm2)
    a9 = _a9(h1e, h2e, hdm1, hdm2, hdm3)
    norm = 0.5 * _ein('qpij,baij,pqab->ij', h2e_v, h2e_v, hdm2)
    h = 0.5 * _ein('qpij,baij,pqab->ij', h2e_v, h2e_v, a9)
    diff = e_core[:, None] + e_core[None, :]
    nb, hb, db = _norm_to_energy_bar(norm, h, -diff, 1.0)
    ec_bar -= db.sum(axis=1) + db.sum(axis=0)

    hdm2_bar = np.zeros_like(hdm2)
    hdm3_bar = np.zeros_like(hdm3)
    a9_bar = np.zeros_like(a9)
    ops = (h2e_v, h2e_v, hdm2)
    v_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', nb, ops, 0, h2e_v.shape)
    v_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', nb, ops, 1, h2e_v.shape)
    hdm2_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', nb, ops, 2, hdm2.shape)
    ops = (h2e_v, h2e_v, a9)
    v_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', hb, ops, 0, h2e_v.shape)
    v_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', hb, ops, 1, h2e_v.shape)
    a9_bar += 0.5 * _ein_bar('qpij,baij,pqab->ij', hb, ops, 2, a9.shape)

    _a9_bar(h1e, h2e, hdm2, hdm3, a9_bar, hdm2_bar, hdm3_bar, bar)
    _hdm3_bar(dm1, dm2, dm3, hdm2, hdm3_bar, bar, hdm2_bar)
    _hdm2_bar(dm1, dm2, hdm2_bar, bar)
    return v_bar, ec_bar


def _Sir_bar(dm1, dm2, dm3, h1e, h2e, h1e_v, h2e_v1, h2e_v2, e_core, e_virt, bar):
    """Adjoint of ``nevpt2_sc._Sir``."""
    from oqp.library.nevpt2_sc import _a12, _a13
    v1_bar = np.zeros_like(h2e_v1)
    v2_bar = np.zeros_like(h2e_v2)
    h1v_bar = np.zeros_like(h1e_v)
    ec_bar = np.zeros(len(e_core))
    ev_bar = np.zeros(len(e_virt))
    if len(e_core) == 0 or len(e_virt) == 0:
        return v1_bar, v2_bar, h1v_bar, ec_bar, ev_bar

    norm_terms = (
        ('rpiq,raib,qpab->ir', (h2e_v1, h2e_v1, dm2), +2.0),
        ('rpiq,rabi,qpab->ir', (h2e_v1, h2e_v2, dm2), -1.0),
        ('rpqi,raib,qpab->ir', (h2e_v2, h2e_v1, dm2), -1.0),
        ('raqi,rabi,qb->ir', (h2e_v2, h2e_v2, dm1), +2.0),
        ('rpqi,rabi,qbap->ir', (h2e_v2, h2e_v2, dm2), -1.0),
        ('rpqi,raai,qp->ir', (h2e_v2, h2e_v2, dm1), +1.0),
        ('rpiq,ri,qp->ir', (h2e_v1, h1e_v, dm1), +4.0),
        ('rpqi,ri,qp->ir', (h2e_v2, h1e_v, dm1), -2.0),
        ('ri,ri->ir', (h1e_v, h1e_v), +2.0),
    )
    norm = sum(sign * _ein(sub, *ops) for sub, ops, sign in norm_terms)
    a12 = _a12(h1e, h2e, dm1, dm2, dm3)
    a13 = _a13(h1e, h2e, dm1, dm2, dm3)
    h_terms = (
        ('rpiq,raib,pqab->ir', (h2e_v1, h2e_v1, a12), +2.0),
        ('rpiq,rabi,pqab->ir', (h2e_v1, h2e_v2, a12), -1.0),
        ('rpqi,raib,pqab->ir', (h2e_v2, h2e_v1, a12), -1.0),
        ('rpqi,rabi,pqab->ir', (h2e_v2, h2e_v2, a13), +1.0),
    )
    h = sum(sign * _ein(sub, *ops) for sub, ops, sign in h_terms)
    diff = e_core[:, None] - e_virt[None, :]
    nb, hb, db = _norm_to_energy_bar(norm, h, -diff, 1.0)
    ec_bar -= db.sum(axis=1)
    ev_bar += db.sum(axis=0)

    a12_bar = np.zeros_like(a12)
    a13_bar = np.zeros_like(a13)

    def _route(op, contrib):
        nonlocal v1_bar, v2_bar, h1v_bar, a12_bar, a13_bar
        if op is h2e_v1:
            v1_bar = v1_bar + contrib
        elif op is h2e_v2:
            v2_bar = v2_bar + contrib
        elif op is h1e_v:
            h1v_bar = h1v_bar + contrib
        elif op is dm1:
            bar.dm1 += contrib
        elif op is dm2:
            bar.dm2 += contrib
        elif op is a12:
            a12_bar = a12_bar + contrib
        elif op is a13:
            a13_bar = a13_bar + contrib
        else:  # pragma: no cover - defensive
            raise AssertionError("unexpected Sir operand")

    for sub, ops, sign in norm_terms:
        for k, op in enumerate(ops):
            _route(op, sign * _ein_bar(sub, nb, ops, k, op.shape))
    for sub, ops, sign in h_terms:
        for k, op in enumerate(ops):
            _route(op, sign * _ein_bar(sub, hb, ops, k, op.shape))

    _a12_bar(h1e, h2e, dm2, dm3, a12_bar, bar)
    _a13_bar(h1e, h2e, dm1, dm2, dm3, a13_bar, bar)
    return v1_bar, v2_bar, h1v_bar, ec_bar, ev_bar


# ----------------------------------------------------------- block assembly
def _blocks_bar(h1e_mo, eri_mo, ncore, nact, block_bars, bar):
    """Push the per-subspace integral adjoints back through ``nevpt2_sc._blocks``.

    ``_blocks`` builds every subspace tensor from the core-dressed matrix
    ``heff = h1e_mo + vhf`` and permuted slices of the chemist ``eri_mo``; the
    adjoints below invert exactly those slices and permutations.
    """
    nbf = h1e_mo.shape[0]
    nocc = ncore + nact
    C = slice(0, ncore)
    A = slice(ncore, nocc)
    V = slice(nocc, nbf)

    hbar = np.zeros((nbf, nbf))
    gbar = np.zeros((nbf,) * 4)
    heff_bar = np.zeros((nbf, nbf))

    def unphys(block_bar):
        """Inverse of ``block.transpose(0,2,1,3)``."""
        return block_bar.transpose(0, 2, 1, 3)

    # active blocks used directly as h1e/h2e
    heff_bar[A, A] += bar.h1e
    gbar[A, A, A, A] += bar.h2e.transpose(np.argsort((0, 2, 1, 3)))

    # Sr: v = phys(eri[V,A,A,A]);  h1v = heff[V,A] - ein('mbbn->mn', v)
    v_bar, h1v_bar = block_bars['Sr']
    heff_bar[V, A] += h1v_bar
    v_bar = v_bar - _ein_bar('mbbn->mn', h1v_bar,
                             (np.zeros_like(v_bar),), 0, v_bar.shape)
    gbar[V, A, A, A] += unphys(v_bar)

    # Si: v = phys(eri[A,C,A,A]);  h1v = heff[A,C]
    v_bar, h1v_bar = block_bars['Si']
    heff_bar[A, C] += h1v_bar
    gbar[A, C, A, A] += unphys(v_bar)

    gbar[C, V, C, V] += block_bars['Sijrs']
    gbar[V, C, A, C] += unphys(block_bars['Sijr'])
    gbar[V, C, V, A] += unphys(block_bars['Srsi'])
    gbar[V, A, V, A] += unphys(block_bars['Srs'])
    gbar[A, C, A, C] += unphys(block_bars['Sij'])

    v1_bar, v2_bar, h1v_bar = block_bars['Sir']
    gbar[V, C, A, A] += unphys(v1_bar)
    gbar[V, A, A, C] += unphys(v2_bar)
    heff_bar[V, C] += h1v_bar

    # heff = h1e_mo + 2 ein('pqii->pq', eri[:,:,C,C]) - ein('piiq->pq', eri[:,C,C,:])
    hbar += heff_bar
    ncore_idx = np.arange(ncore)
    gbar[:, :, ncore_idx, ncore_idx] += 2.0 * heff_bar[:, :, None]
    gbar[:, ncore_idx, ncore_idx, :] -= heff_bar[:, None, :]
    return hbar, gbar


# --------------------------------------------------------------------- driver
def sc_nevpt2_energy_adjoints(h1e_mo, eri_mo, eps, ncore, nact, active_nelec,
                              ci_vector, dms=None):
    """Exact derivatives of the SC-NEVPT2 correlation energy.

    Returns ``(e2, comp, hbar, gbar, epsbar, dmbars)`` where ``dmbars`` is the
    tuple ``(dm1bar, dm2bar, dm3bar, dm4bar)`` of derivatives with respect to
    the active spin-free density matrices.  The caller converts ``dmbars`` into
    a CI-coefficient derivative, because that conversion needs the determinant
    basis the density matrices were built from.

    ``hbar``/``gbar`` are the derivatives with respect to the semicanonical-basis
    MO integrals at FIXED ``eps``; ``epsbar`` is the derivative with respect to
    the orbital energies at fixed integrals.  Keeping the two separate is what
    lets the caller apply the generalized-Fock chain rule once, for both the
    diagonal (``eps``) and the off-diagonal (semicanonical-constraint) response.
    """
    from oqp.library.nevpt2_sc import (
        _Sr, _Si, _Sijrs, _Sijr, _Srsi, _Srs, _Sij, _Sir)

    nbf = h1e_mo.shape[0]
    if dms is None:
        dms = make_rdms(ci_vector, nact, active_nelec, upto=4)
    dm1, dm2, dm3, dm4 = dms
    B = _blocks(h1e_mo, eri_mo, ncore, nact, eps)
    h1e, h2e = B['h1e'], B['h2e']
    ec, ev = B['e_core'], B['e_virt']
    f3 = _f3ca_f3ac(h2e, dm4)

    bar = _Bar(nact, nbf)
    epsbar = np.zeros(nbf)
    nocc = ncore + nact
    comp = {}
    block_bars = {}

    v, h1v = B['Sr']
    comp['Sr'] = _Sr(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ev, f3=f3)[1]
    vb, hb, evb = _Sr_bar(dm1, dm2, dm3, h1e, h2e, h1v, v, ev, f3, bar)
    block_bars['Sr'] = (vb, hb)
    epsbar[nocc:] += evb

    v, h1v = B['Si']
    comp['Si'] = _Si(dm1, dm2, dm3, dm4, h1e, h2e, h1v, v, ec, f3=f3)[1]
    vb, hb, ecb = _Si_bar(dm1, dm2, dm3, h1e, h2e, h1v, v, ec, f3, bar)
    block_bars['Si'] = (vb, hb)
    epsbar[:ncore] += ecb

    comp['Sijrs'] = _Sijrs(ec, ev, B['Sijrs'])[1]
    gb, ecb, evb = _Sijrs_bar(ec, ev, B['Sijrs'])
    block_bars['Sijrs'] = gb
    epsbar[:ncore] += ecb
    epsbar[nocc:] += evb

    comp['Sijr'] = _Sijr(dm1, dm2, h1e, h2e, B['Sijr'], ec, ev)[1]
    vb, ecb, evb = _Sijr_bar(dm1, dm2, h1e, h2e, B['Sijr'], ec, ev, bar)
    block_bars['Sijr'] = vb
    epsbar[:ncore] += ecb
    epsbar[nocc:] += evb

    comp['Srsi'] = _Srsi(dm1, dm2, h1e, h2e, B['Srsi'], ec, ev)[1]
    vb, ecb, evb = _Srsi_bar(dm1, dm2, h1e, h2e, B['Srsi'], ec, ev, bar)
    block_bars['Srsi'] = vb
    epsbar[:ncore] += ecb
    epsbar[nocc:] += evb

    comp['Srs'] = _Srs(dm1, dm2, dm3, h1e, h2e, B['Srs'], ev)[1]
    vb, evb = _Srs_bar(dm1, dm2, dm3, h1e, h2e, B['Srs'], ev, bar)
    block_bars['Srs'] = vb
    epsbar[nocc:] += evb

    comp['Sij'] = _Sij(dm1, dm2, dm3, h1e, h2e, B['Sij'], ec)[1]
    vb, ecb = _Sij_bar(dm1, dm2, dm3, h1e, h2e, B['Sij'], ec, bar)
    block_bars['Sij'] = vb
    epsbar[:ncore] += ecb

    v1, v2, h1v = B['Sir']
    comp['Sir'] = _Sir(dm1, dm2, dm3, h1e, h2e, h1v, v1, v2, ec, ev)[1]
    v1b, v2b, hb, ecb, evb = _Sir_bar(dm1, dm2, dm3, h1e, h2e, h1v, v1, v2,
                                      ec, ev, bar)
    block_bars['Sir'] = (v1b, v2b, hb)
    epsbar[:ncore] += ecb
    epsbar[nocc:] += evb

    # f3ca/f3ac were shared by Sr and Si; unwind them once, after both.
    _f3ca_f3ac_bar(h2e, dm4, bar.f3ca, bar.f3ac, bar)
    hbar, gbar = _blocks_bar(h1e_mo, eri_mo, ncore, nact, block_bars, bar)

    e2 = float(sum(comp.values()))
    return e2, comp, hbar, gbar, epsbar, (bar.dm1, bar.dm2, bar.dm3, bar.dm4)
