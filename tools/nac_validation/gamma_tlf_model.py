"""Closed-form model of the TLF orbital kernel gamma~.

Replicates compute_states_overlap's 7-term assembly in numpy with
LINEARIZED determinant blocks,

    s_ij(i1,i2) = delta + sg_ij * theta(i2,i1)   (i1 != i2)
    s_ab(a1,a2) = delta + sg_ab * theta(a1,a2)   (a1 != a2)
    s_ia(i ,a ) = delta_socc + sg_ia * theta(i,a)

(diagonal first-order parts are symmetric and drop from the
antisymmetrized kernel), then extracts the kernel analytically as the
coefficient of theta_pq and validates against the oracle tensor
/tmp/nactest/tlf_gamma.npz (same-process amplitudes included there).
The three signs are scanned to find the convention matching the oracle.
"""
import numpy as np
import itertools

RS = 1.0 / np.sqrt(2.0)

d = np.load('/tmp/nactest/tlf_gamma.npz')
gam_o = d['gam']                  # (nbf,nbf,nst,nst) oracle kernel
noca, nocb, nbf = int(d['noca']), int(d['nocb']), int(d['nbf'])
X0_raw = d['X0']
nst = gam_o.shape[2]
nvirb = nbf - nocb
nij = noca * nvirb


def unfold(bv, st):
    """folded -> (noca, nvirb) amplitude matrix, mult=1 conventions."""
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    x = np.zeros((noca, nvirb))
    for i in range(1, noca + 1):
        for jj in range(nocb + 1, nbf + 1):
            ij = (jj - nocb - 1) * noca + i
            if ij == ijlr1:
                x[i - 1, jj - nocb - 1] = bv[ijlr1 - 1, st - 1] * RS
            elif ij == ijlr2:
                x[i - 1, jj - nocb - 1] = -bv[ijlr1 - 1, st - 1] * RS
            else:
                x[i - 1, jj - nocb - 1] = bv[ij - 1, st - 1]
    return x


bvF = X0_raw.reshape(-1).reshape((nst, nij)).T
C = [unfold(bvF, s + 1) for s in range(nst)]   # unfolded amplitudes per state

# masks in the (noca, nvirb) amplitude grid
iocc = np.arange(noca)
avir = np.arange(nvirb)
OO = (iocc[:, None] >= nocb) & (avir[None, :] < 2)     # the 4 socc^2 slots
GEN = ~OO

# socc indicator: occ index i in socc <-> i >= nocb; virtb index a in socc <-> a < 2
# s_ia(i,a) at coincidence: delta for i = global(a) i.e. i = nocb + a, a < 2
S_IA0 = np.zeros((noca, nvirb))
for a in range(2):
    S_IA0[nocb + a, a] = 1.0


def kernel_pair(cI, cJ, sg_ij, sg_ab, sg_ia):
    """Analytic first-order kernel G[p,q] (nbf x nbf, MO indices) such that
    S_IJ^(1) = sum_pq G_pq theta_pq, for the 7-term TLF assembly."""
    G = np.zeros((nbf, nbf))
    co = cI * GEN     # x_old with OO removed (for the generic-term sums)
    cn = cJ * GEN

    # ---- term A: sum_{i,a} alpham(i,a) betam(i,a)
    # alpham = sum_b co(i,b) s_ab(b,a); betam = sum_j cn(j,a) s_ij(i,j)
    # first order: d s_ab(b,a) = sg_ab theta(B,Aglob); d s_ij(i,j) = sg_ij theta(j,i)
    # (A) d from s_ab, s_ij = I:  sum_{i,a,b} co(i,b) cn(i,a) sg_ab th(b,a)
    T = co.T @ cn                            # T[b,a] = sum_i co(i,b) cn(i,a)
    G[nocb:, nocb:] += sg_ab * T
    # (B) d from s_ij, s_ab = I:  sum_{i,j,a} co(i,a) cn(j,a) sg_ij th(j,i)
    T2 = cn @ co.T                           # T2[j,i] = sum_a cn(j,a) co(i,a)
    G[:noca, :noca] += sg_ij * T2

    # ---- term B: sum_{i1,i2} gammam(i1,i2) deltam(i1,i2)
    # gammam = sum_b co(i1,b) s_ia(i2,b); deltam = sum_a cn(i2,a) s_ia(i1,a)
    # zeroth: gammam0(i1,i2) = sum_b co(i1,b) S_IA0(i2,b), etc.
    g0 = co @ S_IA0.T                        # g0[i1,i2]
    d0 = cn @ S_IA0.T                        # d0[i2,i1] -> careful: deltam(i1,i2)=sum_a cn(i2,a) S_IA0(i1,a)
    d0 = (cn @ S_IA0.T).T                    # d0[i1,i2] = sum_a cn(i2,a) S_IA0(i1,a)
    # first order from gammam: d gammam(i1,i2) = sg_ia sum_b co(i1,b) th(i2, Bglob)
    # contracted with deltam0:  sum_{i1,i2,b} co(i1,b) d0[i1,i2] sg_ia th(i2,b+nocb)
    W = co.T @ d0                            # W[b, i2] = sum_i1 co(i1,b) d0[i1,i2]
    G[:noca, nocb:] += sg_ia * W.T
    # first order from deltam: d deltam(i1,i2) = sg_ia sum_a cn(i2,a) th(i1, Aglob)
    # contracted with gammam0: sum co(i1,b)S_IA0(i2,b) * cn(i2,a) sg_ia th(i1,a+nocb)
    W2 = g0 @ cn                             # W2[i1, a] = sum_i2 g0[i1,i2] cn(i2,a)
    G[:noca, nocb:] += sg_ia * W2

    # ---- term C: OO x OO, s_ij s_ab only, weight 1
    # sum_{p,q,r,s in socc} cI(p,qv) s_ij(p,r) cJ(r,sv) s_ab(qv,sv)
    pq = slice(nocb, noca)                   # socc occ indices
    qv = slice(0, 2)                         # socc virtb indices
    cIo = cI[pq, qv]                         # 2x2 OO amplitudes of I
    cJo = cJ[pq, qv]
    # d s_ij(p,r): sg_ij th(r,p);  d s_ab(q,s): sg_ab th(qg, sg)
    # from s_ij: sum_{p,r} [sum_q cIo(p,q) cJo(r,q)] sg_ij th(r,p)
    TC = cJo @ cIo.T                         # TC[r,p]
    G[nocb:noca, nocb:noca] += sg_ij * TC
    # from s_ab: sum_{q,s} [sum_p cIo(p,q) cJo(p,s)] sg_ab th(qg, sg)
    TC2 = cIo.T @ cJo                        # TC2[q,s]
    G[nocb:noca, nocb:noca] += sg_ab * TC2

    # ---- terms D and E: OO x generic cross, weight 1/sqrt2,
    #      products s_ij s_ab + s_ia s_ia
    for (cA, cB, sw) in ((cI, cJ, +1), (cJ, cI, -1)):
        # sw: the antisymmetrization I<->J swaps below handled at the end via
        # the (S_IJ - S_JI) construction; here build S_IJ only -> sw unused.
        pass
    # D: old OO (cI), new generic (cJ)
    cIo_full = np.zeros_like(cI); cIo_full[pq, qv] = cI[pq, qv]
    cJg = cJ * GEN
    # direct: d[s_ij(p,r) s_ab(q,s)]: from s_ij: th(r,p)*delta_qs... at coincidence
    #   s_ab(q,s)=delta(qg,sg) -> r free over occ, s.t. s=q
    # sum_{p,q in socc; r} cIo(p,q) cJg(r,q) sg_ij th(r,p) / sqrt2
    TD = cJg[:, 0:2] @ cIo_full[pq, qv].T    # TD[r, p_rel] = sum_q cJg(r,q) cIo(p,q)
    G[:noca, nocb:noca] += RS * sg_ij * TD
    #   from s_ab: delta_pr -> sum_{p,q in socc; s} cIo(p,q) cJg(p,s) sg_ab th(qg, sg)
    TD2 = cIo_full[pq, qv].T @ cJg[pq, :]    # TD2[q_rel, s]
    G[nocb:noca, nocb:].flat  # noop for clarity
    G[nocb:nocb + 2, nocb:] += RS * sg_ab * TD2
    # exchange: d[s_ia(p,s) s_ia(r,q)]:
    #   from first factor: s_ia(r,q)|0 = delta(r = qg, both socc) ->
    #     sum_{p,q socc; s} cIo(p,q) cJg(qg, s) sg_ia th(p, sg)
    TE = cIo_full[pq, qv] @ cJg[pq, :]       # TE[p_rel, s] = sum_q cIo(p,q) cJg(nocb+q, s)
    G[nocb:noca, nocb:] += RS * sg_ia * TE
    #   from second factor: s_ia(p,s)|0 = delta(p = sg) ->
    #     sum_{p,q socc; r} cIo(p,q) cJg(r, p_rel) sg_ia th(r, qg)
    TE2 = cJg[:, 0:2] @ cIo_full[pq, qv]     # TE2[r, q_rel] = sum_p cJg(r,p_rel) cIo(p,q)
    G[:noca, nocb:nocb + 2] += RS * sg_ia * TE2
    # E: old generic (cI), new OO (cJ) -- mirror of D
    cJo_full = np.zeros_like(cJ); cJo_full[pq, qv] = cJ[pq, qv]
    cIg = cI * GEN
    # direct from s_ij: sum_{r,s socc; p} cIg(p,s) cJo(r,s) sg_ij th(r,p)
    TF = cJo_full[pq, qv] @ cIg[:, 0:2].T    # TF[r_rel, p]
    G[nocb:noca, :noca] += RS * sg_ij * TF
    # direct from s_ab: sum_{r,s socc; q} cIg(r, q)... delta_pr with p generic occ=r in socc
    TF2 = cIg[pq, :].T @ cJo_full[pq, qv]    # TF2[q, s_rel]
    G[nocb:, nocb:nocb + 2] += RS * sg_ab * TF2
    # exchange first factor: s_ia(r,q)|0 = delta(r=qg) ->
    #   sum_{r,s socc; ... } cIg(p, s')... mirror of TE/TE2 with roles swapped:
    #   sum_{r,s socc; p} cIg(p, s) cJo(r, s)... d s_ia(p, sg)... careful:
    #   term E product: cIg(p,q') cJo(r,s) [s_ij(p,r) s_ab(q',s) + s_ia(p,s) s_ia(r,q')]
    #   exchange from d s_ia(p,s): s_ia(r,q')|0 = delta(rg = q') ->
    #     sum_{r,s socc; p} cIg(p, r_rel) cJo(r,s) sg_ia th(p, sg)
    TG = cIg[:, 0:2] @ cJo_full[pq, qv]      # TG[p, s_rel]
    G[:noca, nocb:nocb + 2] += RS * sg_ia * TG
    #   exchange from d s_ia(r,q'): s_ia(p,s)|0 = delta(p = sg) ->
    #     sum_{r,s socc; q'} cIg(sg, q') cJo(r,s) sg_ia th(r, q'g)
    TG2 = cJo_full[pq, qv] @ cIg[pq, :]      # TG2[r_rel, q']
    G[nocb:noca, nocb:] += RS * sg_ia * TG2
    return G


def kernel_antisym(I, J, sgs):
    GIJ = kernel_pair(C[I - 1], C[J - 1], *sgs)
    GJI = kernel_pair(C[J - 1], C[I - 1], *sgs)
    return GIJ - GJI


best = None
for sgs in itertools.product((1, -1), repeat=3):
    err = 0.0
    for (I, J) in ((1, 2), (1, 3), (2, 3)):
        Gm = kernel_antisym(I, J, sgs)
        Go = gam_o[:, :, I - 1, J - 1]
        err += min(np.linalg.norm(Gm - Go), np.linalg.norm(Gm + Go))
    if best is None or err < best[0]:
        best = (err, sgs)
print('best signs (s_ij, s_ab, s_ia):', best[1], ' total err =', best[0])

sgs = best[1]
for (I, J) in ((1, 2), (1, 3), (2, 3)):
    Gm = kernel_antisym(I, J, sgs)
    Go = gam_o[:, :, I - 1, J - 1]
    s = 1.0 if np.linalg.norm(Gm - Go) < np.linalg.norm(Gm + Go) else -1.0
    print(f'({I},{J}): |model - oracle| = {np.linalg.norm(s*Gm - Go):.6f}   '
          f'|oracle| = {np.linalg.norm(Go):.6f}')
