"""C1 Delta-structure fit on the ethylene v7o npz. Candidates:
gamma:Sk (pair/transpose), transport tA/tB (occ/virt Sk-actions via the
offline formula-kernel G_met), Sx-actions, d_num, and per-block Delta
decompositions. C1 kills the symmetry degeneracies.
Run: python eth_delta_fit.py <eth_v7o.npz> <eth_v7h-like-not-needed> <dnum.npz> (uses npz only)
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'repo_nac', 'tools', 'nac_lagrangian'))
import nac_formula_kernel as FK

d = np.load(sys.argv[1])
ctxf = np.load(sys.argv[2])
dn = np.load(sys.argv[3])
dcv = dn['dcv' if 'dcv' in dn.files else dn.files[0]]
U1, U2, w1, w2 = d['Ux1'], d['Ux2'], d['w1'], d['w2']
noca, nocb = int(ctxf['noca']), int(ctxf['nocb'])
gam0 = ctxf['gam']; Sk_an = ctxf['Sk_an']; Sx = ctxf['Sx']
Xf0 = ctxf['Xf']; wprobe = ctxf['wprobe']; c0 = int(ctxf['c0'])
nstate_ = 3
# phase alignment ctx -> v7o via the w-probe
s = np.zeros(nstate_)
for J in range(nstate_):
    s[J] = np.sign(float(np.dot(w1[c0, J, :], wprobe[J, :])))
print('phase alignment s_J =', s)
Xf = Xf0 * s[None, :]
gam = gam0.copy()
for I in range(nstate_):
    for J in range(nstate_):
        gam[I, J] = gam0[I, J] * s[I] * s[J]
ncoord, nbf = U1.shape[0], U1.shape[1]
nstate = 3
nvirb = nbf - nocb
nij = noca * nvirb
RS = 1.0 / np.sqrt(2.0)
noc = noca - 1
ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
ijlr2 = (noca - nocb - 1) * noca + noca

def unfold_vec(v):
    x = np.zeros((noca, nvirb))
    for i in range(1, noca + 1):
        for jj in range(nocb + 1, nbf + 1):
            ij = (jj - nocb - 1) * noca + i
            if ij == ijlr1:
                x[i - 1, jj - nocb - 1] = v[ijlr1 - 1] * RS
            elif ij == ijlr2:
                x[i - 1, jj - nocb - 1] = -v[ijlr1 - 1] * RS
            else:
                x[i - 1, jj - nocb - 1] = v[ij - 1]
    return x

genmask = np.ones((noca, nvirb))
genmask[nocb:noca, 0:2] = 0.0
Xt0 = [unfold_vec(Xf[:, s]) for s in range(nstate)]
ctx = dict(nstate=nstate, noca=noca, nocb=nocb, nbf=nbf, nvirb=nvirb,
           nij=nij, noc=noc, RS=RS, genmask=genmask, Xt=Xt0)
sij0 = FK.s_ij_of(ctx, np.eye(nbf))
sab0 = FK.s_ab_of(ctx, np.eye(nbf))
sia0 = FK.s_ia_of(ctx, np.eye(nbf))
EPSA = 1e-6

def ampdir(J, dXt):
    Xp = [x.copy() for x in Xt0]
    Xm = [x.copy() for x in Xt0]
    Xp[J] = Xt0[J] + EPSA * dXt
    Xm[J] = Xt0[J] - EPSA * dXt
    Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
    Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
    return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

def space_of(i):
    return 0 if i < nocb else (1 if i < noca else 2)

natom = ncoord // 3
for I in range(nstate):
    for J in range(I + 1, nstate):
        Mt = d[f'MT_{I}{J}'] + d[f'MTG_{I}{J}']
        y = d[f'ytil_{I}{J}']
        T1 = d[f'T1_{I}{J}']
        ex = np.array([float(np.dot(y, w1[c, J])) for c in range(ncoord)]) - T1
        pr = np.array([float(np.sum(Mt * U1[c])) for c in range(ncoord)])
        D = ex - pr
        tsum = D.reshape(natom, 3).sum(axis=0)
        # candidate basis
        X = Xt0[J]
        cands = {}
        cands['gSk'] = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                                 for c in range(ncoord)])
        cands['gSk_T'] = np.array([float(np.sum(gam[J, I] * Sk_an[c]))
                                   for c in range(ncoord)])
        cands['dnum'] = dcv[I, J].reshape(-1)
        tA = np.zeros(ncoord); tB = np.zeros(ncoord)
        sA = np.zeros(ncoord); sB = np.zeros(ncoord)
        for c in range(ncoord):
            Skm = Sk_an[c]; Sxm = Sx[c]
            tA[c] = ampdir(J, Skm[:noca, :noca].T @ X)[I]
            tB[c] = ampdir(J, X @ Skm[nocb:, nocb:])[I]
            sA[c] = ampdir(J, Sxm[:noca, :noca].T @ X)[I]
            sB[c] = ampdir(J, X @ Sxm[nocb:, nocb:])[I]
        cands['trA'] = tA; cands['trB'] = tB
        cands['sxA'] = sA; cands['sxB'] = sB
        print(f'\n===== pair ({I+1},{J+1}) |D|={np.linalg.norm(D):.5f} '
              f'tsum={np.round(tsum,5)} =====')
        for k, v in cands.items():
            nv = np.linalg.norm(v)
            if nv < 1e-12:
                print(f'  {k}: ZERO')
                continue
            r = float(np.dot(D, v)) / (nv * nv)
            rr = np.linalg.norm(D - r * v) / (np.linalg.norm(D) + 1e-300)
            cc = float(np.dot(D, v)) / (np.linalg.norm(D) * nv + 1e-300)
            print(f'  {k}: cos={cc:+.4f} r={r:+.5f} rel-res={rr:.4f}')
        keys = [k for k in cands if np.linalg.norm(cands[k]) > 1e-12]
        A = np.stack([cands[k] for k in keys], axis=1)
        coef, *_ = np.linalg.lstsq(A, D, rcond=None)
        fit = A @ coef
        print(f'  LSQ[{keys}]: coef={np.round(coef, 4)} '
              f'rel-res={np.linalg.norm(D - fit) / np.linalg.norm(D):.4f}')
