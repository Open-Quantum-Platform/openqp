"""v4a ADJOINT GATE: verify the production per-pair adjoint form of the
certified 7.27 assembly, using the saved exact w_ref data.

  G_met[I,J] : the formula-metric gradient vector  (ampdir is linear:
               ampdir_J(v)[I] = G_met[I,J] . v), extracted by unit sweep
  ytil[I,J]  = sum_{k != kJ} V_k (V_k^T G_met[I,J]) / (om_J - w_k)
  ADJOINT IDENTITY (must be machine-exact, gauge-free):
     ytil[I,J] . w^c  ==  G_met[I,J] . PT(w^c)      for every c
  Then d^c = antisym[ ytil.w_ref^c + gam:(Sk_an + Ux)^c ] vs d_num
  (per-pair phase ambiguity vs the saved w_ref is reported via the
   pair-product rule).

Run: python v4a_adjoint.py <input.inp> <v3g.npz> <dnum.npz>
"""
import os
import sys
import numpy as np

EPSA = 1e-5


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, g_npz, dnum_npz = sys.argv[1], sys.argv[2], sys.argv[3]
    G = np.load(g_npz)
    w_ref, Ux, Sk_FD = G['w_ref'], G['Ux_FD'], G['Sk_FD']
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v4a.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E0 = list(mol.energies)
    Om = [E0[k + 1] - E0[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb, nvirb, RS = ctx['noca'], ctx['nocb'], ctx['nvirb'], ctx['RS']
    X0_raw = ctx['X0_raw']
    Xt0 = ctx['Xt']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    C = W0.T
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
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

    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    nbf2 = nbf * nbf
    dsk_raw = np.array(mol.data['OQP::dbg_dsket'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        Sk_an[c] = C.T @ dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T @ C

    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'],
                        copy=True).ravel()[:nij].copy()

    print(f'building full A ({nij}x{nij})...', flush=True)
    A = np.zeros((nij, nij))
    for k in range(nij):
        e = np.zeros(nij)
        e[k] = 1.0
        A[:, k] = matvec(e)
    mol.data['OQP::td_bvec_mo'] = X0_raw
    phys = [k for k in range(nij) if k != ijlr2 - 1]
    Ap = 0.5 * (A[np.ix_(phys, phys)] + A[np.ix_(phys, phys)].T)
    w_full, Vp = np.linalg.eigh(Ap)
    V = np.zeros((nij, nij - 1))
    V[phys, :] = Vp
    ks = [int(np.argmax(np.abs(V.T @ Xf[:, s]))) for s in range(nstate)]

    # ---- G_met extraction (unit sweep over the raw amplitude space) ------
    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))

    def ampdir_vec(J, v):
        dXt = unfold_vec(v)
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    print('extracting G_met (unit sweep)...', flush=True)
    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir_vec(J, e)
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]
    # linearity check on a random direction
    rng_v = np.zeros(nij)
    rng_v[phys] = np.sin(np.arange(nij - 1) + 1.0)
    lin_err = 0.0
    for J in range(nstate):
        gd = ampdir_vec(J, rng_v)
        for I in range(nstate):
            if I != J:
                lin_err = max(lin_err, abs(gd[I] - float(np.dot(
                    G_met[I, J], rng_v))))
    print(f'G_met linearity check: max err = {lin_err:.3e}', flush=True)

    # ---- adjoint solves + identity ---------------------------------------
    ytil = np.zeros((nstate, nstate, nij))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            coef = V.T @ G_met[I, J]
            acc = np.zeros(nij)
            for k in range(nij - 1):
                if k == ks[J]:
                    continue
                den = Om[J] - w_full[k]
                if abs(den) < 1e-9:
                    continue
                acc += V[:, k] * (coef[k] / den)
            ytil[I, J] = acc

    def PT_of_w(J, w):
        w = w.copy()
        w[ijlr2 - 1] = 0.0
        coef = V.T @ w
        acc = np.zeros(nij)
        for k in range(nij - 1):
            if k == ks[J]:
                continue
            den = Om[J] - w_full[k]
            if abs(den) < 1e-9:
                continue
            acc += V[:, k] * (coef[k] / den)
        return acc

    print('\n===== ADJOINT IDENTITY (must be machine-exact) =====')
    aid = 0.0
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            for c in range(ncoord):
                lhs = float(np.dot(ytil[I, J], w_ref[c, J]))
                rhs = float(np.dot(G_met[I, J], PT_of_w(J, w_ref[c, J])))
                aid = max(aid, abs(lhs - rhs))
    print(f'max |ytil.w - G_met.PT(w)| = {aid:.3e}', flush=True)

    # ---- full d via the adjoint form -------------------------------------
    damp = np.zeros((ncoord, nstate, nstate))
    dorb = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            for c in range(ncoord):
                damp[c, I, J] = float(np.dot(ytil[I, J], w_ref[c, J]))
                dorb[c, I, J] = float(np.sum(gam[I, J] * (Sk_an[c] + Ux[c])))
    dp = damp + dorb
    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    print('\n===== ADJOINT-FORM d vs d_num (pair-phase up to product rule) =====')
    signs = []
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            dn = dcv_n[I, J].reshape(-1)
            v = dpa[:, I, J]
            cc = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                         * np.linalg.norm(v) + 1e-300)
            signs.append(np.sign(cc))
            print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(dn):.6f} '
                  f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                  f'maxdiff={min(np.abs(v-dn).max(), np.abs(v+dn).max()):.3e} (sign-resolved)')
    print(f'pair-sign product = {np.prod(signs):+.0f} (must be +1)')


if __name__ == '__main__':
    main()
