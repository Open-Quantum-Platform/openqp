"""v6 (H2O): dissect what the polarized production-gradient chain
actually computes, against the exact in-process referee.

  A = g_pol^c[I,J]      (4-term polarization of the normal-mode chain)
  B = ytil . w_ref^c    (exact moving-frame ytil^T dA X_J)
  A - B = the chain's elimination bookkeeping (+ any eigen-substitution
          garbage). Printed next to:
     wsx_e (engine Fock-W), -Tr[WAO_pol S^x] (polarized full W), and
     zB_L (seam with L-only RHS) to identify the decomposition.

Run: python v6_h2o.py <input.inp>
"""
import os
import sys
import numpy as np

EPSA = 1e-5
H = 1e-3


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v6.log'))
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
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
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

    def refold(x):
        v = np.zeros(nij)
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    v[ijlr1 - 1] = x[i - 1, jj - nocb - 1] / RS
                elif ij == ijlr2:
                    pass
                else:
                    v[ij - 1] = x[i - 1, jj - nocb - 1]
        return v

    def sg_apply(v, sg):
        x = unfold_vec(v)
        x = sg[:noca, None] * x * sg[None, nocb:]
        return refold(x)

    # dSfull export for the W.S^x contraction
    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    nbf2 = nbf * nbf
    dsf_raw = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    dsf = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        dsf[c] = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T

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

    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))

    def ampdir_unf(J, dXt):
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    print('G_met sweep...', flush=True)
    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir_unf(J, unfold_vec(e))
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]
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

    # ---- w_ref sweep -----------------------------------------------------
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False

    def displaced_wref(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            vin = sg_apply(Xf[:, s], sg)
            Ax = matvec(vin)
            Ax_ref[s] = sg_apply(Ax, sg)
        return Ax_ref

    print('displaced w_ref sweep...', flush=True)
    w_ref = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += H
        Axp = displaced_wref(cp)
        cm = xyz0.copy(); cm[c] -= H
        Axm = displaced_wref(cm)
        w_ref[c] = (Axp - Axm) / (2 * H)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- normal-mode polarized chain (g and WAO) -------------------------
    def chain(J, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[J * nij:(J + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(J + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        wao = np.array(mol.data['OQP::WAO'], copy=True).ravel().copy()
        oqp.tdhf_mrsf_gradient(mol)
        g = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return g, wao

    print('polarized chain runs...', flush=True)
    g00, wao00 = chain(0, np.zeros(nij))
    res = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            gm, wm = chain(J, ytil[I, J] + Xf[:, J])
            gy, wy = chain(J, ytil[I, J])
            gx, wx = chain(J, Xf[:, J])
            g_pol = 0.5 * (gm - gy - gx + g00)
            wao_pol = 0.5 * (wm - wy - wx + wao00)
            res[(I, J)] = (g_pol, wao_pol)
            print(f'  ({I+1},{J+1}) done', flush=True)

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]
                M[q, p] = pk[idx]
                idx += 1
        return M

    print('\n===== v6 DISSECTION: g_pol vs ytil.w_ref =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            g_pol, wao_pol = res[(I, J)]
            yw = np.array([float(np.dot(ytil[I, J], w_ref[c, J]))
                           for c in range(ncoord)])
            diff = g_pol - yw
            # candidate explanation: -Tr[W_pol S^x]
            nw = len(wao_pol)
            wsx_p = np.zeros(ncoord)
            if nw == nbf * (nbf + 1) // 2 * 2:
                Wm = unpack_sym(wao_pol[:nw // 2]) + unpack_sym(wao_pol[nw // 2:])
            elif nw == nbf * (nbf + 1) // 2:
                Wm = unpack_sym(wao_pol)
            else:
                Wm = None
                print(f'  WAO length {nw} unexpected (nbf2={nbf*(nbf+1)//2})')
            if Wm is not None:
                for c in range(ncoord):
                    wsx_p[c] = -float(np.sum(Wm * dsf[c]))
            cc = float(np.dot(diff, wsx_p)) / (np.linalg.norm(diff)
                                               * np.linalg.norm(wsx_p) + 1e-300)
            print(f'pair ({I+1},{J+1}): |g_pol|={np.linalg.norm(g_pol):.6f} '
                  f'|ytil.w_ref|={np.linalg.norm(yw):.6f} '
                  f'|diff|={np.linalg.norm(diff):.6f} '
                  f'|WpolSx|={np.linalg.norm(wsx_p):.6f} '
                  f'cos(diff,WpolSx)={cc:+.4f} '
                  f'|diff-WpolSx|={np.linalg.norm(diff-wsx_p):.6f}')
    np.savez(inp.replace('.inp', '_v6.npz'),
             w_ref=w_ref, ytil=ytil, G_met=G_met,
             g_pol=np.array([[res.get((i, j), (np.zeros(ncoord),))[0]
                              if (i, j) in res else np.zeros(ncoord)
                              for j in range(nstate)] for i in range(nstate)]))


if __name__ == '__main__':
    main()
