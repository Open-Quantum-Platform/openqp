"""v6b (H2O): CHANNEL-SPLIT polarization -- locate the spurious content
of the polarized production chain, pair (1,2) [broken] vs (2,3) [good].

Per slot v in {0, ytil, X, ytil+X}: solve z once, then evaluate the
gradient three ways (records re-pushed before each):
   g_full   : as-is
   g_noW    : WAO zeroed          -> W-channel = g_full - g_noW
   g_skel   : WAO and td_p zeroed -> z/P-channel = g_noW - g_skel
Polarize each channel; engines (slot-injection amp+esum) referee the
skeleton channel; ytil.w_ref referees the exact total.

Run: python v6b_h2o.py <input.inp>
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
    from scipy.sparse.linalg import minres, LinearOperator
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v6b.log'))
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

    PAIRS = [(0, 1), (1, 2)]     # (1,2) broken, (2,3) good
    print('G_met sweeps (targets of the two pairs)...', flush=True)
    G_need = {}
    for I, J in PAIRS:
        gv = np.zeros(nij)
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            gv[k] = ampdir_unf(J, unfold_vec(e))[I]
        G_need[(I, J)] = gv

    def ytil_minres(I, J):
        XJ = Xf[:, J]
        phys = [k for k in range(nij) if k != ijlr2 - 1]

        def op(v):
            vv = np.zeros(nij)
            vv[phys] = v
            vv -= XJ * float(np.dot(XJ, vv))
            av = Om[J] * vv - matvec(vv)
            av -= XJ * float(np.dot(XJ, av))
            return av[phys]

        rhs = G_need[(I, J)].copy()
        rhs -= XJ * float(np.dot(XJ, rhs))
        y, info = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                         rhs[phys], rtol=1e-9, maxiter=3000)
        out = np.zeros(nij)
        out[phys] = y
        out -= XJ * float(np.dot(XJ, out))
        print(f'  MINRES ({I+1},{J+1}): info={info} '
              f'resid={np.linalg.norm(op(y)-rhs[phys]):.2e}', flush=True)
        return out

    ytil = {}
    for I, J in PAIRS:
        ytil[(I, J)] = ytil_minres(I, J)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- exact referee: w_ref sweep --------------------------------------
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
            Ax_ref[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        return Ax_ref

    print('w_ref sweep...', flush=True)
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

    # ---- channel-split chain ---------------------------------------------
    def chain3(J, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[J * nij:(J + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(J + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        tdp = np.array(mol.data['OQP::td_p'], copy=True)
        wao = np.array(mol.data['OQP::WAO'], copy=True)
        oqp.tdhf_mrsf_gradient(mol)
        g_full = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::WAO'] = np.zeros_like(wao)
        mol.data['OQP::td_p'] = tdp.copy()
        oqp.tdhf_mrsf_gradient(mol)
        g_noW = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::WAO'] = np.zeros_like(wao)
        mol.data['OQP::td_p'] = np.zeros_like(tdp)
        oqp.tdhf_mrsf_gradient(mol)
        g_skel = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::WAO'] = wao
        mol.data['OQP::td_p'] = tdp
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return g_full, g_noW, g_skel

    # slot-injection engines for the skeleton referee
    def engines(I, J, yv):
        rr = X0_raw.copy().reshape(-1)
        rr[I * nij:(I + 1) * nij] = yv
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_nac_amp(mol)
        a = np.array(mol.data['OQP::nac_amp'], copy=True
                     ).reshape((nstate, nstate, natom, 3))
        oqp.mrsf_nac_esum(mol, I + 1, J + 1)
        es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return a[J, I].reshape(-1) + es

    print('channel-split polarization...', flush=True)
    z0 = chain3(0, np.zeros(nij))
    for I, J in PAIRS:
        y = ytil[(I, J)]
        cm = chain3(J, y + Xf[:, J])
        cy = chain3(J, y)
        cx = chain3(J, Xf[:, J])
        pol = [0.5 * (cm[k] - cy[k] - cx[k] + z0[k]) for k in range(3)]
        g_full_p, g_noW_p, g_skel_p = pol
        Wch = g_full_p - g_noW_p
        Pch = g_noW_p - g_skel_p
        eng = engines(I, J, y)
        yw = np.array([float(np.dot(y, w_ref[c, J])) for c in range(ncoord)])
        print(f'\n===== pair ({I+1},{J+1}) channel split =====')
        print(f'  exact ytil.w_ref : {np.round(yw, 6)}')
        print(f'  g_full_pol       : {np.round(g_full_p, 6)}')
        print(f'  skel_pol         : {np.round(g_skel_p, 6)}')
        print(f'  engines(referee) : {np.round(eng, 6)}')
        print(f'  z/P-channel_pol  : {np.round(Pch, 6)}')
        print(f'  W-channel_pol    : {np.round(Wch, 6)}')
        print(f'  |skel-eng|={np.abs(g_skel_p-eng).max():.3e} '
              f'|full-exact|={np.abs(g_full_p-yw).max():.3e}', flush=True)


if __name__ == '__main__':
    main()
