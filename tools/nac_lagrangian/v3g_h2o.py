"""v3g (H2O): the exact-w referee. Hypothesis: my w = w_skel + w_rot
misses the Fock DENSITY-RESPONSE channel G[dD] (mrsf_matvec_apply is
frozen-Fock: it never rebuilds F_AO from the rotated density).

 (W1) w_ref = FD of the FULL displaced matvec A(x')X_J (displaced SCF,
      sg gauge round-trip) -- the exact moving-frame dA X_J.
 (W2) w_ref vs w_mine: norms/cos + eigencomponent structure of the diff.
 (W3) sum rule X_J^T w vs FD d(omega_J).
 (W4) PT(w_ref) vs dX_FD  -- certifies the PT formula itself.
 (W5) final: antisym[ampdir(PT(w_ref)) + gam:(Sk_an + Ux_FD)] vs d_num.

Run: python v3g_h2o.py <input.inp> <dnum.npz> <skeldir>
"""
import os
import sys
import numpy as np

TROT = 1e-5
EPSA = 1e-5
H = 1e-3


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, dnum_npz, skeldir = sys.argv[1], sys.argv[2], sys.argv[3]
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v3g.log'))
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

    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False

    def displaced(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        om_d = [mol.energies[k + 1] - mol.energies[0] for k in range(nstate)]
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
        M_f = M_f * sg[None, :]
        Xm_ = np.array(mol.data['OQP::td_bvec_mo'], copy=True
                       ).reshape(-1).reshape((nstate, nij)).T.copy()
        Xd = [unfold_vec(Xm_[:, s]) for s in range(nstate)]
        sgo = sg[:noca]
        sgv = sg[nocb:]
        Xd = [sgo[:, None] * x * sgv[None, :] for x in Xd]
        for s in range(nstate):
            if np.sum(Xt0[s] * Xd[s]) < 0:
                Xd[s] = -Xd[s]
        # exact w referee: full displaced matvec on the sg-transported
        # reference amplitudes, result transported back
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            vin = sg_apply(Xf[:, s], sg)
            Ax = matvec(vin)
            Ax_ref[s] = sg_apply(Ax, sg)
        # Sk with frozen C
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Xd, Skf, Ax_ref, om_d

    print('displaced sweep (full matvec referee)...', flush=True)
    T_FD = np.zeros((ncoord, nbf, nbf))
    Sk_FD = np.zeros((ncoord, nbf, nbf))
    dX_FD = np.zeros((ncoord, nstate, noca, nvirb))
    w_ref = np.zeros((ncoord, nstate, nij))
    dom_FD = np.zeros((ncoord, nstate))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += H
        Mp, Xp, Skp, Axp, omp = displaced(cp)
        cm = xyz0.copy(); cm[c] -= H
        Mm, Xm2, Skm, Axm, omm = displaced(cm)
        T_FD[c] = (Mp - Mm) / (2 * H)
        Sk_FD[c] = (Skp - Skm) / (2 * H)
        w_ref[c] = (Axp - Axm) / (2 * H)
        dom_FD[c] = (np.array(omp) - np.array(omm)) / (2 * H)
        for s in range(nstate):
            dX_FD[c, s] = (Xp[s] - Xm2[s]) / (2 * H)
    Ux_FD = T_FD - Sk_FD

    # restore, rebuild my w
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw
    w_skel = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        fp = np.load(os.path.join(skeldir, f'p{c}.npz'))
        fm = np.load(os.path.join(skeldir, f'p{c + ncoord}.npz'))
        w_skel[c] = (fp['Ax'] - fm['Ax']) / (2 * H)
    print('rotation FD sweep...', flush=True)
    w_rot = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        for sgn in (1.0, -1.0):
            Wp = (np.eye(nbf) + sgn * TROT * Ux_FD[c].T) @ W0
            mol.data['OQP::VEC_MO_A'] = Wp
            mol.data['OQP::VEC_MO_B'] = Wp
            for s in range(nstate):
                w_rot[c, s] += sgn * matvec(Xf[:, s]) / (2 * TROT)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw
    w_mine = w_skel + w_rot

    print('\n===== W2/W3: w_ref vs w_mine, sum rule =====')
    for J in range(nstate):
        for c in (0, 2, 5):
            a, b = w_ref[c, J], w_mine[c, J]
            cc = float(np.dot(a, b)) / (np.linalg.norm(a)
                                        * np.linalg.norm(b) + 1e-300)
            d = a - b
            d[ijlr2 - 1] = 0.0
            cf = V.T @ d
            top = np.argsort(-np.abs(cf))[:3]
            print(f'J={J+1} c={c}: |w_ref|={np.linalg.norm(a):.4f} '
                  f'|w_mine|={np.linalg.norm(b):.4f} cos={cc:+.5f} '
                  f'|diff|={np.linalg.norm(d):.4f} top: '
                  + ', '.join(f'k={k}(w={w_full[k]:.3f},dc={cf[k]:+.3f})'
                              for k in top))
        sr_ref = np.array([float(np.dot(Xf[:, J], w_ref[c, J]))
                           for c in range(ncoord)])
        sr_mine = np.array([float(np.dot(Xf[:, J], w_mine[c, J]))
                            for c in range(ncoord)])
        print(f'J={J+1} sumrule: max|X^Tw_ref - dom|={np.abs(sr_ref-dom_FD[:,J]).max():.3e} '
              f'max|X^Tw_mine - dom|={np.abs(sr_mine-dom_FD[:,J]).max():.3e}')
    print('', flush=True)

    def PT_of(wsrc):
        dX = np.zeros((ncoord, nstate, nij))
        for J in range(nstate):
            for c in range(ncoord):
                w = wsrc[c, J].copy()
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
                dX[c, J] = acc
        return dX

    dX_PTref = PT_of(w_ref)
    print('===== W4: PT(w_ref) vs dX_FD =====')
    for J in range(nstate):
        wa = 0.0
        for c in range(ncoord):
            fd = dX_FD[c, J]
            pt = unfold_vec(dX_PTref[c, J])
            if np.linalg.norm(fd) < 1e-6:
                continue
            cc = float(np.sum(fd * pt)) / (np.linalg.norm(fd)
                                           * np.linalg.norm(pt) + 1e-300)
            wa = max(wa, np.abs(pt - fd).max())
            if c in (0, 5):
                print(f'J={J+1} c={c}: |FD|={np.linalg.norm(fd):.6f} '
                      f'|PT|={np.linalg.norm(pt):.6f} cos={cc:+.6f}')
        print(f'J={J+1}: max|PT-FD| = {wa:.5f}')
    print('', flush=True)

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

    def judge(tag, dxsrc, Tsrc):
        damp = np.zeros((ncoord, nstate, nstate))
        dorb = np.zeros((ncoord, nstate, nstate))
        for J in range(nstate):
            for c in range(ncoord):
                g = ampdir_unf(J, dxsrc(c, J))
                for I in range(nstate):
                    if I != J:
                        damp[c, I, J] = g[I]
        for I in range(nstate):
            for J in range(nstate):
                if I == J:
                    continue
                for c in range(ncoord):
                    dorb[c, I, J] = float(np.sum(gam[I, J] * Tsrc[c]))
        dp = damp + dorb
        dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
        print(f'===== {tag} =====')
        for I in range(nstate):
            for J in range(nstate):
                if I >= J:
                    continue
                dn = dcv_n[I, J].reshape(-1)
                v = dpa[:, I, J]
                cc = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                             * np.linalg.norm(v) + 1e-300)
                print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(dn):.6f} '
                      f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                      f'maxdiff={np.abs(v-dn).max():.3e}')
        print('', flush=True)

    judge('W5: antisym[ampdir(PT(w_ref)) + gam:(Sk_an+Ux_FD)]',
          lambda c, J: unfold_vec(dX_PTref[c, J]), Sk_an + Ux_FD)
    np.savez(inp.replace('.inp', '_v3g.npz'), w_ref=w_ref, w_mine=w_mine,
             dX_FD=dX_FD, dX_PTref=dX_PTref, Ux_FD=Ux_FD, Sk_FD=Sk_FD,
             T_FD=T_FD, dom_FD=dom_FD, w_full=w_full)


if __name__ == '__main__':
    main()
