"""v9: replace the closed-form Mt_G by the MEASURED vector G-channel:
  Gch[c] = ytil . [A(C, F + eps*G[dD(Ux_c)]) - A(C, F)] X_J / eps
J1' = staged:Ux + Gch  vs exact (ytil.w_ref - T1);  then the corrected
full assembly vs d_num. Uses the v7o npz (MT staged matrices + Ux + w1
+ ytil + T1) from ITS OWN process? NO -- phases! Everything in-process.
Run: python v9_gvec.py <inp> <dnum.npz>
"""
import os, sys
import numpy as np

def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    from scipy.sparse.linalg import minres, LinearOperator
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp, dnum_npz = sys.argv[1], sys.argv[2]
    dn = np.load(dnum_npz)
    dcv = dn['dcv' if 'dcv' in dn.files else dn.files[0]]
    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v9.log'))
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
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
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
            flatF[(I + J * nstate) * nbf * nbf:
                  (I + J * nstate + 1) * nbf * nbf] = gam[I, J].T.reshape(-1)
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
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

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

    print('G_met + MINRES ytil...', flush=True)
    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir(J, unfold_vec(e))
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]
    phys = [k for k in range(nij) if k != ijlr2 - 1]
    ytil = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            XJ = Xf[:, J]

            def op(v):
                vv = np.zeros(nij)
                vv[phys] = v
                vv -= XJ * float(np.dot(XJ, vv))
                av = Om[J] * vv - matvec(vv)
                av -= XJ * float(np.dot(XJ, av))
                return av[phys]

            rhs = G_met[I, J].copy()
            rhs -= XJ * float(np.dot(XJ, rhs))
            y, _ = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                          rhs[phys], rtol=1e-9, maxiter=3000)
            out = np.zeros(nij)
            out[phys] = y
            out -= XJ * float(np.dot(XJ, out))
            ytil[(I, J)] = out
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # T1 engines
    T1 = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            rr = X0_raw.copy().reshape(-1)
            rr[I * nij:(I + 1) * nij] = ytil[(I, J)]
            mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
            T1[(I, J)] = a[J, I].reshape(-1) + es
            mol.data['OQP::td_bvec_mo'] = X0_raw
    print('T1 done.', flush=True)

    # displaced sweep: w_ref (sign-fixed C) + Ux + w_skel(1-iter F[D_ref])
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = 1e-3

    def disp(coord):
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
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Wd = sg[:, None] * Wd
        mol.data['OQP::VEC_MO_A'] = Wd
        mol.data['OQP::VEC_MO_B'] = Wd.copy()
        Ax = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax[s] = matvec(Xf[:, s])
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf, Ax

    print('displaced sweep...', flush=True)
    Ux = np.zeros((ncoord, nbf, nbf))
    w_ref = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += HD
        Mp, Skp, Axp = disp(cp)
        cm = xyz0.copy(); cm[c] -= HD
        Mm, Skm, Axm = disp(cm)
        Ux[c] = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
        w_ref[c] = (Axp - Axm) / (2 * HD)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    # staged C-channel per pair (directional along Ux_c) + measured G-channel
    nbf_tri = nbf * (nbf + 1) // 2

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]; M[q, p] = pk[idx]; idx += 1
        return M

    def pack_sym(M):
        pk = np.zeros(nbf_tri)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                pk[idx] = M[p, q]; idx += 1
        return pk

    D0a = unpack_sym(SAVE0['OQP::DM_A'].ravel())
    D0b = unpack_sym(SAVE0['OQP::DM_B'].ravel())
    F0a = unpack_sym(SAVE0['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(SAVE0['OQP::FOCK_B'].ravel())
    for attr in ('maxit', 'scf_maxit'):
        try:
            setattr(mol.data._data.control, attr, 1)
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 1
    occ_a = np.zeros(nbf); occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf); occ_b[:nocb] = 1.0
    EPSG = 1e-4
    TROT = 1e-5

    print('per-coordinate staged + measured-G channels...', flush=True)
    UCH = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            UCH[(I, J)] = np.zeros(ncoord)
    base = {}
    for J in range(nstate):
        base[J] = matvec(Xf[:, J])
    for c in range(ncoord):
        # G[dD_c] build once per coordinate
        Ma = Ux[c] * occ_a[None, :] + Ux[c].T * occ_a[:, None]
        Mb = Ux[c] * occ_b[None, :] + Ux[c].T * occ_b[:, None]
        dDa = C @ Ma @ C.T
        dDb = C @ Mb @ C.T
        mol.data['OQP::DM_A'] = pack_sym(D0a + EPSG * dDa).reshape(SAVE0['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = pack_sym(D0b + EPSG * dDb).reshape(SAVE0['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'], copy=True).ravel())
        Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'], copy=True).ravel())
        for k, v in SAVE0.items():
            mol.data[k] = v.copy()
        Ga = (Fa - F0a) / EPSG
        Gb = (Fb - F0b) / EPSG
        # push FOCK + epsF*G, matvec per state (frozen C)
        epsF = 1e-4
        mol.data['OQP::FOCK_A'] = pack_sym(F0a + epsF * Ga).reshape(SAVE0['OQP::FOCK_A'].shape)
        mol.data['OQP::FOCK_B'] = pack_sym(F0b + epsF * Gb).reshape(SAVE0['OQP::FOCK_B'].shape)
        gv = {}
        for J in range(nstate):
            gv[J] = (matvec(Xf[:, J]) - base[J]) / epsF
        for k, v in SAVE0.items():
            mol.data[k] = v.copy()
        # staged C-channel
        st = {}
        for J in range(nstate):
            acc = np.zeros(nij)
            for sgn in (1.0, -1.0):
                Wp = (np.eye(nbf) + sgn * TROT * Ux[c].T) @ W0
                mol.data['OQP::VEC_MO_A'] = Wp
                mol.data['OQP::VEC_MO_B'] = Wp
                acc += sgn * matvec(Xf[:, J]) / (2 * TROT)
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            st[J] = acc
        for I in range(nstate):
            for J in range(nstate):
                if I == J:
                    continue
                UCH[(I, J)][c] = float(np.dot(ytil[(I, J)], st[J] + gv[J]))
        if c == 0:
            print('  c=0 done', flush=True)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    print('\n===== v9 J1p: [staged + measured-G] vs exact =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            ex = np.array([float(np.dot(ytil[(I, J)], w_ref[c, J]))
                           for c in range(ncoord)]) - T1[(I, J)]
            d = UCH[(I, J)] - ex
            print(f'({I+1},{J+1}): |Uch|={np.linalg.norm(UCH[(I,J)]):.5f} '
                  f'|exact|={np.linalg.norm(ex):.5f} maxdiff={np.abs(d).max():.3e}')

    print('\n===== v9 FULL ASSEMBLY vs d_num =====')
    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            gU = np.array([float(np.sum(gam[I, J] * Ux[c])) for c in range(ncoord)])
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c])) for c in range(ncoord)])
            dp[:, I, J] = T1[(I, J)] + UCH[(I, J)] + gU + gsk
    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    for I in range(nstate):
        for J in range(I + 1, nstate):
            v = dpa[:, I, J]
            ref = dcv[I, J].reshape(-1)
            cc = float(np.dot(ref, v)) / (np.linalg.norm(ref) * np.linalg.norm(v) + 1e-300)
            md = min(np.abs(v - ref).max(), np.abs(v + ref).max())
            print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(ref):.6f} '
                  f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                  f'sign-res maxdiff={md:.3e}')

if __name__ == '__main__':
    main()
