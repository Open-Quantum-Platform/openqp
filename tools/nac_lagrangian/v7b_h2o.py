"""v7b: v7a + the same-space antisym channel T4 (FD-U referee).

  d[I,J]^c = [amp2e+esum](ytil,X_J)^c            (T1, slot injection)
           + zB[Lt_a(ytil,X_J) + gamma_a]^c      (T2, seam)
           + T3^c                                 (FULL sym channel)
           + gamma:Sk_an^c
  T3^c = staged directional derivative of E(ytil,X_J;C) along
         V_c = -1/2 S^x_MO (frozen-Fock matvec, 2 calls/coord)
       + sum_sigma Tr[dD_sigma(V_c) . G_sigma[P~]]   (Fock density
         response; G built in-process: hf_energy at maxit=1 on D+eps*P~)

Known omission: same-space antisym canonical channel (T4, v7b).
Run: python v7b_h2o.py <input.inp> <dnum.npz>
"""
import os
import sys
import numpy as np

EPSA = 1e-5
TROT = 1e-5
EPSG = 1e-4


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    from scipy.sparse.linalg import minres, LinearOperator
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, dnum_npz = sys.argv[1], sys.argv[2]
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    os.environ['NAC_DUMP_PIJ'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7b.log'))
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
    e0a_r = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b_r = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca
    nbf_tri = nbf * (nbf + 1) // 2

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

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]
                M[q, p] = pk[idx]
                idx += 1
        return M

    def pack_sym(M):
        pk = np.zeros(nbf_tri)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                pk[idx] = M[p, q]
                idx += 1
        return pk

    # gamma + dS exports
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
    dsf_raw = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    Sx_MO = np.zeros((ncoord, nbf, nbf))
    dsf_AO = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        dsk = dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        dsf_AO[c] = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        Sk_an[c] = C.T @ dsk @ C
        Sx_MO[c] = C.T @ dsf_AO[c] @ C

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

    print('G_met sweeps...', flush=True)
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

    phys = [k for k in range(nij) if k != ijlr2 - 1]

    def ytil_minres(I, J):
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
        y, info = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                         rhs[phys], rtol=1e-9, maxiter=3000)
        out = np.zeros(nij)
        out[phys] = y
        out -= XJ * float(np.dot(XJ, out))
        return out

    print('MINRES ytil solves...', flush=True)
    ytil = {}
    for I in range(nstate):
        for J in range(nstate):
            if I != J:
                ytil[(I, J)] = ytil_minres(I, J)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T1 engines + pij densities --------------------------------------
    def inject(I, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[I * nij:(I + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)

    T1 = {}
    PIJ = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            inject(I, ytil[(I, J)])
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
            pa = np.array(mol.data['OQP::dbg_pij_a'], copy=True
                          ).reshape(-1).reshape(nbf, nbf).T
            pb = np.array(mol.data['OQP::dbg_pij_b'], copy=True
                          ).reshape(-1).reshape(nbf, nbf).T
            T1[(I, J)] = a[J, I].reshape(-1) + es
            PIJ[(I, J)] = (pa.copy(), pb.copy())
            mol.data['OQP::td_bvec_mo'] = X0_raw
    print('T1 engines done.', flush=True)

    # ---- T3 staged part ---------------------------------------------------
    print('T3 staged sym sweeps...', flush=True)
    T3s = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[(I, J)]
            t3 = np.zeros(ncoord)
            for c in range(ncoord):
                V = -0.5 * Sx_MO[c]
                acc = 0.0
                for sgn in (1.0, -1.0):
                    Wp = (np.eye(nbf) + sgn * TROT * V.T) @ W0
                    mol.data['OQP::VEC_MO_A'] = Wp
                    mol.data['OQP::VEC_MO_B'] = Wp
                    Ax = matvec(Xf[:, J])
                    acc += sgn * float(np.dot(y, Ax)) / (2 * TROT)
                t3[c] = acc
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            T3s[(I, J)] = t3
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T3 G-channel -----------------------------------------------------
    SAVE_KEYS = ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
                 'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
                 'OQP::E_MO_B']
    saved = {k: np.array(mol.data[k], copy=True) for k in SAVE_KEYS}
    F0a = unpack_sym(saved['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(saved['OQP::FOCK_B'].ravel())
    D0a = unpack_sym(saved['OQP::DM_A'].ravel())
    D0b = unpack_sym(saved['OQP::DM_B'].ravel())
    ok_maxit = False
    for attr in ('maxit', 'scf_maxit', 'maxit_scf'):
        try:
            setattr(mol.data._data.control, attr, 1)
            ok_maxit = True
            print(f'set control.{attr} = 1')
            break
        except Exception:
            continue
    mol.config['scf']['maxit'] = 1
    print(f'maxit control set: {ok_maxit}', flush=True)

    def gbuild(Pa, Pb):
        mol.data['OQP::DM_A'] = pack_sym(D0a + EPSG * Pa).reshape(
            saved['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = pack_sym(D0b + EPSG * Pb).reshape(
            saved['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'],
                                 copy=True).ravel())
        Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'],
                                 copy=True).ravel())
        for k in SAVE_KEYS:
            mol.data[k] = saved[k].copy()
        return (Fa - F0a) / EPSG, (Fb - F0b) / EPSG

    print('T3 G-channel builds...', flush=True)
    occ_a = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf)
    occ_b[:nocb] = 1.0
    T3g = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Pa, Pb = PIJ[(I, J)]
            Ga, Gb = gbuild(Pa, Pb)
            t3 = np.zeros(ncoord)
            for c in range(ncoord):
                V = -0.5 * Sx_MO[c]
                Ma = V * (occ_a[:, None] + occ_a[None, :])
                Mb = V * (occ_b[:, None] + occ_b[None, :])
                dDa = C @ Ma @ C.T
                dDb = C @ Mb @ C.T
                t3[c] = float(np.sum(Ga * dDa) + np.sum(Gb * dDb))
            T3g[(I, J)] = t3
            print(f'  ({I+1},{J+1}) |G-ch|={np.linalg.norm(t3):.6f} '
                  f'|staged|={np.linalg.norm(T3s[(I,J)]):.6f}', flush=True)

    # ---- displaced sweep for Ux_FD (T4 referee U) -------------------------
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = 1.0e-3

    def displaced_M(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a_r
        mol.data['OQP::E_MO_B_old'] = e0b_r
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf

    print('displaced Ux sweep...', flush=True)
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += HD
        Mp, Skp = displaced_M(cp)
        cm = xyz0.copy(); cm[c] -= HD
        Mm, Skm = displaced_M(cm)
        Ux_FD[c] = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    # ---- same-space antisym Mt sweep + T4 ---------------------------------
    ss_pairs = []
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        for p in range(lo, hi):
            for q in range(lo, p):
                ss_pairs.append((p, q))
    print(f'same-space antisym sweep: {len(ss_pairs)} generators x pairs...',
          flush=True)
    T4 = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[(I, J)]
            S_swp = np.zeros(len(ss_pairs))
            for ip, (p, q) in enumerate(ss_pairs):
                K = np.zeros((nbf, nbf))
                K[p, q] = 1.0
                K[q, p] = -1.0
                acc = 0.0
                for sgn in (1.0, -1.0):
                    Wp = (np.eye(nbf) - sgn * TROT * K) @ W0
                    mol.data['OQP::VEC_MO_A'] = Wp
                    mol.data['OQP::VEC_MO_B'] = Wp
                    Ax = matvec(Xf[:, J])
                    acc += sgn * float(np.dot(y, Ax)) / (2 * TROT)
                S_swp[ip] = acc
            mol.data['OQP::VEC_MO_A'] = W0
            mol.data['OQP::VEC_MO_B'] = Wb0
            t4 = np.zeros(ncoord)
            for c in range(ncoord):
                Ua = 0.5 * (Ux_FD[c] - Ux_FD[c].T)
                acc = 0.0
                for ip, (p, q) in enumerate(ss_pairs):
                    acc += (S_swp[ip] + 2.0 * gam[I, J][p, q]) * Ua[p, q]
                t4[c] = acc
            T4[(I, J)] = t4
            print(f'  ({I+1},{J+1}) |T4|={np.linalg.norm(t4):.6f}', flush=True)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- T2 seam ----------------------------------------------------------
    os.environ['NAC_DUMP_RHS'] = '1'

    def rhs_of(target, amp_vec):
        rr = X0_raw.copy().reshape(-1)
        rr[(target - 1) * nij:target * nij] = amp_vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(target)
        oqp.tdhf_mrsf_z_vector(mol)
        out = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return out

    def unpack_rot_to_mat(v):
        Lm = np.zeros((nbf, nbf))
        idx = 0
        for ii in range(nocb, noca):
            for jj in range(nocb):
                Lm[ii, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for jj in range(nocb):
                Lm[kk, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for ii in range(nocb, noca):
                Lm[kk, ii] = v[idx]; idx += 1
        return Lm

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    print('T2 seams...', flush=True)
    T2 = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[(I, J)]
            qy = rhs_of(J + 1, y)
            qx = rhs_of(J + 1, Xf[:, J])
            qm = rhs_of(J + 1, y + Xf[:, J])
            Lpol = 0.5 * (qm - qy - qx)
            Lmat = unpack_rot_to_mat(Lpol) + gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            if not bool(mol.mol_energy.Z_Vector_converged):
                oqp.set_mrsf_nac_cphf(mol, 0, 0)
                T2[(I, J)] = None
                print(f'  ({I+1},{J+1}): z NOT converged', flush=True)
                continue
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2[(I, J)] = gZ - gS
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass

    # ---- assemble ----------------------------------------------------------
    print('\n===== v7b ASSEMBLY (+same-space T4) vs d_num =====')
    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J or T2[(I, J)] is None:
                continue
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            dp[:, I, J] = (T1[(I, J)] + T2[(I, J)] + T3s[(I, J)]
                           + T3g[(I, J)] + T4[(I, J)] + gsk)
    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            dn = dcv_n[I, J].reshape(-1)
            v = dpa[:, I, J]
            cc = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                         * np.linalg.norm(v) + 1e-300)
            md = min(np.abs(v - dn).max(), np.abs(v + dn).max())
            print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(dn):.6f} '
                  f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                  f'sign-resolved maxdiff={md:.3e}')
    np.savez(inp.replace('.inp', '_v7b.npz'), dp=dp, dpa=dpa)


if __name__ == '__main__':
    main()
