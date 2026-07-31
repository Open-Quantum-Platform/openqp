"""v7g: U-channel forensics with the EXACT FD referee:
  exact_U = ytil.w_ref - T1 + gamma:Ux_FD
  prod_U  = -seam(Lt+gamma) + T3pp
plus the gamma-side split and per-piece dumps.
  seam(e_pq) = -U^x_FULL_pq  =>  T2 = -seam(Lt+gamma)
  elimination V-mask = ov-engine weights: 1/2 same-space, 1 on the
  (lo,hi) cross rows; T3 = staged(V) + Tr[dD(V) G[P~]] + gamma:V.

  d[I,J]^c = [amp2e+esum](ytil,X_J)^c            (T1, slot injection)
           + zB[Lt_a(ytil,X_J) + gamma_a]^c      (T2, seam)
           + T3^c                                 (FULL sym channel)
           + gamma:Sk_an^c
  T3^c = staged directional derivative of E(ytil,X_J;C) along
         V_c = -1/2 S^x_MO (frozen-Fock matvec, 2 calls/coord)
       + sum_sigma Tr[dD_sigma(V_c) . G_sigma[P~]]   (Fock density
         response; G built in-process: hf_energy at maxit=1 on D+eps*P~)

Known omission: same-space antisym canonical channel (T4, v7b).
Run: python v7g_h2o.py <input.inp> <dnum.npz>
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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7g.log'))
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

    # elimination direction V_c with the ov-engine weights:
    # V[r,s] = -S^x[r,s] * (1/2 if same-space else 1 if space(r)<space(s) else 0)
    def space_of(i):
        return 0 if i < nocb else (1 if i < noca else 2)
    Wmask = np.zeros((nbf, nbf))
    for rr_ in range(nbf):
        for ss_ in range(nbf):
            sr, sc = space_of(rr_), space_of(ss_)
            if sr == sc:
                Wmask[rr_, ss_] = 0.5
            elif sr < sc:
                Wmask[rr_, ss_] = 1.0
    Vmask = [ -Sx_MO[c] * Wmask for c in range(ncoord) ]

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

    # ---- displaced w_ref + Ux sweep (exact referee) -----------------------
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

    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
              'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = 1.0e-3

    def displaced_all(coord):
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
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax_ref[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf, Ax_ref

    print('displaced w_ref/Ux sweep...', flush=True)
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    w_ref = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += HD
        Mp, Skp, Axp = displaced_all(cp)
        cm = xyz0.copy(); cm[c] -= HD
        Mm, Skm, Axm = displaced_all(cm)
        Ux_FD[c] = (Mp - Mm) / (2 * HD) - (Skp - Skm) / (2 * HD)
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
                V = Vmask[c]
                acc = 0.0
                for sgn in (1.0, -1.0):
                    Wp = (np.eye(nbf) + sgn * TROT * V.T) @ W0
                    mol.data['OQP::VEC_MO_A'] = Wp
                    mol.data['OQP::VEC_MO_B'] = Wp
                    Ax = matvec(Xf[:, J])
                    acc += sgn * float(np.dot(y, Ax)) / (2 * TROT)
                # gamma's own elimination content along the same V
                t3[c] = acc + float(np.sum(gam[I, J] * V))
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
                V = Vmask[c]
                Ma = V * occ_a[None, :] + V.T * occ_a[:, None]
                Mb = V * occ_b[None, :] + V.T * occ_b[:, None]
                dDa = C @ Ma @ C.T
                dDb = C @ Mb @ C.T
                t3[c] = float(np.sum(Ga * dDa) + np.sum(Gb * dDb))
            T3g[(I, J)] = t3
            print(f'  ({I+1},{J+1}) |G-ch|={np.linalg.norm(t3):.6f} '
                  f'|staged|={np.linalg.norm(T3s[(I,J)]):.6f}', flush=True)

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
    T2g_only = {}
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
            # gamma-only seam
            mol.data['OQP::nac_orbgrad_L'] = gam[I, J].ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            gZ2 = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS2 = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2g_only[(I, J)] = gZ2 - gS2
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass

    # ---- U-channel forensics ----------------------------------------------
    print('\n===== v7g U-CHANNEL FORENSICS =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J or T2[(I, J)] is None:
                continue
            y = ytil[(I, J)]
            yw = np.array([float(np.dot(y, w_ref[c, J])) for c in range(ncoord)])
            gU = np.array([float(np.sum(gam[I, J] * Ux_FD[c]))
                           for c in range(ncoord)])
            exact_U = yw - T1[(I, J)] + gU
            prod_U = -T2[(I, J)] + T3s[(I, J)] + T3g[(I, J)]
            # gamma-side split
            gV = np.array([float(np.sum(gam[I, J] * Vmask[c]))
                           for c in range(ncoord)])
            gU_prod = -T2g_only[(I, J)] + gV
            print(f'pair ({I+1},{J+1}):')
            print(f'  exact_U : {np.round(exact_U, 5)}')
            print(f'  prod_U  : {np.round(prod_U, 5)}')
            print(f'  diff    : {np.round(prod_U - exact_U, 5)}')
            print(f'  [gamma-side] exact g:U : {np.round(gU, 5)}')
            print(f'  [gamma-side] prod      : {np.round(gU_prod, 5)}')
            print(f'  pieces: -seam_full={np.round(-T2[(I,J)],5)}')
            print(f'          -seam_g   ={np.round(-T2g_only[(I,J)],5)}')
            print(f'          T3staged+gV={np.round(T3s[(I,J)],5)}')
            print(f'          T3G       ={np.round(T3g[(I,J)],5)}', flush=True)


if __name__ == '__main__':
    main()
