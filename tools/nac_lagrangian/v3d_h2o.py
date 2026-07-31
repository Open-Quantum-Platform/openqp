"""v3d TRIANGULATION (H2O): pin down, in ONE run,
  (T1) the Sk/S^x export convention   [FD cross-overlap vs C^T dsk C]
  (T2) PT-dX vs actual FD-dX          [per state/coord, eigencomponent split]
  (T3) my ampdir == A8's d_amp        [metric-function replication]
  (T4) v4 reproduction: ampdir(dX_FD) + gam:(Sk_FD + Ux) vs d_num
       and the analytic-Sk variant    [isolates the Sk substitution]

Run: python v3d_h2o.py <input.inp> <sweep.npz> <dnum.npz> <skeldir>
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

    inp, sweep_npz, dnum_npz, skeldir = (sys.argv[1], sys.argv[2],
                                         sys.argv[3], sys.argv[4])
    sw = np.load(sweep_npz)
    Ux_ref = sw['Ux']
    d_amp_ref = sw['d_amp']
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v3d.log'))
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
                    v[ijlr1 - 1] = (x[i - 1, jj - nocb - 1]) / RS
                elif ij == ijlr2:
                    pass
                else:
                    v[ij - 1] = x[i - 1, jj - nocb - 1]
        return v

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
    Sk_an = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        dsk = dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        Sk_an[c] = C.T @ dsk @ C

    # full A, physical subspace
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
    print(f'|A-A^T|max={np.abs(A-A.T).max():.3e}  '
          f'|A[:,lr2]|={np.linalg.norm(A[:, ijlr2-1]):.3e}  '
          f'|A[lr2,:]|={np.linalg.norm(A[ijlr2-1, :]):.3e}')
    phys = [k for k in range(nij) if k != ijlr2 - 1]
    Ap = 0.5 * (A[np.ix_(phys, phys)] + A[np.ix_(phys, phys)].T)
    w_full, Vp = np.linalg.eigh(Ap)
    V = np.zeros((nij, nij - 1))
    V[phys, :] = Vp
    for s in range(nstate):
        res = np.linalg.norm(Ap @ Xf[phys, s] - Om[s] * Xf[phys, s])
        print(f'|A_phys X_{s+1} - om X_{s+1}| = {res:.3e}')
    ks = [int(np.argmax(np.abs(V.T @ Xf[:, s]))) for s in range(nstate)]
    print('eigidx:', ks, 'low eigs:', np.round(w_full[:6], 6), flush=True)

    # w vectors: skel (reuse worker npz) + rot (in-process)
    w_skel = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        fp = np.load(os.path.join(skeldir, f'p{c}.npz'))
        fm = np.load(os.path.join(skeldir, f'p{c + ncoord}.npz'))
        w_skel[c] = (fp['Ax'] - fm['Ax']) / (2 * H)
    print('rotation FD sweep...', flush=True)
    w_rot = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        for sgn in (1.0, -1.0):
            Wp = (np.eye(nbf) + sgn * TROT * Ux_ref[c].T) @ W0
            mol.data['OQP::VEC_MO_A'] = Wp
            mol.data['OQP::VEC_MO_B'] = Wp
            for s in range(nstate):
                w_rot[c, s] += sgn * matvec(Xf[:, s]) / (2 * TROT)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw

    dX_PT = np.zeros((ncoord, nstate, nij))
    for J in range(nstate):
        for c in range(ncoord):
            w = w_skel[c, J] + w_rot[c, J]
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
            dX_PT[c, J] = acc

    # formula metric pieces
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

    # ---- displaced sweep (A8-style, in-process) --------------------------
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
        if np.any(sg < 0):
            print('  NOTE: orbital sign flips at displaced geometry')
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
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Xd, Skf

    print('displaced sweep...', flush=True)
    Sk_FD = np.zeros((ncoord, nbf, nbf))
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    dX_FD = np.zeros((ncoord, nstate, noca, nvirb))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += H
        Mp, Xp, Skp = displaced(cp)
        cm = xyz0.copy(); cm[c] -= H
        Mm, Xm2, Skm = displaced(cm)
        T = (Mp - Mm) / (2 * H)
        Sk_FD[c] = (Skp - Skm) / (2 * H)
        Ux_FD[c] = T - Sk_FD[c]
        for s in range(nstate):
            dX_FD[c, s] = (Xp[s] - Xm2[s]) / (2 * H)

    # ---- T1: Sk convention -----------------------------------------------
    print('\n===== T1: Sk_FD vs C^T dsk C (export) =====')
    n1 = np.linalg.norm(Sk_FD.ravel())
    for tag, v in (('as-exported', Sk_an),
                   ('transposed', np.transpose(Sk_an, (0, 2, 1)))):
        dv = np.abs(v - Sk_FD).max()
        cc = float(np.sum(v * Sk_FD)) / (np.linalg.norm(v.ravel()) * n1 + 1e-300)
        print(f'  {tag}: cos={cc:+.6f} maxdiff={dv:.3e} '
              f'(|Sk_FD|max={np.abs(Sk_FD).max():.3e})')
    print(f'  Ux consistency: max|Ux_FD - Ux_npz| = '
          f'{np.abs(Ux_FD - Ux_ref).max():.3e}')

    # ---- T2: PT-dX vs FD-dX ----------------------------------------------
    print('\n===== T2: dX_PT vs dX_FD =====')
    for J in range(nstate):
        for c in range(ncoord):
            fd = dX_FD[c, J]
            if np.linalg.norm(fd) < 1e-6:
                continue
            pt = unfold_vec(dX_PT[c, J])
            cc = float(np.sum(fd * pt)) / (np.linalg.norm(fd)
                                           * np.linalg.norm(pt) + 1e-300)
            line = (f'J={J+1} c={c}: |FD|={np.linalg.norm(fd):.6f} '
                    f'|PT|={np.linalg.norm(pt):.6f} cos={cc:+.4f}')
            diff = refold(fd) - dX_PT[c, J]
            cf = V.T @ diff
            top = np.argsort(-np.abs(cf))[:3]
            line += '  top-mismatch eig: ' + ', '.join(
                f'k={k}(w={w_full[k]:.4f},dc={cf[k]:+.4f})' for k in top)
            print(line)
    print('', flush=True)

    # ---- T3: ampdir replication ------------------------------------------
    print('===== T3: ampdir(dX_FD) vs frozen d_amp =====')
    amp_fd = np.zeros((nstate, nstate, ncoord))
    for J in range(nstate):
        for c in range(ncoord):
            g = ampdir_unf(J, dX_FD[c, J])
            for I in range(nstate):
                if I != J:
                    amp_fd[I, J, c] = g[I]
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            ref = d_amp_ref[:, I, J]
            v = amp_fd[I, J]
            cc = float(np.dot(ref, v)) / (np.linalg.norm(ref)
                                          * np.linalg.norm(v) + 1e-300)
            print(f'pair ({I+1},{J+1}): |ref|={np.linalg.norm(ref):.6f} '
                  f'|mine|={np.linalg.norm(v):.6f} cos={cc:+.6f} '
                  f'maxdiff={np.abs(v-ref).max():.3e}')

    # ---- T4: final assemblies --------------------------------------------
    print('\n===== T4: d_pred vs d_num =====')
    for lbl, ampsrc, Sksrc in (('FDdX+FDSk ', amp_fd, Sk_FD),
                               ('FDdX+anSk ', amp_fd, Sk_an)):
        for I in range(nstate):
            for J in range(nstate):
                if I >= J:
                    continue
                dn = dcv_n[I, J].reshape(-1)
                csf = np.array([np.sum(gam[I, J] * (Sksrc[c] + Ux_FD[c]))
                                for c in range(ncoord)])
                dp = ampsrc[I, J] + csf
                cc = float(np.dot(dn, dp)) / (np.linalg.norm(dn)
                                              * np.linalg.norm(dp) + 1e-300)
                print(f'[{lbl}] pair ({I+1},{J+1}): |d_num|={np.linalg.norm(dn):.6f} '
                      f'|pred|={np.linalg.norm(dp):.6f} cos={cc:+.6f} '
                      f'maxdiff={np.abs(dp-dn).max():.3e}')
        print('', flush=True)


if __name__ == '__main__':
    main()
