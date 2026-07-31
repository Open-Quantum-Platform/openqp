"""v3f (H2O): localize the non-physical (gauge) content of the FD orbital
response and re-judge the continuous-gauge assembly with a cleaned U.

 (P1) MO-energy spectrum + same-space near-degeneracy map
 (P2) Ux at two step sizes: the physical U^x is h-independent; arbitrary
      degenerate-sector rotations are not. |Ux(h1)-Ux(h2)| maps the junk.
 (P3) cleaned U := Ux(h1) with the inconsistent entries zeroed;
      w_rot along U_clean; dX_PT; judge
      antisym[ampdir(dX_PT) + gam:(Sk_an + U_clean)] vs d_num.

Run: python v3f_h2o.py <input.inp> <dnum.npz> <skeldir>
"""
import os
import sys
import numpy as np

TROT = 1e-5
EPSA = 1e-5
H1 = 1e-3
H2 = 5e-4


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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v3f.log'))
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
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True).ravel()
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca

    print('P1: |W0-Wb0|max =', np.abs(W0 - Wb0).max())
    print('P1: eps_a =', np.round(e0a, 6))
    spaces = [range(0, nocb), range(nocb, noca), range(noca, nbf)]
    for sp in spaces:
        sp = list(sp)
        for i in range(len(sp)):
            for j in range(i + 1, len(sp)):
                if abs(e0a[sp[i]] - e0a[sp[j]]) < 1e-5:
                    print(f'P1: near-degenerate pair ({sp[i]},{sp[j]}) '
                          f'de={e0a[sp[i]]-e0a[sp[j]]:.2e}')
    print('', flush=True)

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
        mol.data['OQP::E_MO_A_old'] = e0a.reshape(-1)
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
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Xd, Skf

    def sweep(h):
        T = np.zeros((ncoord, nbf, nbf))
        Sk = np.zeros((ncoord, nbf, nbf))
        dX = np.zeros((ncoord, nstate, noca, nvirb))
        for c in range(ncoord):
            cp = xyz0.copy(); cp[c] += h
            Mp, Xp, Skp = displaced(cp)
            cm = xyz0.copy(); cm[c] -= h
            Mm, Xm2, Skm = displaced(cm)
            T[c] = (Mp - Mm) / (2 * h)
            Sk[c] = (Skp - Skm) / (2 * h)
            for s in range(nstate):
                dX[c, s] = (Xp[s] - Xm2[s]) / (2 * h)
        return T, Sk, dX

    print('displaced sweep h1...', flush=True)
    T1, Sk1, dX1 = sweep(H1)
    print('displaced sweep h2...', flush=True)
    T2, Sk2, dX2 = sweep(H2)
    U1 = T1 - Sk1
    U2 = T2 - Sk2

    print('\n===== P2: Ux(h1) vs Ux(h2) block map =====')
    dU = np.abs(U1 - U2).max(axis=0)      # (nbf, nbf) worst over coords
    blocks = {'dd': (range(0, nocb), range(0, nocb)),
              'ds': (range(nocb, noca), range(0, nocb)),
              'dv': (range(noca, nbf), range(0, nocb)),
              'ss': (range(nocb, noca), range(nocb, noca)),
              'sv': (range(noca, nbf), range(nocb, noca)),
              'vv': (range(noca, nbf), range(noca, nbf))}
    for k, (rows, cols) in blocks.items():
        sub = dU[np.ix_(list(rows), list(cols))]
        print(f'  {k}: max|U(h1)-U(h2)|={sub.max():.3e}  '
              f'|U(h1)|max={np.abs(U1.max(axis=0)[np.ix_(list(rows), list(cols))]).max():.3e}')
    idx = np.argsort(-dU.ravel())[:8]
    for f in idx:
        p, q = divmod(int(f), nbf)
        print(f'  worst ({p},{q}): dU={dU[p,q]:.3e} eps_p={e0a[p]:.4f} '
              f'eps_q={e0a[q]:.4f}')
    print('', flush=True)

    # cleaned U: zero entries where the two step sizes disagree materially
    mask = dU > max(0.02 * np.abs(U1).max(), 5 * np.median(dU) + 1e-12)
    print(f'cleaning mask: {int(mask.sum())} of {nbf*nbf} entries zeroed')
    U_cl = U1.copy()
    for c in range(ncoord):
        U_cl[c][mask] = 0.0

    # restore reference, rot FD along U_cl
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw
    print('rotation FD sweep along U_clean...', flush=True)
    w_skel = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        fp = np.load(os.path.join(skeldir, f'p{c}.npz'))
        fm = np.load(os.path.join(skeldir, f'p{c + ncoord}.npz'))
        w_skel[c] = (fp['Ax'] - fm['Ax']) / (2 * H1)
    w_rot = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        for sgn in (1.0, -1.0):
            Wp = (np.eye(nbf) + sgn * TROT * U_cl[c].T) @ W0
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

    print('\n===== dX_PT(U_clean) vs dX_FD(h1) =====')
    for J in range(nstate):
        na = 0.0
        for c in range(ncoord):
            fd = dX1[c, J]
            pt = unfold_vec(dX_PT[c, J])
            if np.linalg.norm(fd) < 1e-6:
                continue
            cc = float(np.sum(fd * pt)) / (np.linalg.norm(fd)
                                           * np.linalg.norm(pt) + 1e-300)
            na = max(na, np.abs(pt - fd).max())
            if c in (0, 5):
                print(f'J={J+1} c={c}: |FD|={np.linalg.norm(fd):.6f} '
                      f'|PT|={np.linalg.norm(pt):.6f} cos={cc:+.5f}')
        print(f'J={J+1}: max|PT-FD| over coords = {na:.4f}')
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

    judge('J1 sanity v4-repro (FD h1 pair)', lambda c, J: dX1[c, J], T1)
    judge('J2 continuous: amp(PT)+gam:(Sk_an+U_clean)',
          lambda c, J: unfold_vec(dX_PT[c, J]), Sk_an + U_cl)
    judge('J3 mixed-check: amp(FD)+gam:(Sk_an+U_clean)',
          lambda c, J: dX1[c, J], Sk_an + U_cl)


if __name__ == '__main__':
    main()
