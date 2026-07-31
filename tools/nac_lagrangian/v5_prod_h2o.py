"""v5 PRODUCTION-STRUCTURE GATE: the no-FD-per-coordinate assembly.

  d^c[I,J] = amp2e(ytil,X_J)^c + esum(ytil,X_J)^c + wsx(ytil,X_J)^c
           + zB[L_a(ytil,X_J) + gamma_a]^c + gamma:Sk_an^c
  then antisym over (I,J).

Ingredients:
  ytil      : (om_J - A)^{-1}|perp G_met  (eigen for the gate; MINRES
              certified on one pair)
  engines   : SLOT INJECTION -- push ytil into state-I's bvec slot; the
              bilinear Fortran engines then compute the (ytil, X_J) pair
              objects with no new code
  seam      : combined-RHS z-vector via the orbgrad hook (Handy-Schaefer;
              the Fock density response lives in the LHS -- the channel
              my manual w missed)
  residual  : whatever is left = same-space canonical + 2e-W content
              (measured, then implemented if material)

Gate G-A: ytil.w_skel^c == [amp2e+esum](ytil,X_J)^c  (engine generality)
Run: python v5_prod_h2o.py <input.inp> <dnum.npz>
"""
import os
import sys
import subprocess
import numpy as np

EPSA = 1e-5
H = 1e-3
NPAR = 12
WOMP = '2'


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, dnum_npz = sys.argv[1], sys.argv[2]
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v5.log'))
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

    # ---- G_met + ytil ----------------------------------------------------
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

    # MINRES certification on pair (1,2): solve (om_J - A)|perp y = G_met
    try:
        from scipy.sparse.linalg import minres, LinearOperator
        I0, J0 = 0, 1
        XJ = Xf[:, J0]

        def op(v):
            vv = np.zeros(nij)
            vv[phys] = v
            vv -= XJ * float(np.dot(XJ, vv))
            av = Om[J0] * vv - A @ vv
            av -= XJ * float(np.dot(XJ, av))
            return av[phys]

        rhs = G_met[I0, J0].copy()
        rhs -= XJ * float(np.dot(XJ, rhs))
        y, info = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                         rhs[phys], rtol=1e-10, maxiter=2000)
        yfull = np.zeros(nij)
        yfull[phys] = y
        yfull -= XJ * float(np.dot(XJ, yfull))
        print(f'MINRES vs eigen ytil (1,2): info={info} '
              f'maxdiff={np.abs(yfull - ytil[I0, J0]).max():.3e}', flush=True)
    except Exception as ex:
        print('MINRES check skipped:', ex, flush=True)

    # ---- regenerate skel workers in THIS process's phases ---------------
    skeldir = inp.replace('.inp', '') + '_skel'
    ref_npz = os.path.join(skeldir, 'ref.npz')
    np.savez(ref_npz, C_a=W0, C_b=Wb0, X0_raw=X0_raw,
             nstate=nstate, nij=nij)
    jobs = [(os.path.join(skeldir, f'p{idx}.inp'),
             os.path.join(skeldir, f'p{idx}.npz'))
            for idx in range(2 * ncoord)]
    print(f'rerunning {len(jobs)} skel workers (v5 phases)...', flush=True)
    env = dict(os.environ, OMP_NUM_THREADS=WOMP)
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'skel_gate.py')
    procs = []
    for pinp, out in jobs:
        while len([p for p in procs if p.poll() is None]) >= NPAR:
            import time
            time.sleep(0.5)
        procs.append(subprocess.Popen(
            [sys.executable, helper, '--worker', pinp, ref_npz, out],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()
    w_skel = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        fp = np.load(os.path.join(skeldir, f'p{c}.npz'))
        fm = np.load(os.path.join(skeldir, f'p{c + ncoord}.npz'))
        w_skel[c] = (fp['Ax'] - fm['Ax']) / (2 * H)

    # ---- slot-injection engines -----------------------------------------
    def inject(I, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[I * nij:(I + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)

    amp_e = {}
    esum_e = {}
    wsx_e = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            inject(I, ytil[I, J])
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            amp_e[(I, J)] = a[J, I].reshape(-1).copy()
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            esum_e[(I, J)] = np.array(mol.data['OQP::nac_esum'],
                                      copy=True).reshape(-1)
            wsx_e[(I, J)] = np.array(mol.data['OQP::nac_wsx'],
                                     copy=True).reshape(-1)
            mol.data['OQP::td_bvec_mo'] = X0_raw
    print('slot-injection engines done.', flush=True)

    print('\n===== G-A: engine generality: ytil.w_skel vs amp+esum(ytil,X) =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            lhs = np.array([float(np.dot(ytil[I, J], w_skel[c, J]))
                            for c in range(ncoord)])
            rhs = amp_e[(I, J)] + esum_e[(I, J)]
            cc = float(np.dot(lhs, rhs)) / (np.linalg.norm(lhs)
                                            * np.linalg.norm(rhs) + 1e-300)
            print(f'pair ({I+1},{J+1}): |ytil.w_skel|={np.linalg.norm(lhs):.6f} '
                  f'|engines|={np.linalg.norm(rhs):.6f} cos={cc:+.6f} '
                  f'maxdiff={np.abs(lhs-rhs).max():.3e}')
    print('', flush=True)

    # ---- combined-RHS z-vector seam --------------------------------------
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

    zB = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            qy = rhs_of(J + 1, ytil[I, J])
            qx = rhs_of(J + 1, Xf[:, J])
            qm = rhs_of(J + 1, ytil[I, J] + Xf[:, J])
            Lpol = 0.5 * (qm - qy - qx)
            Lmat = unpack_rot_to_mat(Lpol) + gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            if not bool(mol.mol_energy.Z_Vector_converged):
                oqp.set_mrsf_nac_cphf(mol, 0, 0)
                zB[(I, J)] = None
                print(f'  ({I+1},{J+1}): z NOT converged', flush=True)
                continue
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            zB[(I, J)] = gZ - gS
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass
    mol.data['OQP::td_bvec_mo'] = X0_raw
    print('seam solves done.', flush=True)

    # ---- assembly + judge -------------------------------------------------
    print('\n===== v5 PRODUCTION-STRUCTURE ASSEMBLY vs d_num =====')
    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J or zB[(I, J)] is None:
                continue
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            dp[:, I, J] = (amp_e[(I, J)] + esum_e[(I, J)] + wsx_e[(I, J)]
                           + zB[(I, J)] + gsk)
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
    np.savez(inp.replace('.inp', '_v5.npz'), dp=dp, dpa=dpa, ytil=ytil,
             G_met=G_met, w_skel=w_skel,
             amp_e=np.array([[amp_e.get((i, j), np.zeros(ncoord))
                              for j in range(nstate)] for i in range(nstate)]),
             esum_e=np.array([[esum_e.get((i, j), np.zeros(ncoord))
                               for j in range(nstate)] for i in range(nstate)]),
             wsx_e=np.array([[wsx_e.get((i, j), np.zeros(ncoord))
                              for j in range(nstate)] for i in range(nstate)]),
             zB=np.array([[zB.get((i, j)) if zB.get((i, j)) is not None
                           else np.zeros(ncoord)
                           for j in range(nstate)] for i in range(nstate)]))


if __name__ == '__main__':
    main()
