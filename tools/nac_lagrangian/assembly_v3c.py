"""ASSEMBLY GATE v3c -- the v4-certified pair (formula-metric amplitude
term + gamma^formula:T) with the FD ingredients replaced ANALYTICALLY:

  d_IJ^c = ampdir_J(dX_J^c)[I]  +  gamma^IJ : (Sk^c + U^x_c)

  dX_J^c = full-eigenbasis PT:  sum_{K != J} V_K (V_K^T w_J^c)/(om_J - w_K)
  w_J^c  = (dA X_J)^c = frozen-MO displaced FD (skel workers, Ax export)
                       + rotation-response FD along U^x_c (matvec staging)
  Sk^c   = C^T dSket_AO C   (analytic export)
  U^x_c  = frozen A8/A10 FD orbital response (production: z-vector chain)

Diagnostics: dsf==dsk+dsk^T, C1 both signs, raw-A symmetry/eigen residuals,
amp term vs frozen d_amp (the REAL C2), final SIGNED d vs d_num.

Run: python assembly_v3c.py <input.inp> <sweep.npz> <dnum.npz> <skeldir>
"""
import os
import sys
import subprocess
import numpy as np

TROT = 1e-5     # rotation-FD step
EPSA = 1e-5     # amplitude directional-FD step
NPAR = 12
WOMP = '2'


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, sweep_npz, dnum_npz, skeldir = (sys.argv[1], sys.argv[2],
                                         sys.argv[3], sys.argv[4])
    sw = np.load(sweep_npz)
    Ux = sw['Ux']
    d_amp_ref = sw['d_amp']
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v3c.log'))
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
    print(f'lr-slot check: X0[lr1]={Xf[ijlr1-1,:]} X0[lr2]={Xf[ijlr2-1,:]}',
          flush=True)

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

    # ---- gamma + dS exports ----------------------------------------------
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
    Sk_MO = np.zeros((ncoord, nbf, nbf))
    Sx_MO = np.zeros((ncoord, nbf, nbf))
    dev_sym = 0.0
    for c in range(ncoord):
        dsk = dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        dsf = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        dev_sym = max(dev_sym, np.abs(dsf - (dsk + dsk.T)).max())
        Sk_MO[c] = C.T @ dsk @ C
        Sx_MO[c] = C.T @ dsf @ C
    print(f'DIAG: max|dsf - (dsk+dsk^T)| = {dev_sym:.3e}')
    d_p = max(np.abs(0.5 * (Ux[c] + Ux[c].T) + 0.5 * Sx_MO[c]).max()
              for c in range(ncoord))
    d_m = max(np.abs(0.5 * (Ux[c] + Ux[c].T) - 0.5 * Sx_MO[c]).max()
              for c in range(ncoord))
    print(f'DIAG C1: |sym(Ux)+0.5Sx|={d_p:.3e}  |sym(Ux)-0.5Sx|={d_m:.3e} '
          f'(scale {np.abs(Sx_MO).max():.3e})', flush=True)

    # ---- full A build ----------------------------------------------------
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

    print(f'building full A ({nij}x{nij}) by matvec columns...', flush=True)
    A = np.zeros((nij, nij))
    for k in range(nij):
        e = np.zeros(nij)
        e[k] = 1.0
        A[:, k] = matvec(e)
    mol.data['OQP::td_bvec_mo'] = X0_raw
    print(f'  |A - A^T|max = {np.abs(A - A.T).max():.3e}')
    for s in range(nstate):
        res = np.linalg.norm(A @ Xf[:, s] - Om[s] * Xf[:, s])
        print(f'  |A X_{s+1} - om X_{s+1}| = {res:.3e}')
    As = 0.5 * (A + A.T)
    w_full, V = np.linalg.eigh(As)

    # ---- (dA X_s)_skel from displaced workers ----------------------------
    ref_npz = os.path.join(skeldir, 'ref.npz')
    if not os.path.exists(ref_npz):
        np.savez(ref_npz, C_a=W0, C_b=Wb0, X0_raw=X0_raw,
                 nstate=nstate, nij=nij)
    jobs = []
    for idx in range(2 * ncoord):
        out = os.path.join(skeldir, f'p{idx}.npz')
        need = True
        if os.path.exists(out):
            try:
                need = 'Ax' not in np.load(out).files
            except Exception:
                need = True
        if need:
            jobs.append((os.path.join(skeldir, f'p{idx}.inp'), out))
    print(f'rerunning {len(jobs)} displaced workers (Ax export)...', flush=True)
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
    H = 1.0e-3
    w_skel = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        fp = np.load(os.path.join(skeldir, f'p{c}.npz'))
        fm = np.load(os.path.join(skeldir, f'p{c + ncoord}.npz'))
        w_skel[c] = (fp['Ax'] - fm['Ax']) / (2 * H)

    # ---- (dA X_s)_rot: directional matvec FD along Ux^c ------------------
    print('rotation-response FD sweep...', flush=True)
    w_rot = np.zeros((ncoord, nstate, nij))
    for c in range(ncoord):
        for sgn in (1.0, -1.0):
            Wp = (np.eye(nbf) + sgn * TROT * Ux[c].T) @ W0
            mol.data['OQP::VEC_MO_A'] = Wp
            mol.data['OQP::VEC_MO_B'] = Wp
            for s in range(nstate):
                Ax = matvec(Xf[:, s])
                w_rot[c, s] += sgn * Ax / (2 * TROT)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- PT amplitude derivative + formula-metric amp term ---------------
    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))

    def ampdir(J, dXJ_vec):
        dXt = unfold_vec(dXJ_vec)
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    ks = [int(np.argmax(np.abs(V.T @ Xf[:, s]))) for s in range(nstate)]
    print('state eigenindices:', ks,
          'eig vs Om:', [f'{w_full[ks[s]]:.8f}/{Om[s]:.8f}'
                         for s in range(nstate)], flush=True)

    amp_pred = np.zeros((nstate, nstate, ncoord))
    for J in range(nstate):
        for c in range(ncoord):
            w = w_skel[c, J] + w_rot[c, J]
            coef = V.T @ w
            dX = np.zeros(nij)
            for k in range(nij):
                if k == ks[J]:
                    continue
                den = Om[J] - w_full[k]
                if abs(den) < 1e-10:
                    continue
                dX += V[:, k] * (coef[k] / den)
            g = ampdir(J, dX)
            for I in range(nstate):
                if I != J:
                    amp_pred[I, J, c] = g[I]

    print('\n===== C2 (real): amp term vs frozen d_amp =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            ref = d_amp_ref[:, I, J]
            v = amp_pred[I, J]
            cc = float(np.dot(ref, v)) / (np.linalg.norm(ref)
                                          * np.linalg.norm(v) + 1e-300)
            print(f'pair ({I+1},{J+1}): |d_amp|={np.linalg.norm(ref):.6f} '
                  f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.6f} '
                  f'maxdiff={np.abs(v-ref).max():.3e}')

    print('\n===== C3: FULL d_pred vs d_num (SIGNED) =====')
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            dn = dcv_n[I, J].reshape(-1)
            csf = np.array([np.sum(gam[I, J] * (Sk_MO[c] + Ux[c]))
                            for c in range(ncoord)])
            dp = amp_pred[I, J] + csf
            csf_r = np.array([np.sum(gam[J, I] * (Sk_MO[c] + Ux[c]))
                              for c in range(ncoord)])
            dp_r = amp_pred[J, I] + csf_r
            da = 0.5 * (dp - dp_r)
            for tag, v in (('direct', dp), ('antisym', da)):
                cc = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                             * np.linalg.norm(v) + 1e-300)
                print(f'pair ({I+1},{J+1}) [{tag:7s}]: |d_num|={np.linalg.norm(dn):.6f} '
                      f'|d_pred|={np.linalg.norm(v):.6f} cos={cc:+.6f} '
                      f'maxdiff={np.abs(v-dn).max():.3e}')
            print(f'  d_num: {np.round(dn, 5)}')
            print(f'  pred : {np.round(dp, 5)}', flush=True)


if __name__ == '__main__':
    main()
