"""DIAGNOSTIC v2: + wsx term (-Tr[W^IJ S^x], the large-esum canceller),
+ empirical certification of the 3-D gamma push/read conventions via the
engine's own OQP::nac_trden_mo echo, + d_num caching.
Terms: [amp2e/gap, esum/gap, wsx/gap, zB(L)/gap, zB(gap*gam)/gap, ov].
"""
import os
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_d2.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E = list(mol.energies)
    Om = [E[k + 1] - E[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    noca, nocb, nbf, nij = ctx['noca'], ctx['nocb'], ctx['nbf'], ctx['nij']
    X0_raw = ctx['X0_raw']

    cache = inp.replace('.inp', '_dnum.npz')
    if os.path.exists(cache):
        dcv_n = np.load(cache)['dcv']
        print('d_num loaded from cache', flush=True)
    else:
        nac = NAC(mol)
        nacv_n, dcv_n, flags = nac.numerical_nac()
        print('numerical flags:', set(flags), flush=True)
        np.savez(cache, dcv=dcv_n)
        mol.data['OQP::td_bvec_mo'] = X0_raw

    oqp.mrsf_nac_amp(mol)
    amp2e = np.array(mol.data['OQP::nac_amp'], copy=True
                     ).reshape((nstate, nstate, natom, 3))
    esum_v, wsx_v = {}, {}
    for i in range(1, nstate + 1):
        for j in range(1, nstate + 1):
            if i == j:
                continue
            oqp.mrsf_nac_esum(mol, i, j)
            esum_v[(i, j)] = np.array(mol.data['OQP::nac_esum'],
                                      copy=True).reshape(-1)
            wsx_v[(i, j)] = np.array(mol.data['OQP::nac_wsx'],
                                     copy=True).reshape(-1)

    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    try:
        del mol.data['OQP::nac_gamma_tlf']
    except Exception:
        pass
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    ovF = np.array(mol.data['OQP::nac_overlap'], copy=True
                   ).reshape(-1).reshape((nstate, nstate, ncoord))
    ov = np.transpose(ovF, (1, 0, 2))

    # ---- CONVENTION CERTIFICATION: what gamma did the engine actually use?
    ech_raw = np.array(mol.data['OQP::nac_trden_mo'], copy=True)
    print('trden echo raw numpy shape:', ech_raw.shape, flush=True)
    ech = ech_raw.reshape(-1)
    print('\n===== gamma push/read certification =====')
    for a in range(nstate):
        for b in range(nstate):
            if a == b:
                continue
            # try echo layout: F-flat block (ist=a+1, jst=b+1)
            blk = ech[(a + b * nstate) * nbf * nbf:(a + b * nstate + 1) * nbf * nbf]
            Gecho = blk.reshape((nbf, nbf)).T   # F-flat -> matrix
            best = None
            for I in range(nstate):
                for J in range(nstate):
                    if I == J:
                        continue
                    for tag, M in (('g', gam[I, J]), ('gT', gam[I, J].T)):
                        dv = np.abs(Gecho - M).max()
                        if best is None or dv < best[0]:
                            best = (dv, I, J, tag)
            print(f'echo block(ist={a+1},jst={b+1}): best match '
                  f'gam[{best[1]},{best[2]}]{best[3]} maxdiff={best[0]:.2e}',
                  flush=True)

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

    def zB_of(Lmat, target):
        mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        mol.data.set_tdhf_target(target)
        oqp.set_mrsf_nac_cphf(mol, target, target % nstate + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        if not bool(mol.mol_energy.Z_Vector_converged):
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            return None
        gZ = grad_now()
        mol.data['OQP::td_p'] = np.zeros_like(
            np.array(mol.data['OQP::td_p'], copy=True))
        mol.data['OQP::WAO'] = np.zeros_like(
            np.array(mol.data['OQP::WAO'], copy=True))
        gS = grad_now()
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return gZ - gS

    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    rhs_single = {s: rhs_of(s + 1, Xf[:, s]) for s in range(nstate)}

    print('\n===== TERM DIAGNOSTIC v2 (upper orientation) =====', flush=True)
    for I in range(nstate):
        for J in range(I + 1, nstate):
            gap = Om[J] - Om[I]
            dn = dcv_n[I, J].reshape(-1)
            Lpol = 0.5 * (rhs_of(I + 1, Xf[:, I] + Xf[:, J])
                          - rhs_single[I] - rhs_of(I + 1, Xf[:, J]))
            zL = zB_of(unpack_rot_to_mat(Lpol), I + 1)
            zG = zB_of(gap * gam[I, J], I + 1)
            if zL is None or zG is None:
                print(f'pair ({I+1},{J+1}): z NOT converged', flush=True)
                continue
            terms = {
                'amp': amp2e[J, I].reshape(-1) / gap,
                'esum': esum_v[(I + 1, J + 1)] / gap,
                'wsx': wsx_v[(I + 1, J + 1)] / gap,
                'zL': zL / gap,
                'zG': zG / gap,
                'ov': ov[I, J],
            }
            print(f'\npair ({I+1},{J+1}) gap={gap:.6f} |d_num|={np.linalg.norm(dn):.6f}')
            for k, v in terms.items():
                c = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                            * np.linalg.norm(v) + 1e-300)
                print(f'  {k:5s}: |v|={np.linalg.norm(v):.6f} cos(d_num)={c:+.4f}')
            A = np.stack(list(terms.values()), axis=1)
            coef, *_ = np.linalg.lstsq(A, dn, rcond=None)
            fit = A @ coef
            print(f'  LSQ coef [amp,esum,wsx,zL,zG,ov] = '
                  + np.array2string(coef, precision=4, suppress_small=True)
                  + f'  resid={np.linalg.norm(dn-fit):.3e}')
            # sign-scan over z terms only, unit coefficients elsewhere
            base = terms['amp'] + terms['esum'] + terms['wsx'] + terms['ov']
            for sL in (1, -1):
                for sG in (1, -1):
                    s = base + sL * terms['zL'] + sG * terms['zG']
                    cc = float(np.dot(dn, s)) / (np.linalg.norm(dn)
                                                 * np.linalg.norm(s) + 1e-300)
                    print(f'  unit-sum sL={sL:+d} sG={sG:+d}: |v|={np.linalg.norm(s):.6f} '
                          f'cos={cc:+.6f} resid={np.linalg.norm(dn-s):.3e}')
            print('', flush=True)


if __name__ == '__main__':
    main()
