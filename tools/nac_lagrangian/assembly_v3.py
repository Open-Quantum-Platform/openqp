"""ASSEMBLY GATE v3 -- the exact-term Lagrangian assembly with certified
ingredients, no fitted constants:

  d_IJ = [skel_IJ + M^IJ : U^x]/gap_IJ  +  gam^IJ : (Sk_MO + U^x)

  skel = amp2e + esum            (frozen-MO FD-certified pair)
  M    = La + Ls                 (full unrestricted orbital derivative of
                                  E_IJ = X_I^T A X_J, from the frozen
                                  generator sweeps; La antisym, Ls sym)
  U^x  = frozen FD orbital response (A8/A10 npz; production replaces this
                                  by z-vector + canonical chain)
  Sk_MO= C^T dSket_AO C          (NEW analytic export, NAC_DUMP_DS)
  gam  = closed-form gamma^formula

Stacked checks:
  C1  sym(U^x) == -1/2 S^x_MO(full)      [export/Ux cross-certification]
  C2  [skel + M:U^x]/gap == d_amp (npz)  [the PT identity, per coordinate]
  C3  d_pred vs d_num (frozen)           [final, SIGNED]

Run: python assembly_v3.py <input.inp> <sweep.npz> <dnum.npz>
"""
import os
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, sweep_npz, dnum_npz = sys.argv[1], sys.argv[2], sys.argv[3]
    sw = np.load(sweep_npz)
    La, Ls, Ux = sw['La'], sw['Ls'], sw['Ux']
    d_amp_ref = sw['d_amp']                     # (ncoord, nstate, nstate)
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f[dn_f.files[0] if 'dcv' not in dn_f.files else 'dcv']

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v3.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E0 = list(mol.energies)
    Om = [E0[k + 1] - E0[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    nbf = ctx['nbf']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    C = W0.T                                    # C[ao, mo]

    # skeleton engines
    oqp.mrsf_nac_amp(mol)
    amp2e = np.array(mol.data['OQP::nac_amp'], copy=True
                     ).reshape((nstate, nstate, natom, 3))
    esum_v = {}
    for i in range(1, nstate + 1):
        for j in range(1, nstate + 1):
            if i == j:
                continue
            oqp.mrsf_nac_esum(mol, i, j)
            esum_v[(i, j)] = np.array(mol.data['OQP::nac_esum'],
                                      copy=True).reshape(-1)

    # gamma + dSket/dSfull exports (gamma push only to satisfy the engine;
    # the ov output itself is unused here)
    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    # raw buffer = Fortran flat of (nbf^2, 3N): contiguous nbf^2 blocks/coord
    dsk_raw = np.array(mol.data['OQP::dbg_dsket'], copy=True).reshape(-1)
    dsf_raw = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    nbf2 = nbf * nbf
    Sk_MO = np.zeros((ncoord, nbf, nbf))
    Sx_MO = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        dsk = dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        dsf = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        Sk_MO[c] = C.T @ dsk @ C
        Sx_MO[c] = C.T @ dsf @ C

    # ---- C1: sym(Ux) vs -1/2 Sx_MO ---------------------------------------
    print('===== C1: sym(U^x) == -1/2 S^x_MO =====')
    worst = 0.0
    for c in range(ncoord):
        sym = 0.5 * (Ux[c] + Ux[c].T)
        dv = np.abs(sym + 0.5 * Sx_MO[c]).max()
        worst = max(worst, dv)
    print(f'max deviation over {ncoord} coords: {worst:.3e}')
    # scale context
    print(f'scale |S^x|max = {np.abs(Sx_MO).max():.3e}', flush=True)

    # ---- build M and the amplitude part ----------------------------------
    M = np.zeros_like(La)
    for I in range(nstate):
        for J in range(nstate):
            M[I, J] = La[I, J] + Ls[I, J]
            np.fill_diagonal(M[I, J], np.diag(Ls[I, J]))

    print('\n===== C2: [skel + M:U^x]/gap vs d_amp (frozen FD) =====')
    amp_pred = np.zeros((nstate, nstate, ncoord))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            gap = Om[J] - Om[I]
            sk = amp2e[J, I].reshape(-1) + esum_v[(I + 1, J + 1)]
            mu = np.array([np.sum(M[I, J] * Ux[c]) for c in range(ncoord)])
            amp_pred[I, J] = (sk + mu) / gap
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
            # reverse orientation independently
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
