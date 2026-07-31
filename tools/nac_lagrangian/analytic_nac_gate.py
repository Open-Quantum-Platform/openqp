"""PRODUCTION-CANDIDATE GATE: the fully-analytic MRSF NAC assembly
(combined-RHS prescription, MRSF_NAC_DERIVATION.md 4.3 + 7.23):

  h_IJ  = amp2e + esum + z.B^x,   A^orb z = [L_polarized + gap*gamma_a]|rot
  d_IJ  = antisym{ h_IJ/gap + gamma_a : Sk-contraction (mrsf_nac_overlap) }

No finite differences anywhere except inside the certified Fortran
engines. Compared SIGNED against the pristine production numerical NAC.

Run:  python analytic_nac_gate.py <input.inp>
"""
import os
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_an.log'))
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

    # ---- 1) pristine numerical reference ----------------------------------
    nac = NAC(mol)
    nacv_n, dcv_n, flags = nac.numerical_nac()
    print('numerical flags:', set(flags), flush=True)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- 2) analytic skeleton engines -------------------------------------
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
    print('skeleton engines done.', flush=True)

    # ---- 3) closed-form gamma^formula + Sk contraction --------------------
    gam = FK.gamma_closed(ctx)          # antisym in (p,q) by construction
    print('gamma^closed done.', flush=True)
    # push as OQP::nac_gamma_tlf (F-layout: flat[k + nbf^2*ist + ...],
    # k = p + q*nbf -> C-flat of gamma^T per pair)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            blockG = gam[I, J]
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                blockG.T.reshape(-1)
    try:
        del mol.data['OQP::nac_gamma_tlf']
    except Exception:
        pass
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    ovraw = np.array(mol.data['OQP::nac_overlap'], copy=True)
    ovF = ovraw.reshape(-1).reshape((nstate, nstate, ncoord))
    ov = np.transpose(ovF, (1, 0, 2))    # ov[I-1, J-1, :] = Sk-contraction
    print('gamma:Sk contraction done.', flush=True)

    # ---- 4) combined-RHS z-vector per pair --------------------------------
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

    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    rhs_single = {s: rhs_of(s + 1, Xf[:, s]) for s in range(nstate)}
    zB = {}
    for I in range(nstate):
        for J in range(I + 1, nstate):
            gap = Om[J] - Om[I]
            Lpol = 0.5 * (rhs_of(I + 1, Xf[:, I] + Xf[:, J])
                          - rhs_single[I] - rhs_of(I + 1, Xf[:, J]))
            Lmat = unpack_rot_to_mat(Lpol) + gap * gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(I + 1)
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
            print(f'  ({I+1},{J+1}): |zB|={np.linalg.norm(gZ-gS):.6f}', flush=True)
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass

    # ---- 5) assembly + SIGNED comparison ----------------------------------
    print('\n========== ANALYTIC-NAC PRODUCTION-CANDIDATE GATE ==========')
    d_pred = np.zeros((nstate, nstate, ncoord))
    for I in range(nstate):
        for J in range(I + 1, nstate):
            gap = Om[J] - Om[I]
            if zB[(I, J)] is None:
                continue
            h = (amp2e[J, I].reshape(-1) + esum_v[(I + 1, J + 1)]
                 + zB[(I, J)])
            d_pred[I, J] = h / gap + ov[I, J]
            # the (J,I) orientation from the independently-computed pieces
            h_JI = (amp2e[I, J].reshape(-1) + esum_v[(J + 1, I + 1)]
                    + zB[(I, J)] * (-1.0))     # z is pair-symmetric in Ltot? report both
            d_pred[J, I] = h_JI / (-gap) + ov[J, I]
    d_a = 0.5 * (d_pred - np.transpose(d_pred, (1, 0, 2)))
    for I in range(nstate):
        for J in range(I + 1, nstate):
            dn = dcv_n[I, J].reshape(-1)
            dp = d_a[I, J]
            c = float(np.dot(dn, dp)) / (np.linalg.norm(dn)
                                         * np.linalg.norm(dp) + 1e-300)
            print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(dn):.6f} '
                  f'|d_pred|={np.linalg.norm(dp):.6f} SIGNED cos={c:+.6f} '
                  f'max|diff|={np.abs(dp-dn).max():.3e}')
            print(f'  num : {np.round(dn, 5)}')
            print(f'  pred: {np.round(dp, 5)}')


if __name__ == '__main__':
    main()
