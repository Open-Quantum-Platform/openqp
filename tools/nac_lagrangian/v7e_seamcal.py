"""v7e SEAM CALIBRATION: inject UNIT-L matrices into the orbgrad hook.
For L = e_(p,q) (single 1.0), z = A_orb^{-1} pack(L) and the seam
z.B^x must equal the true CPHF response U^x_(p,q) coordinate by
coordinate. Compare against the FD orbital response (antisym part and
full), for several rotation pairs and both orientations.

Run: python v7e_seamcal.py <input.inp>
"""
import os
import sys
import numpy as np

H = 1e-3


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7e.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom

    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb = ctx['noca'], ctx['nocb']
    X0_raw = ctx['X0_raw']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a_r = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b_r = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)

    # ---- FD orbital response ---------------------------------------------
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
              'OQP::E_MO_B']}

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

    print('FD orbital-response sweep...', flush=True)
    Ux_FD = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        cp = xyz0.copy(); cp[c] += H
        Mp, Skp = displaced_M(cp)
        cm = xyz0.copy(); cm[c] -= H
        Mm, Skm = displaced_M(cm)
        Ux_FD[c] = (Mp - Mm) / (2 * H) - (Skp - Skm) / (2 * H)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw

    # ---- unit-L seam probes ----------------------------------------------
    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    def seam_of(Lmat):
        mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        mol.data.set_tdhf_target(1)
        oqp.set_mrsf_nac_cphf(mol, 1, 2)
        oqp.tdhf_mrsf_z_vector(mol)
        conv = bool(mol.mol_energy.Z_Vector_converged)
        gZ = grad_now()
        mol.data['OQP::td_p'] = np.zeros_like(
            np.array(mol.data['OQP::td_p'], copy=True))
        mol.data['OQP::WAO'] = np.zeros_like(
            np.array(mol.data['OQP::WAO'], copy=True))
        gS = grad_now()
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return gZ - gS, conv

    probes = [
        (nocb, 0),           # socc1 <- doc1   (ds block)
        (nocb + 1, nocb - 1),  # socc2 <- doc4 (ds block)
        (noca, 0),           # virt1 <- doc1   (dv block)
        (noca + 2, nocb - 1),  # virt3 <- doc4 (dv block)
        (noca, nocb),        # virt1 <- socc1  (sv block)
        (noca + 1, nocb + 1),  # virt2 <- socc2 (sv block)
    ]
    print('\n===== UNIT-L SEAM CALIBRATION =====')
    print('seam(L=e_pq) should equal U^x_pq(c); compare vs FD U (antisym/full)')
    for (p, q) in probes:
        L = np.zeros((nbf, nbf))
        L[p, q] = 1.0
        s, conv = seam_of(L)
        ua = np.array([0.5 * (Ux_FD[c][p, q] - Ux_FD[c][q, p])
                       for c in range(ncoord)])
        uf = np.array([Ux_FD[c][p, q] for c in range(ncoord)])
        r_a = (float(np.dot(s, ua)) / (np.dot(ua, ua) + 1e-300))
        r_f = (float(np.dot(s, uf)) / (np.dot(uf, uf) + 1e-300))
        ca = float(np.dot(s, ua)) / (np.linalg.norm(s)
                                     * np.linalg.norm(ua) + 1e-300)
        print(f'L=e({p},{q}) conv={conv}: |seam|={np.linalg.norm(s):.6f} '
              f'|Ua|={np.linalg.norm(ua):.6f} |Uf|={np.linalg.norm(uf):.6f} '
              f'cos(s,Ua)={ca:+.5f} fit s=r*Ua: r={r_a:+.4f}  '
              f'fit s=r*Uf: r={r_f:+.4f}')
        # also the transpose injection
        L = np.zeros((nbf, nbf))
        L[q, p] = 1.0
        s2, conv2 = seam_of(L)
        print(f'   e({q},{p}): |seam|={np.linalg.norm(s2):.6f} '
              f'(vs -e({p},{q}): maxdiff={np.abs(s2 + s).max():.2e})')
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass


if __name__ == '__main__':
    main()
