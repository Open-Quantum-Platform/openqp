"""v17: THE ONE-LINE FIX -- Sk from the RAW PRODUCT C0^T S_cross C0
(numpy, from the overlap_ao record) instead of the routine output.
 (M1) reconstruction: |C0(1 + h u+) - C+_signfixed|  (raw coefficients)
 (M2) pure C-channel h-scaling: for h in {1e-3, 5e-4}:
        lhs(h) = [A(C'_sf, INTS_d) - A(C0, INTS_d)] X / h   (one-sided)
        vs staged(u+(h)) at reference integrals
      If |lhs - staged| halves with h -> mixed d2A/dxdtheta (benign
      finite-h referee artifact) -> the paradox dissolves.
Run: python v13_final.py <inp> <coord>
"""
import os, sys
import numpy as np

def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp = sys.argv[1]
    c0 = int(sys.argv[2])
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v17.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    X0_raw = ctx['X0_raw']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

    J = 2
    XJ = Xf[:, J]
    base = matvec(XJ)
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False

    def restore():
        mol.update_system(xyz0)
        oqp.library.ints_1e(mol)
        for k, v in SAVE0.items():
            mol.data[k] = v.copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass

    def one_sided(h):
        cp = xyz0.copy(); cp[c0] += h
        mol.update_system(cp)
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
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        S_np = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
        S_f = S_np.reshape(-1).reshape((nbf, nbf)).T
        Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Wd_sf = sg[:, None] * Wd
        mol.data['OQP::VEC_MO_A'] = Wd_sf
        mol.data['OQP::VEC_MO_B'] = Wd_sf.copy()
        Ax = matvec(XJ)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        Ax0 = matvec(XJ)
        C0m = W0.T
        Skf = C0m.T @ S_f @ C0m       # THE ONE-LINE FIX (raw product)
        return M_f, Skf, Ax, Ax0, Wd_sf

    print('===== v13: the two named measurements =====')
    for h in (1e-3, 5e-4):
        M_f, Skf, Ax, Ax0, Wd_sf = one_sided(h)
        u = (M_f - Skf) / h
        # M1 reconstruction (raw): C0(1+h u) vs C'_sf  -> W-side:
        # (C0(1+hu))^T = (1 + h u^T) W0
        W_rec = (np.eye(nbf) + h * u.T) @ W0
        rec = np.abs(W_rec - Wd_sf).max()
        lhs = (Ax - Ax0) / h
        restore()
        t = 1e-5
        Wp = (np.eye(nbf) + t * u.T) @ W0
        mol.data['OQP::VEC_MO_A'] = Wp
        mol.data['OQP::VEC_MO_B'] = Wp
        a = matvec(XJ)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        staged = (a - base) / t
        d = lhs - staged
        print(f'h={h:.0e}: M1 reconstruction max|C0(1+hu) - C_sf| = {rec:.4e}')
        print(f'         M2 |lhs - staged| = {np.linalg.norm(d):.5f} '
              f'maxdiff={np.abs(d).max():.4e}')
        restore()

if __name__ == '__main__':
    main()
