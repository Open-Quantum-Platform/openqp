"""v14: STATE-SWAP bisection. Inside the displaced state (before any
restore), evaluate the staged pair [A(C0(1+h u)) - A(C0)]/h and compare
with lhs = [A(C+_sf) - A(C0)]/h computed in the SAME state.
Expected == to ~3e-4 (reconstruction). Then restore to the reference
state and evaluate the same staged pair -> the h^0 state shift.
Run: python v14_swap.py <inp> <coord>
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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v14.log'))
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
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    base_ref = matvec(XJ)
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False
    HD = 1e-3

    cp = xyz0.copy(); cp[c0] += HD
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
    Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wd_sf = sg[:, None] * Wd
    # frozen-C reference within the displaced state
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    Ax0_d = matvec(XJ)
    # true displaced C
    mol.data['OQP::VEC_MO_A'] = Wd_sf
    mol.data['OQP::VEC_MO_B'] = Wd_sf.copy()
    Ax_d = matvec(XJ)
    # u and the staged pair INSIDE the displaced state
    oqp.get_structures_ao_overlap(mol)  # needs VEC_MO = ? careful: Sk with W0
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    oqp.get_structures_ao_overlap(mol)
    Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
    u = (M_f - Skf) / HD
    W_rec = (np.eye(nbf) + HD * u.T) @ W0
    mol.data['OQP::VEC_MO_A'] = W_rec
    mol.data['OQP::VEC_MO_B'] = W_rec.copy()
    Ax_rec_d = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0

    lhs = (Ax_d - Ax0_d) / HD
    staged_d = (Ax_rec_d - Ax0_d) / HD
    print('===== v14 swap verdicts =====')
    print(f'IN-DISPLACED-STATE: |lhs - staged_d| = '
          f'{np.linalg.norm(lhs - staged_d):.6f} '
          f'maxdiff={np.abs(lhs - staged_d).max():.4e}')

    # restore to reference state; same staged pair
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    mol.data['OQP::VEC_MO_A'] = W_rec
    mol.data['OQP::VEC_MO_B'] = W_rec.copy()
    Ax_rec_r = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    base_r = matvec(XJ)
    staged_r = (Ax_rec_r - base_r) / HD
    print(f'REFERENCE-STATE  : |lhs - staged_r| = '
          f'{np.linalg.norm(lhs - staged_r):.6f} '
          f'maxdiff={np.abs(lhs - staged_r).max():.4e}')
    print(f'STATE SHIFT      : |staged_d - staged_r| = '
          f'{np.linalg.norm(staged_d - staged_r):.6f}')
    print(f'(|base_ref - base_r| sanity = {np.abs(base_ref - base_r).max():.2e})')

if __name__ == '__main__':
    main()
