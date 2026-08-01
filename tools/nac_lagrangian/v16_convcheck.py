"""v16 (THE CONVENTION CHECK): reconstruct M and Sk from the run's own
overlap_ao record and the coefficient records; find where
M - Sk != C0^T S_cross (C'sg - C0). The failing entries name the
convention behind V.
Run: python v16_convcheck.py <inp> <coord>
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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v16.log'))
    r.run()
    mol = r.mol
    ctx = FK.build_context(mol)
    nbf = ctx['nbf']
    noca, nocb = ctx['noca'], ctx['nocb']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C0 = W0.T           # C0[ao, mo]
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
    Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    oqp.get_structures_ao_overlap(mol)
    M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    S_np = np.array(mol.data['OQP::overlap_ao_non_orthogonal'], copy=True)
    M_f = M_np.reshape(-1).reshape((nbf, nbf)).T     # Fortran matrix
    S_f = S_np.reshape(-1).reshape((nbf, nbf)).T     # Fortran matrix
    # Sk call
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    oqp.get_structures_ao_overlap(mol)
    Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    Sk_f = Sk_np.reshape(-1).reshape((nbf, nbf)).T

    Cd = Wd.T          # displaced C'[ao, mo] (unfixed signs)
    # candidate reconstructions of M from S_f and coefficients
    cands = {
        'C0T S Cd':  C0.T @ S_f @ Cd,
        'C0T ST Cd': C0.T @ S_f.T @ Cd,
        'CdT S C0':  Cd.T @ S_f @ C0,
        'CdT ST C0': Cd.T @ S_f.T @ C0,
    }
    print('===== which reconstruction matches the M record? =====')
    for k, v in cands.items():
        print(f'  {k}: maxdiff={np.abs(v - M_f).max():.3e}')
    print('===== and the Sk record? =====')
    for k, v in (('C0T S C0', C0.T @ S_f @ C0), ('C0T ST C0', C0.T @ S_f.T @ C0)):
        print(f'  {k}: maxdiff={np.abs(v - Sk_f).max():.3e}')

    # the identity check with the established orientation (try both)
    dC = Cd - C0
    for k, S in (('S', S_f), ('ST', S_f.T)):
        lhsA = (M_f - Sk_f)
        # NOTE M may be built with Cd (uncorrected signs); compare without sg
        rhsA = C0.T @ S @ dC
        print(f'identity [{k}]: max|M - Sk - C0T S dC| = '
              f'{np.abs(lhsA - rhsA).max():.3e}')
    # locate worst violations for the better orientation
    best = min((('S', S_f), ('ST', S_f.T)),
               key=lambda kv: np.abs((M_f - Sk_f) - C0.T @ kv[1] @ dC).max())
    S = best[1]
    D = (M_f - Sk_f) - C0.T @ S @ dC
    idx = np.argsort(-np.abs(D).ravel())[:10]
    print(f'worst identity violations [{best[0]}]:')
    for f in idx:
        p, q = divmod(int(f), nbf)
        print(f'  ({p},{q}): D={D[p,q]:+.3e}  '
              f'(p:{"doc" if p<nocb else ("socc" if p<noca else "virt")}, '
              f'q:{"doc" if q<nocb else ("socc" if q<noca else "virt")})')

if __name__ == '__main__':
    main()
