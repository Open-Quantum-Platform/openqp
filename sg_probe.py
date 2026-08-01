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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_sg.log'))
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
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False
    HD = 1e-3
    for c0 in (0, 5):
        for sgn, lab in ((1, '+h'), (-1, '-h')):
            cp = xyz0.copy(); cp[c0] += sgn * HD
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
            neg = [i for i in range(nbf) if sg[i] < 0]
            socc = list(range(nocb, noca))
            print(f'c0={c0} {lab}: flipped orbitals (0-based) = {neg}  '
                  f'socc={socc}  socc-flips={[i for i in neg if i in socc]}',
                  flush=True)

if __name__ == '__main__':
    main()
