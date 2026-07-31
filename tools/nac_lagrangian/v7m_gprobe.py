"""v7m: G-build fidelity probe. (i) eps-linearity (eps vs eps/2);
(ii) self-adjointness Tr[Q1 G[Q2]] == Tr[Q2 G[Q1]]; (iii) G[0] == 0.
Run: python v7m_gprobe.py <input.inp>
"""
import os, sys
import numpy as np

def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7m.log'))
    r.run()
    mol = r.mol
    ctx = FK.build_context(mol)
    nbf = ctx['nbf']
    nbf_tri = nbf * (nbf + 1) // 2

    def unpack_sym(pk):
        M = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                M[p, q] = pk[idx]; M[q, p] = pk[idx]; idx += 1
        return M

    def pack_sym(M):
        pk = np.zeros(nbf_tri)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                pk[idx] = M[p, q]; idx += 1
        return pk

    KEYS = ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
            'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']
    saved = {k: np.array(mol.data[k], copy=True) for k in KEYS}
    F0a = unpack_sym(saved['OQP::FOCK_A'].ravel())
    F0b = unpack_sym(saved['OQP::FOCK_B'].ravel())
    D0a = unpack_sym(saved['OQP::DM_A'].ravel())
    D0b = unpack_sym(saved['OQP::DM_B'].ravel())
    for attr in ('maxit', 'scf_maxit'):
        try:
            setattr(mol.data._data.control, attr, 1)
            break
        except Exception:
            pass
    mol.config['scf']['maxit'] = 1

    def gbuild(Pa, Pb, eps):
        mol.data['OQP::DM_A'] = pack_sym(D0a + eps * Pa).reshape(saved['OQP::DM_A'].shape)
        mol.data['OQP::DM_B'] = pack_sym(D0b + eps * Pb).reshape(saved['OQP::DM_B'].shape)
        oqp.hf_energy(mol)
        Fa = unpack_sym(np.array(mol.data['OQP::FOCK_A'], copy=True).ravel())
        Fb = unpack_sym(np.array(mol.data['OQP::FOCK_B'], copy=True).ravel())
        for k in KEYS:
            mol.data[k] = saved[k].copy()
        return (Fa - F0a) / eps, (Fb - F0b) / eps

    rng = np.random.RandomState(7)
    Q1 = rng.randn(nbf, nbf); Q1 = 0.05 * (Q1 + Q1.T)
    Q2 = rng.randn(nbf, nbf); Q2 = 0.05 * (Q2 + Q2.T)

    Ga1, Gb1 = gbuild(Q1, Q1, 1e-4)
    Ga1h, Gb1h = gbuild(Q1, Q1, 5e-5)
    print(f'linearity: |G(eps)-G(eps/2)|max = {max(np.abs(Ga1-Ga1h).max(), np.abs(Gb1-Gb1h).max()):.3e} '
          f'(|G|max={np.abs(Ga1).max():.4f})')
    Ga2, Gb2 = gbuild(Q2, Q2, 1e-4)
    t12 = float(np.sum(Q1 * Ga2) + np.sum(Q1 * Gb2))
    t21 = float(np.sum(Q2 * Ga1) + np.sum(Q2 * Gb1))
    print(f'self-adjointness: Tr[Q1 G[Q2]]={t12:.6f} vs Tr[Q2 G[Q1]]={t21:.6f} '
          f'diff={abs(t12-t21):.3e}')
    Gz_a, Gz_b = gbuild(np.zeros((nbf, nbf)), np.zeros((nbf, nbf)), 1e-4)
    print(f'G[0]: max = {max(np.abs(Gz_a).max(), np.abs(Gz_b).max()):.3e}')

if __name__ == '__main__':
    main()
