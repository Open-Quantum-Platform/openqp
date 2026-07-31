"""Audit probe: check the sign/orientation conventions of numerical_nac().

Canonical convention under test:
    dcv[i,j] = d_ij = <Psi_i|d/dR Psi_j>   antisymmetric
    nacv[i,j] = h_ij = (E_j - E_i) * d_ij  symmetric

NOTE: numerical_nac() fans out with multiprocessing (spawn on macOS), so
everything must live under the __main__ guard.
"""
import sys


def main():
    import oqp                   # must precede numpy (ILP64 LAPACK interposition)
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    import numpy as np

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_probe.log'))
    r.run()
    mol = r.mol

    nac = NAC(mol)
    nacv, dcv, flags = nac.numerical_nac()
    ns = nac.nstate
    E = np.array(mol.energies[1:1 + ns])

    print('\n=== flags:', flags)
    print('=== state energies:', E)

    def asym(A):
        return np.linalg.norm(A + np.transpose(A, (1, 0, 2, 3)))

    def sym(A):
        return np.linalg.norm(A - np.transpose(A, (1, 0, 2, 3)))

    print(f'\n|dcv + dcv^T|  = {asym(dcv):.3e}   (0 => dcv ANTIsymmetric, canonical)')
    print(f'|dcv - dcv^T|  = {sym(dcv):.3e}')
    print(f'|nacv + nacv^T| = {asym(nacv):.3e}  (0 => nacv antisymmetric, NON-canonical)')
    print(f'|nacv - nacv^T| = {sym(nacv):.3e}   (0 => nacv symmetric, canonical)')

    for i in range(ns):
        for j in range(ns):
            if i >= j:
                continue
            gap_can = E[j] - E[i]
            h_can = gap_can * dcv[i, j]
            h_code = nacv[i, j]
            num = np.dot(h_code.ravel(), h_can.ravel())
            den = np.linalg.norm(h_code) * np.linalg.norm(h_can) + 1e-30
            print(f'\npair ({i+1},{j+1})  E_j-E_i = {gap_can:+.6f}')
            print(f'  |d_ij| = {np.linalg.norm(dcv[i,j]):.6f}')
            print(f'  code nacv[i,j] vs canonical (E_j-E_i)*d_ij : cos = {num/den:+.6f}, '
                  f'ratio = {np.linalg.norm(h_code)/(np.linalg.norm(h_can)+1e-30):.6f}')


if __name__ == '__main__':
    main()
