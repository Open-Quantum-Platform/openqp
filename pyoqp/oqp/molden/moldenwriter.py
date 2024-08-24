"""
Molden compatibility module
writer class
"""

import numpy as np

from oqp.periodic_table import ELEMENTS_NAME, SYMBOL_MAP


class MoldenWriter:
    """
     Class to write molecule information in Molden format

     Molecular orbitals sequence for various angular momentum:
     5D: D 0, D+1, D-1, D+2, D-2
     6D: xx, yy, zz, xy, xz, yz

     7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
     10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz

     9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
     15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
          xxyy xxzz yyzz xxyz yyxz zzxy
    """

    ORDERS = {0: [0],
              1: [0, 1, 2],
              2: [0, 1, 2, 3, 4, 5],
              3: [0, 1, 2, 5, 3, 4, 7, 8, 6, 9],
              4: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
              }
    SHELL_TYPES = 'xspdfg'
    NORMS = np.sqrt(np.pi * np.sqrt(np.pi) * np.array([1.0, 0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875]))

    NBFS = list(int((x + 1) * (x + 2) / 2) for x in range(5))

    def __init__(self, file):
        self.file = file
        print('[Molden Format]', file=self.file)

    def write_atoms(self, num_atoms, q_n, xyz, header=True, angstrom=False):
        """Write atomic coordinates section"""
        if header:
            units = 'Ang' if angstrom else 'AU'
            print(f'[Atoms] {units}', file=self.file)

        for i in range(num_atoms):
            atomic_number = int(q_n[i])
            name = ELEMENTS_NAME[atomic_number]
            x_c, y_c, z_c = xyz[i]
            print(f'{name} {i + 1} {atomic_number} {x_c:12.8f} {y_c:12.8f} {z_c:12.8f}',
                  file=self.file)

    def write_basis(self, num_atoms, basis, header=True):
        """Write GTO basis set section"""

        if header:
            print('[GTO]', file=self.file)

        molden_bas = tuple([] for _ in range(num_atoms))
        id_prim0 = 0
        for (sh_at, sh_typ, sh_nc) in zip(basis['centers'], basis['types'], basis['ncontr']):
            molden_bas[sh_at].append({
                'typ': MoldenWriter.SHELL_TYPES[sh_typ],
                'ang': sh_typ - 1,
                'nc': sh_nc,
                'alpha': (basis['alpha'])[id_prim0:id_prim0 + sh_nc],
                'coef': (basis['coef'])[id_prim0:id_prim0 + sh_nc],
            })
            id_prim0 += sh_nc

        for (i, atom) in enumerate(molden_bas):
            print(i + 1, file=self.file)
            for shell in atom:
                print(f'{shell["typ"]} {shell["nc"]}', file=self.file)
                ang = shell['ang']
                norm = MoldenWriter.NORMS[ang]
                for (alpha, coef) in zip(shell['alpha'], shell['coef']):
                    coef *= norm * (2 * alpha) ** -(ang / 2.0 + 0.75)
                    print(f'{alpha:22.12e} {coef:22.12e}', file=self.file)
            print("", file=self.file)

    def write_mo(self, basis, orbitals, energies, occupancies, spin, header=True):
        """Write molecular orbitals section"""

        if header:
            print('[MO]', file=self.file)

        reorder = []
        cur_bf = 0
        for i in range(basis['nsh']):
            ang = basis['types'][i] - 1
            neword = list(i + cur_bf for i in MoldenWriter.ORDERS[ang])
            reorder += neword
            cur_bf += MoldenWriter.NBFS[ang]

        for (eorb, orb, occup) in zip(energies, orbitals, occupancies):
            print(f'Ene= {eorb:15.8f}', file=self.file)
            print(f'Spin= {spin}', file=self.file)
            print(f'Occup= {occup:15.8f}', file=self.file)
            for (i, coef) in enumerate(orb[reorder]):
                print(f'{i + 1} {coef:15.8f}', file=self.file)


def write_frequency(mol, freqs, modes):
    atoms = mol.get_atoms()
    natom = len(atoms)
    xyz = mol.get_system().reshape((natom, 3))
    nmode = len(modes)
    modes = modes.reshape((nmode, natom, 3))
    freqs = freqs.reshape((-1, 1))

    frequency = '\n'.join(['%10.2f' % x for x in freqs])
    coord = ''
    for n, c in enumerate(xyz):
        x, y, z = c
        s = atoms[n]
        coord += '%-5s %24.8f %24.8f %24.8f\n' % (ELEMENTS_NAME[SYMBOL_MAP[int(s)]], x, y, z)

    vibs = ''
    for i in range(nmode):
        vibs += '  vibration  %5s\n%s\n' % (
            i + 1, '\n'.join([' '.join(['%24.8f' % y for y in x]) for x in modes[i]]))

    molden = """ [MOLDEN FORMAT]
 [N_FREQ]
%s
 [NATOM]
%s
 [FREQ]
%s
 [FR-COORD]
%s
 [FR-NORM-COORD]
%s

""" % (nmode, natom, frequency, coord, vibs)

    return molden
