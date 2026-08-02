"""
Molden compatibility module
writer class
"""

from io import StringIO

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
    # Pure spherical-harmonic reorder: OpenQP internal order is m = -l..+l
    # (CCA); Molden wants m = 0, +1, -1, +2, -2, +3, -3, +4, -4. Entries are
    # zero-based positions into the internal (m=-l..+l) vector.
    SPH_ORDERS = {0: [0],
                  1: [0, 1, 2],
                  2: [2, 3, 1, 4, 0],
                  3: [3, 4, 2, 5, 1, 6, 0],
                  4: [4, 5, 3, 6, 2, 7, 1, 8, 0],
                  }
    SHELL_TYPES = 'spdfg'
    NORMS = np.sqrt(np.pi * np.sqrt(np.pi) * np.array([1.0, 0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875]))

    NBFS = list(int((x + 1) * (x + 2) / 2) for x in range(5))
    SPH_NBFS = list(2 * x + 1 for x in range(5))

    @staticmethod
    def _is_spherical(basis):
        """Infer whether the basis uses pure spherical shells, from the total
        basis-function count (auto-selection is uniform per basis set)."""
        angs = basis['angs']
        ncart = int(sum(int((a + 1) * (a + 2) // 2) for a in angs))
        nsph = int(sum(2 * int(a) + 1 for a in angs))
        nbf = int(basis['nbf'])
        if nbf == nsph and nsph != ncart:
            return True
        return False

    @staticmethod
    def supports_portable_ordering(basis):
        """Return whether this writer can represent every shell safely."""
        orders = MoldenWriter.SPH_ORDERS if MoldenWriter._is_spherical(basis) else MoldenWriter.ORDERS
        return all(
            int(angular_momentum) in orders
            and int(angular_momentum) < len(MoldenWriter.SHELL_TYPES)
            and int(angular_momentum) < len(MoldenWriter.NORMS)
            for angular_momentum in basis['angs']
        )

    def write_spherical_markers(self, basis):
        """Emit Molden [5D]/[7F]/[9G] markers when the basis is spherical."""
        if not MoldenWriter._is_spherical(basis):
            return
        present = set(int(a) for a in basis['angs'])
        if 2 in present:
            print('[5D]', file=self.file)
        if 3 in present:
            print('[7F]', file=self.file)
        if 4 in present:
            print('[9G]', file=self.file)

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
        for (sh_at, sh_ang, sh_nc) in zip(basis['centers'], basis['angs'], basis['ncontr']):
            molden_bas[sh_at].append({
                'typ': MoldenWriter.SHELL_TYPES[sh_ang],
                'ang': sh_ang,
                'nc': sh_nc,
                'alpha': (basis['alpha'])[id_prim0:id_prim0 + sh_nc],
                'coef': (basis['coef'])[id_prim0:id_prim0 + sh_nc],
            })
            id_prim0 += sh_nc

        for (i, atom) in enumerate(molden_bas):
            # Keep the legacy placeholders documented by the Molden format.
            # They are optional to Molden itself, but stricter third-party
            # readers (including GUI importers) expect both fields.
            print(f'{i + 1} 0', file=self.file)
            for shell in atom:
                print(f'{shell["typ"]} {shell["nc"]} 1.00', file=self.file)
                ang = shell['ang']
                for (alpha, coef) in zip(shell['alpha'], shell['coef']):
                    coef = MoldenWriter.molden_primitive_coefficient(ang, alpha, coef)
                    print(f'{alpha:22.12e} {coef:22.12e}', file=self.file)
            print("", file=self.file)

    @staticmethod
    def molden_primitive_coefficient(angular_momentum, exponent, coefficient):
        """Return a primitive coefficient in the convention written to Molden."""
        norm = MoldenWriter.NORMS[int(angular_momentum)]
        return float(coefficient) * norm * (2 * float(exponent)) ** -(float(angular_momentum) / 2.0 + 0.75)

    @staticmethod
    def orbital_reorder(basis):
        """Map OpenQP AO ordering to the ordering required by Molden."""
        spherical = MoldenWriter._is_spherical(basis)
        orders = MoldenWriter.SPH_ORDERS if spherical else MoldenWriter.ORDERS
        nbfs = MoldenWriter.SPH_NBFS if spherical else MoldenWriter.NBFS

        reorder = []
        cur_bf = 0
        for i in range(basis['nsh']):
            ang = int(basis['angs'][i])
            reorder.extend(j + cur_bf for j in orders[ang])
            cur_bf += nbfs[ang]
        return reorder

    def write_mo(self, basis, orbitals, energies, occupancies, spin, header=True,
                 symmetries=None, already_reordered=False):
        """Write molecular orbitals section"""

        if header:
            print('[MO]', file=self.file)

        reorder = MoldenWriter.orbital_reorder(basis)

        if symmetries is None:
            # C1 is the portable fallback when no point-group labels were
            # requested.  Emit a Sym record for every MO so strict readers do
            # not see a ragged orbital metadata table when Dyson states follow.
            symmetries = ('A' for _ in energies)

        for (eorb, orb, occup, symmetry) in zip(energies, orbitals, occupancies, symmetries):
            if symmetry:
                print(f'Sym= {symmetry}', file=self.file)
            print(f'Ene= {eorb:15.8f}', file=self.file)
            print(f'Spin= {spin}', file=self.file)
            print(f'Occup= {occup:15.8f}', file=self.file)
            coefficients = np.asarray(orb) if already_reordered else np.asarray(orb)[reorder]
            for (i, coef) in enumerate(coefficients):
                print(f'{i + 1} {coef:15.8f}', file=self.file)

    def write_frequency(self, mol, freqs, modes):
        """Append Molden frequency and normal-coordinate sections."""
        atoms = np.asarray(mol.get_atoms()).reshape(-1)
        natom = len(atoms)
        xyz = np.asarray(mol.get_system()).reshape((natom, 3))
        freqs = np.asarray(freqs).reshape(-1)
        modes = np.asarray(modes).reshape((len(freqs), natom, 3))

        print('[FREQ]', file=self.file)
        for frequency in freqs:
            print(f'{frequency:10.2f}', file=self.file)

        infrared = np.asarray(
            getattr(mol, 'infrared_intensities', np.zeros(0)), dtype=float).reshape(-1)
        raman = np.asarray(
            getattr(mol, 'raman_activities', np.zeros(0)), dtype=float).reshape(-1)
        if len(infrared) >= len(freqs):
            print('[INT]', file=self.file)
            for index in range(len(freqs)):
                print(f'{infrared[index]:16.8f}', file=self.file)
        if len(raman) >= len(freqs):
            # Raman activities are not part of the standard one-value-per-mode
            # Molden [INT] section.  Keep them in a separate extension so
            # strict readers can consume [INT] without shifting mode indices.
            print('[RAMAN]', file=self.file)
            for index in range(len(freqs)):
                print(f'{raman[index]:16.8f}', file=self.file)

        print('[FR-COORD]', file=self.file)
        for atom, coord in zip(atoms, xyz):
            symbol = ELEMENTS_NAME[SYMBOL_MAP[int(atom)]]
            print(f'{symbol:<5s} {coord[0]:24.8f} {coord[1]:24.8f} {coord[2]:24.8f}', file=self.file)

        print('[FR-NORM-COORD]', file=self.file)
        for index, vectors in enumerate(modes, start=1):
            print(f'  vibration  {index:5d}', file=self.file)
            for vector in vectors:
                print(' '.join(f'{component:24.8f}' for component in vector), file=self.file)


def write_frequency(mol, freqs, modes):
    output = StringIO()
    writer = MoldenWriter(output)
    writer.write_frequency(mol, freqs, modes)
    return output.getvalue()
