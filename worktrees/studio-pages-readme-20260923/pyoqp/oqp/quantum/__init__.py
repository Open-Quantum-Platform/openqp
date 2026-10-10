"""Quantum-computing bridge for OpenQP.

Export an OpenQP mean-field calculation as a second-quantized molecular
Hamiltonian / FCIDUMP, ready for quantum-computing electronic-structure
toolkits (Qiskit Nature, OpenFermion, Block2, ...).

Typical use::

    from oqp.quantum import from_openqp

    ham = from_openqp(mol, eri_ao=eri)   # eri = AO two-electron integrals
    ham.to_fcidump("molecule.FCIDUMP")

The integral-transform and FCIDUMP I/O helpers
(:mod:`oqp.quantum.integrals`, :mod:`oqp.quantum.fcidump`) have no dependency
on the compiled ``oqp`` extension and can be used stand-alone.
"""

from oqp.quantum.hamiltonian import MolecularHamiltonian, from_openqp
from oqp.quantum.fcidump import write_fcidump, read_fcidump
from oqp.quantum.integrals import (
    unpack_triangular,
    ao_to_mo_1body,
    ao_to_mo_2body,
)

__all__ = [
    "MolecularHamiltonian",
    "from_openqp",
    "write_fcidump",
    "read_fcidump",
    "unpack_triangular",
    "ao_to_mo_1body",
    "ao_to_mo_2body",
]
