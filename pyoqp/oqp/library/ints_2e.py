import numpy as np

import oqp

__all__ = ["ints_2e", "eri_ao"]


def ints_2e(mol):
    """Compute the two-electron repulsion integrals (ERIs) in the AO basis.

    Populates the ``OQP::ERI_AO`` tag with the full ``nbf**4`` tensor of
    integrals ``(mu nu|la si)`` in chemist notation. This is the conventional
    in-core path (memory grows as ``nbf**4``); it targets small active systems
    for FCIDUMP / quantum-computing export rather than production SCF.
    """
    oqp.int2e(mol)


def eri_ao(mol):
    """Return the AO two-electron integrals as a dense ``(nbf,nbf,nbf,nbf)``
    NumPy array in chemist notation, computing them if necessary.

    The full 8-fold permutational symmetry of the integrals makes the
    C/Fortran storage-order difference irrelevant, so the flat tag can be
    reshaped directly to ``(pq|rs)``.
    """
    basis = mol.data.get_basis()
    if not basis:
        raise RuntimeError("No basis available; apply a basis first.")
    nbf = int(basis["nbf"])

    flat = None
    try:
        cached = np.asarray(mol.data["OQP::ERI_AO"], dtype=float).ravel()
        if cached.size == nbf ** 4:
            flat = cached
    except Exception:
        flat = None

    if flat is None:
        ints_2e(mol)
        flat = np.asarray(mol.data["OQP::ERI_AO"], dtype=float).ravel()

    return flat.reshape(nbf, nbf, nbf, nbf)
