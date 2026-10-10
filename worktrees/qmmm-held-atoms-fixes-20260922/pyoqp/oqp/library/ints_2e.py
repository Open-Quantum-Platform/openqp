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


def eri_ao(mol, recompute=True):
    """Return dense AO two-electron integrals in chemist notation.

    By default this recomputes ``OQP::ERI_AO`` on every call.  A same-size tag
    may belong to an earlier geometry/calculation on the same ``Molecule`` and
    would silently mix stale ERIs with current one-electron data.  Set
    ``recompute=False`` only when the caller has explicitly managed cache
    validity.
    """
    basis = mol.data.get_basis()
    if not basis:
        raise RuntimeError("No basis available; apply a basis first.")
    nbf = int(basis["nbf"])

    if recompute:
        ints_2e(mol)
    else:
        try:
            cached = np.asarray(mol.data["OQP::ERI_AO"], dtype=float).ravel()
            if cached.size == nbf ** 4:
                return cached.reshape(nbf, nbf, nbf, nbf)
        except Exception:
            pass
        ints_2e(mol)

    flat = np.asarray(mol.data["OQP::ERI_AO"], dtype=float).ravel()
    if flat.size != nbf ** 4:
        raise RuntimeError(
            f"OQP::ERI_AO has {flat.size} elements, expected {nbf ** 4} "
            f"for nbf={nbf}.")
    return flat.reshape(nbf, nbf, nbf, nbf)
