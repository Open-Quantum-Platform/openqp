"""Integral utilities for building the second-quantized molecular Hamiltonian.

These helpers are deliberately free of any dependency on the compiled ``oqp``
extension so that they can be unit-tested and reused independently of an
OpenQP build. They operate on plain NumPy arrays.

Conventions
-----------
* One-electron integrals ``h_pq`` are stored as a square ``(norb, norb)``
  matrix in the (real) molecular-orbital basis.
* Two-electron integrals use **chemist notation** ``(pq|rs)`` -- the same
  convention used by the FCIDUMP interchange format, PySCF and most quantum
  chemistry codes. The stored tensor has shape ``(norb, norb, norb, norb)``
  with index order ``[p, q, r, s]`` meaning
  ``(pq|rs) = \\int \\phi_p(1)\\phi_q(1) (1/r12) \\phi_r(2)\\phi_s(2) dr1 dr2``.
"""

import numpy as np


def unpack_triangular(packed, norb, lower=True):
    """Expand a packed triangular matrix into a dense symmetric matrix.

    OpenQP stores symmetric one-electron matrices (``OQP::Hcore``, ``OQP::SM``,
    Fock, density, ...) as the lower triangle in column-major Fortran order,
    which -- read back as a flat C array -- is exactly the row-major lower
    triangle. This mirrors the unpacking already done in ``Molecule``.

    Parameters
    ----------
    packed : array_like
        Flat array of length ``norb * (norb + 1) // 2``.
    norb : int
        Dimension of the resulting square matrix.
    lower : bool
        Whether ``packed`` holds the lower (default) or upper triangle.

    Returns
    -------
    numpy.ndarray
        Dense ``(norb, norb)`` symmetric matrix.
    """
    packed = np.asarray(packed, dtype=float).ravel()
    expected = norb * (norb + 1) // 2
    if packed.size != expected:
        raise ValueError(
            f"packed triangular array has {packed.size} elements, "
            f"expected {expected} for norb={norb}")
    mat = np.zeros((norb, norb), dtype=float)
    if lower:
        rows, cols = np.tril_indices(norb)
    else:
        rows, cols = np.triu_indices(norb)
    mat[rows, cols] = packed
    mat[cols, rows] = packed
    return mat


def ao_to_mo_1body(h_ao, mo_coeff):
    """Transform a one-electron AO matrix into the MO basis.

    ``h_pq = sum_{mu,nu} C_{mu p} h_{mu nu} C_{nu q}``.

    Parameters
    ----------
    h_ao : array_like, shape (nao, nao)
        One-electron integrals in the AO basis (e.g. core Hamiltonian).
    mo_coeff : array_like, shape (nao, nmo)
        MO coefficient matrix (AO rows, MO columns).
    """
    h_ao = np.asarray(h_ao, dtype=float)
    c = np.asarray(mo_coeff, dtype=float)
    return c.T @ h_ao @ c


def ao_to_mo_2body(eri_ao, mo_coeff):
    """Transform two-electron AO integrals ``(mu nu|la si)`` into the MO basis.

    Performs the standard O(N^5) four-index transformation, contracting one
    index at a time, and returns ``(pq|rs)`` in chemist notation.

    Parameters
    ----------
    eri_ao : array_like, shape (nao, nao, nao, nao)
        Two-electron repulsion integrals in the AO basis, chemist notation.
    mo_coeff : array_like, shape (nao, nmo)
        MO coefficient matrix (AO rows, MO columns).

    Returns
    -------
    numpy.ndarray, shape (nmo, nmo, nmo, nmo)
        ``(pq|rs)`` in the MO basis.
    """
    eri = np.asarray(eri_ao, dtype=float)
    c = np.asarray(mo_coeff, dtype=float)
    if eri.ndim != 4:
        raise ValueError("eri_ao must be a rank-4 tensor (nao,nao,nao,nao)")
    # `c` is (nao, nmo): rows are AOs, columns are MOs (same convention as
    # ao_to_mo_1body). Contract the AO axis of `c` with each AO axis of the
    # integral tensor, one index at a time, to keep the cost at O(N^5).
    tmp = np.einsum('ip,ijkl->pjkl', c, eri, optimize=True)
    tmp = np.einsum('jq,pjkl->pqkl', c, tmp, optimize=True)
    tmp = np.einsum('kr,pqkl->pqrl', c, tmp, optimize=True)
    mo = np.einsum('ls,pqrl->pqrs', c, tmp, optimize=True)
    return mo
