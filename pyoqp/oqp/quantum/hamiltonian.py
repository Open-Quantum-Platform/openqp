"""Build the second-quantized molecular Hamiltonian from an OpenQP molecule.

This is the bridge that turns a converged OpenQP mean-field calculation into
the input expected by quantum-computing electronic-structure workflows
(Qiskit Nature, OpenFermion, PennyLane-via-OpenFermion, Block2, ...). The
output is either a :class:`MolecularHamiltonian` of NumPy tensors in the MO
basis, or a FCIDUMP file.

What OpenQP exposes today
-------------------------
The Python data container (``mol.data[...]``) provides everything needed for
the *one-electron* part of the Hamiltonian and all metadata:

* ``OQP::Hcore`` -- core Hamiltonian in the AO basis (packed triangular)
* ``OQP::SM``    -- overlap (packed triangular)
* ``OQP::VEC_MO_A`` / ``OQP::VEC_MO_B`` -- MO coefficients
* ``enuc`` -- nuclear repulsion energy
* ``nelec_A`` / ``nelec_B`` -- electron counts

Two-electron integrals
----------------------
The two-electron repulsion integrals (ERIs) are produced by the Fortran getter
``oqp.int2e(mol)``, which populates the ``OQP::ERI_AO`` tag with the full
``nbf**4`` AO tensor (chemist notation). :func:`from_openqp` calls it
automatically, so a converged calculation yields a complete FCIDUMP with no
external integral source. The ERIs come from the same engine and AO basis as
``OQP::Hcore`` and the MO coefficients, keeping the Hamiltonian consistent.

This is the conventional in-core path (memory ~ ``nbf**4``); it is meant for
small active systems / quantum-computing experiments, not production basis
sets. Callers may still override the source with ``eri_ao=`` or an
``eri_provider`` callable, or skip the two-body part with ``compute_eri=False``.
"""

from dataclasses import dataclass, field

import numpy as np

from oqp.quantum.integrals import (
    unpack_triangular,
    ao_to_mo_1body,
    ao_to_mo_2body,
)
from oqp.quantum.fcidump import write_fcidump


@dataclass
class MolecularHamiltonian:
    """Second-quantized electronic Hamiltonian in the MO basis.

    Attributes
    ----------
    one_body : numpy.ndarray, shape (norb, norb)
        ``h_pq`` one-electron integrals.
    two_body : numpy.ndarray or None, shape (norb, norb, norb, norb)
        ``(pq|rs)`` two-electron integrals (chemist notation). ``None`` when
        ERIs were unavailable (one-electron-only export).
    core_energy : float
        Scalar energy (nuclear repulsion + any frozen-core contribution).
    n_electrons : int
        Number of correlated electrons.
    ms2 : int
        ``2 * S_z`` = (n_alpha - n_beta).
    orbsym : list of int
        Per-orbital symmetry labels (defaults to all 1).
    """

    one_body: np.ndarray
    core_energy: float
    n_electrons: int
    two_body: np.ndarray = None
    ms2: int = 0
    orbsym: list = field(default_factory=list)

    @property
    def n_orbitals(self):
        return self.one_body.shape[0]

    def to_fcidump(self, filename, tol=1.0e-12):
        """Write this Hamiltonian to a FCIDUMP file.

        Requires two-electron integrals to be present.
        """
        if self.two_body is None:
            raise ValueError(
                "Cannot write a FCIDUMP without two-electron integrals. "
                "Pass eri_ao=/eri_provider= to from_openqp (see module docs "
                "for the planned oqp.int2e hook).")
        orbsym = self.orbsym or [1] * self.n_orbitals
        write_fcidump(
            filename, self.one_body, self.two_body, self.core_energy,
            self.n_electrons, ms2=self.ms2, orbsym=orbsym, tol=tol)
        return filename


def _get_nbf(mol):
    basis = mol.data.get_basis()
    if not basis:
        raise RuntimeError("No basis available; run/apply a basis first.")
    return int(basis["nbf"])


def from_openqp(mol, eri_ao=None, eri_provider=None, mo_coeff=None,
                spin="alpha", compute_eri=True):
    """Construct a :class:`MolecularHamiltonian` from an OpenQP ``Molecule``.

    Parameters
    ----------
    mol : oqp.molecule.molecule.Molecule
        A molecule whose SCF has completed (MO coefficients populated).
    eri_ao : array_like, optional
        Two-electron AO integrals ``(mu nu|la si)`` in chemist notation,
        shape ``(nao, nao, nao, nao)``. If given, the two-body MO tensor is
        built and a full FCIDUMP can be written.
    eri_provider : callable, optional
        ``eri_provider(mol) -> eri_ao``; used when ``eri_ao`` is not supplied.
        Lets callers plug in OpenQP's native ERIs (once exposed) or an
        external engine without changing this code.
    mo_coeff : array_like, optional
        Override the MO coefficient matrix (AO rows, MO columns). Defaults to
        the converged restricted/alpha MOs (``OQP::VEC_MO_A``).
    spin : {"alpha", "beta"}
        Which set of converged MOs to use when ``mo_coeff`` is not given.
    compute_eri : bool
        When no ``eri_ao``/``eri_provider`` is given, compute OpenQP's native
        AO ERIs via ``oqp.int2e`` (default). Set ``False`` to build a
        one-electron-only Hamiltonian (no two-body tensor, no FCIDUMP).

    Returns
    -------
    MolecularHamiltonian
    """
    nbf = _get_nbf(mol)

    # --- one-electron (always available) ---------------------------------
    hcore_ao = unpack_triangular(mol.data["OQP::Hcore"], nbf)
    if mo_coeff is None:
        tag = "OQP::VEC_MO_A" if spin == "alpha" else "OQP::VEC_MO_B"
        mo_coeff = mol._mo_coefficients(tag, nbf)
    mo_coeff = np.asarray(mo_coeff, dtype=float)

    h1_mo = ao_to_mo_1body(hcore_ao, mo_coeff)

    # --- metadata ---------------------------------------------------------
    na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
    nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])
    n_electrons = na + nb
    ms2 = na - nb

    # Nuclear-repulsion energy. OpenQP stores it under `vnn` (== `nenergy`);
    # `enuc` exists in the struct but is left at zero, so prefer vnn.
    ecore = 0.0
    for field in ("vnn", "nenergy", "enuc"):
        try:
            val = float(mol.data[field])
        except Exception:
            continue
        if val != 0.0:
            ecore = val
            break

    # --- two-electron -----------------------------------------------------
    # Resolution order: explicit eri_ao -> custom provider -> OpenQP's native
    # AO integrals via oqp.int2e (OQP::ERI_AO). The native path keeps the ERIs
    # consistent with OQP::Hcore and the MO coefficients (same AO ordering and
    # normalization), which is what makes the resulting FCIDUMP correct.
    h2_mo = None
    if eri_ao is None and eri_provider is not None:
        eri_ao = eri_provider(mol)
    if eri_ao is None and compute_eri:
        try:
            from oqp.library.ints_2e import eri_ao as _native_eri_ao
            eri_ao = _native_eri_ao(mol)
        except Exception as exc:  # pragma: no cover - requires compiled oqp
            raise RuntimeError(
                "Could not obtain two-electron integrals from OpenQP "
                "(oqp.int2e / OQP::ERI_AO). Pass eri_ao= or eri_provider=, or "
                "set compute_eri=False for a one-electron-only Hamiltonian. "
                f"Underlying error: {exc}")
    if eri_ao is not None:
        h2_mo = ao_to_mo_2body(eri_ao, mo_coeff)

    return MolecularHamiltonian(
        one_body=h1_mo,
        two_body=h2_mo,
        core_energy=ecore,
        n_electrons=n_electrons,
        ms2=ms2,
        orbsym=[1] * h1_mo.shape[0],
    )
