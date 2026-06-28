"""MRSF excited-state densities read from the OQP tagarray bridge.

The MRSF energy driver (``tdhf_mrsf_energy``) exposes three arrays (written by
the ``misc-excited-analysis`` patch):

* ``OQP::td_trans_density_mo`` -- ``trden(nbf,nbf,nstates*nstates)`` in the
  **alpha-MO basis**.  The off-diagonal blocks ``(ist != jst)`` are the
  state-interaction one-particle transition density matrices (1-TDM)
  :math:`\\gamma^{i\\to j}`; the diagonal blocks ``(ist == jst)`` are the
  *traceless difference* 1-RDMs :math:`\\Delta\\gamma^{n}=\\gamma^{n}-\\gamma^{\\rm ref}`.
  Crucially, because the MRSF ground state ``S0`` is itself a response root
  (not the reference), every block has **only occ-occ and vir-vir parts** --
  the occ-vir block is zero.  This is the genuine state-interaction structure,
  not a reference->amplitude object.
* ``OQP::td_trans_dipole`` -- ``dip(3,nstates,nstates)`` OQP's own transition
  dipoles (a.u.) at the center of mass.
* ``OQP::td_dip_ao`` -- AO electric-dipole integrals (packed lower-triangle).

This module is validated at GATE 2: reconstructing
:math:`\\mu^{i\\to j}=-\\mathrm{Tr}(\\gamma^{i\\to j}_{\\rm AO}\\,r)` from the exposed
1-TDM reproduces ``dip`` to ~1e-15, and ``Tr(gamma^n_AO . S) = N``.
"""
import numpy as np

__all__ = ["MRSFExcitedStates"]


def _f_order(mol, tag, shape):
    """Read a tagarray entry and re-interpret the raw (Fortran column-major)
    buffer with the given Fortran ``shape``."""
    raw = np.array(mol.data[tag], copy=True).ravel(order="C")
    return raw.reshape(shape, order="F")


def _unpack_lt(packed, n):
    """Unpack an OQP packed (lower-triangle) symmetric matrix to dense."""
    m = np.zeros((n, n))
    rows, cols = np.tril_indices(n)
    m[rows, cols] = packed
    m[cols, rows] = packed
    return m


class MRSFExcitedStates:
    """Excited-state density analysis for a finished MRSF calculation.

    Parameters
    ----------
    mol : oqp.molecule.molecule.Molecule
        A molecule whose MRSF energy run has completed in the *current*
        process (so the tagarray still holds the exposed densities).
    """

    def __init__(self, mol):
        self.mol = mol
        self.nbf = int(round(np.asarray(mol.data["OQP::VEC_MO_A"]).size ** 0.5))
        self.na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])   # noca
        self.nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])   # nocb
        self.n_elec = self.na + self.nb
        self.energies = np.asarray(mol.data["OQP::td_energies"]).ravel().copy()
        self.nstates = self.energies.size
        try:
            self.mult = int(mol.config.get("tdhf", {}).get("multiplicity", 1) or 1)
        except Exception:
            self.mult = 1

        nbf = self.nbf
        # AO <-> MO machinery (alpha MOs; ROHF shares them with beta).
        self.C = np.asarray(mol.data["OQP::VEC_MO_A"]).reshape(nbf, nbf).T  # C[ao, mo]
        self.S = _unpack_lt(np.asarray(mol.data["OQP::SM"]).ravel(), nbf)
        # AO electric-dipole integral matrices (symmetric), x/y/z.
        dip_packed = _f_order(mol, "OQP::td_dip_ao", (nbf * (nbf + 1) // 2, 3))
        self.R = [_unpack_lt(dip_packed[:, i], nbf) for i in range(3)]
        # OQP's own transition dipoles.
        self.dip_oqp = _f_order(mol, "OQP::td_trans_dipole", (3, self.nstates, self.nstates))
        # The full state-interaction density tensor (alpha-MO basis).
        self._T = _f_order(mol, "OQP::td_trans_density_mo",
                           (nbf, nbf, self.nstates, self.nstates))

        # Reference ROHF (high-spin) 1-RDM in the MO basis: doubly occupied
        # 0..nb-1, singly occupied nb..na-1 (SOMOs), empty above.
        self.gamma_ref_mo = np.zeros((nbf, nbf))
        for i in range(self.nb):
            self.gamma_ref_mo[i, i] = 2.0
        for i in range(self.nb, self.na):
            self.gamma_ref_mo[i, i] = 1.0

    # ---- transition densities (state interaction, alpha-MO basis) ----
    def tdm_mo(self, i, j):
        """1-TDM gamma^{i->j} in the alpha-MO basis (0-based state indices).

        Only the upper triangle is stored by OQP; gamma^{j->i} == (gamma^{i->j})^T."""
        if i <= j:
            return self._T[:, :, i, j].copy()
        return self._T[:, :, j, i].T.copy()

    def tdm_ao(self, i, j):
        """1-TDM gamma^{i->j} in the AO basis: C gamma_mo C^T."""
        return self.C @ self.tdm_mo(i, j) @ self.C.T

    def diff_density_mo(self, n):
        """Traceless difference density Delta gamma^n = gamma^n - gamma^ref (MO)."""
        return self._T[:, :, n, n].copy()

    def state_density_mo(self, n):
        """Unrelaxed state 1-RDM gamma^n in the MO basis (gamma_ref + Delta)."""
        return self.gamma_ref_mo + self._T[:, :, n, n]

    def state_density_ao(self, n):
        """Unrelaxed state 1-RDM gamma^n in the AO basis."""
        g = self.state_density_mo(n)
        return self.C @ g @ self.C.T

    def mo_density_ao(self, mo_index):
        """AO density of a single (alpha) MO |phi><phi| (for cube export)."""
        c = self.C[:, mo_index]
        return np.outer(c, c)

    # ---- dipoles ----
    def transition_dipole(self, i, j):
        """Reconstruct mu^{i->j} = -Tr(gamma_ao . r) (a.u.), x/y/z."""
        g_ao = self.tdm_ao(i, j)
        return np.array([-np.sum(g_ao * self.R[k]) for k in range(3)])

    def oscillator_strength(self, i, j):
        """f = (2/3) * dE * |mu|^2 (dE = |E_j - E_i| in a.u.)."""
        mu = self.transition_dipole(i, j)
        de = abs(self.energies[j] - self.energies[i])
        return (2.0 / 3.0) * de * float(mu @ mu)

    # ---- spin-flip excitation amplitudes (for NTOs) ----
    def amplitude_matrix(self, n):
        """Spin-adapted MRSF spin-flip amplitude matrix X^{(n)} (noca x nvirb).

        Rows = alpha-occupied MOs 0..na-1 (holes), columns = beta-virtual MOs
        nb..nbf-1 (particles).  Replicates the SOMO sqrt(2) handling of
        ``get_mrsf_transition_density`` so that ``trans_den(X,X)`` reproduces the
        exposed difference density (validated in the GATE 3 self-check)."""
        nbf, noca, nocb = self.nbf, self.na, self.nb
        nvirb = nbf - nocb
        bvec = _f_order(self.mol, "OQP::td_bvec_mo", (noca * nvirb, self.nstates))
        x = np.zeros((noca, nvirb))
        sqrt2 = np.sqrt(2.0)
        # 1-based special configuration indices (singlet mrst=1 / triplet mrst=3)
        ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
        ijlr2 = (noca - nocb - 1) * noca + noca
        ijg = (noca - 1 - nocb - 1) * noca + noca
        ijd = (noca - nocb - 1) * noca + noca - 1
        for i in range(1, noca + 1):
            for j in range(nocb + 1, nbf + 1):
                ij = (j - nocb - 1) * noca + i        # 1-based amplitude index
                col = j - nocb - 1                    # 0-based virtual index
                if self.mult == 1:
                    if ij == ijlr1:
                        x[i - 1, col] = bvec[ijlr1 - 1, n] / sqrt2
                        continue
                    if ij == ijlr2:
                        x[i - 1, col] = -bvec[ijlr1 - 1, n] / sqrt2
                        continue
                else:  # mrst == 3
                    if ij == ijlr1:
                        x[i - 1, col] = bvec[ijlr1 - 1, n] / sqrt2
                        continue
                    if ij == ijlr2:
                        x[i - 1, col] = bvec[ijlr1 - 1, n] / sqrt2
                        continue
                    if ij == ijg or ij == ijd:
                        x[i - 1, col] = 0.0
                        continue
                x[i - 1, col] = bvec[ij - 1, n]
        return x

    def trans_den_from_amplitudes(self, xi, xj):
        """Reproduce OQP's ``get_trans_den`` (alpha-MO basis) from amplitude
        matrices of two roots.  Used to validate ``amplitude_matrix``."""
        nbf, noca, nocb = self.nbf, self.na, self.nb
        nvirb = nbf - nocb
        sqrt2 = np.sqrt(2.0)
        trd = np.zeros((nbf, nbf))
        somo = (noca - 2, noca - 1)            # 0-based SOMO occ indices
        # vir/vir block
        for a in range(nvirb):
            for b in range(nvirb):
                tmp = 0.0
                for k in range(noca):
                    ioo = (a in (0, 1)) and (k in somo)
                    joo = (b in (0, 1)) and (k in somo)
                    scal = sqrt2 if (ioo ^ joo) else 1.0
                    tmp += scal * xj[k, a] * xi[k, b]
                trd[a + nocb, b + nocb] += tmp
        # occ/occ block
        for i in range(noca):
            for j in range(noca):
                tmp = 0.0
                for k in range(nvirb):
                    ioo = (i in somo) and (k in (0, 1))
                    joo = (j in somo) and (k in (0, 1))
                    scal = sqrt2 if (ioo ^ joo) else 1.0
                    tmp += scal * xj[j, k] * xi[i, k]
                trd[i, j] -= tmp
        return trd
