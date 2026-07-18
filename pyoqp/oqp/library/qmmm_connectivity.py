"""QM/MM boundary (dangling-bond) connectivity across chemical bonds.

When a QM/MM partition cuts one or more covalent bonds, each severed QM–MM
bond leaves a *dangling bond* on the QM frontier atom.  The standard remedy is
to cap the dangling bond with a hydrogen **link atom** placed along the broken
QM–MM bond (the IMOMM / scaled-position scheme):

    R_L = R_QM + g * (R_MM - R_QM),      g = (r_H + r_QM) / (r_QM + r_MM)

where ``r_X`` are covalent radii.  Because the link-atom position is a *fixed
linear function* of its two real host atoms, the energy gradient computed for
the link atom is redistributed onto them by the chain rule:

    dE/dR_QM += (1 - g) * dE/dR_L
    dE/dR_MM +=       g * dE/dR_L

so no extra (non-physical) degrees of freedom are introduced.

This module is deliberately free of any OpenMM / OpenQP dependency: it operates
on plain bond pairs, atomic numbers and coordinate arrays so it can be unit
tested in isolation and reused by every QM/MM code path.  The covalent-radii
table matches ``oqp.utils.qmmm`` so the two QM/MM paths place link atoms
identically.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Covalent radii indexed by atomic number (Z); index 0 is a placeholder so the
# list can be indexed directly by Z.  Same values as ``oqp.utils.qmmm``.
COVALENT_RADII = [
    0.0, 0.354, 0.849, 1.336, 1.010, 0.838,
    0.757, 0.700, 0.658, 0.668, 0.920, 1.539, 1.421,
    1.244, 1.117, 1.101, 1.064, 1.044, 1.032, 1.953,
    1.761, 1.513, 1.412, 1.402, 1.345, 1.382, 1.270,
    1.241, 1.164, 1.302, 1.193, 1.260, 1.197, 1.211,
    1.190, 1.192, 1.147, 2.260, 2.052, 1.698, 1.564,
    1.473, 1.467, 1.322, 1.478, 1.332, 1.338, 1.386,
    1.403, 1.459, 1.398, 1.407, 1.386, 1.382, 1.267,
    2.570, 2.277, 1.943, 1.841, 1.823, 1.816, 1.801,
    1.780, 1.771, 1.735, 1.732, 1.710, 1.696, 1.673,
    1.660, 1.637, 1.671, 1.611, 1.511, 1.392, 1.372,
    1.372, 1.371, 1.364, 1.262, 1.340, 1.518, 1.459,
    1.512, 1.500, 1.545, 1.420, 2.880, 2.512, 1.983,
    1.721, 1.711, 1.684, 1.666, 1.657, 1.660, 1.801,
    1.761, 1.750, 1.724, 1.712, 1.689, 1.679, 1.698,
    1.850,
]


def covalent_radius(z: int) -> float:
    """Covalent radius (arb. units, consistent within the table) for atomic
    number ``z``.  Raises ``ValueError`` for elements outside the table."""
    if z < 1 or z >= len(COVALENT_RADII):
        raise ValueError(f"No covalent radius tabulated for atomic number {z}")
    return COVALENT_RADII[z]


def link_g_factor(qm_z: int, mm_z: int) -> float:
    """Scaled-position factor ``g`` for a hydrogen link atom capping a
    QM(``qm_z``)–MM(``mm_z``) bond.

    ``g`` is the fraction of the QM→MM vector at which the capping hydrogen is
    placed.  A physically sensible partition gives ``0 < g < 1``.
    """
    r_h = COVALENT_RADII[1]
    r_qm = covalent_radius(qm_z)
    r_mm = covalent_radius(mm_z)
    denom = r_qm + r_mm
    if denom <= 0.0:
        raise ValueError("Non-positive covalent-radius sum for link-atom scaling")
    return (r_h + r_qm) / denom


@dataclass(frozen=True)
class LinkAtom:
    """A hydrogen link atom capping one severed QM–MM bond.

    Attributes
    ----------
    qm_index : int
        Absolute (topology) index of the QM frontier atom.
    mm_index : int
        Absolute (topology) index of the MM host atom.
    host_row : int
        Row of the QM frontier atom within the ordered QM-atom list
        (``qm_atoms``), used to redistribute the link-atom gradient.
    g : float
        Scaled-position factor (see :func:`link_g_factor`).
    """

    qm_index: int
    mm_index: int
    host_row: int
    g: float


def detect_link_atoms(bonds, qm_atoms, z_of):
    """Detect dangling bonds and build the list of hydrogen link atoms.

    Parameters
    ----------
    bonds : iterable of (int, int)
        Covalent bonds as pairs of absolute atom indices (order irrelevant).
    qm_atoms : sequence of int
        Absolute indices of the QM-region atoms, in the order they are handed
        to the QM engine.
    z_of : callable(int) -> int
        Maps an absolute atom index to its atomic number.

    Returns
    -------
    list[LinkAtom]
        One entry per bond that has exactly one endpoint in the QM region,
        ordered by (qm_index, mm_index) for determinism.

    Raises
    ------
    ValueError
        If a single QM frontier atom is bonded to more than one MM atom
        *through the same partition in a way that is ambiguous* — actually
        multiple link atoms on one QM atom are allowed; the check here only
        rejects a non-positive/degenerate ``g`` (bad QM/MM partition).
    """
    qm_set = set(int(i) for i in qm_atoms)
    row_of = {int(a): r for r, a in enumerate(qm_atoms)}

    seen = set()
    links = []
    for i, j in bonds:
        i, j = int(i), int(j)
        i_qm, j_qm = i in qm_set, j in qm_set
        if i_qm == j_qm:
            # both QM (internal bond) or both MM (irrelevant): no dangling bond
            continue
        qm_index, mm_index = (i, j) if i_qm else (j, i)
        key = (qm_index, mm_index)
        if key in seen:
            continue
        seen.add(key)
        g = link_g_factor(z_of(qm_index), z_of(mm_index))
        if not (0.0 < g < 1.0):
            raise ValueError(
                "Unphysical QM/MM partition: link-atom scaling factor "
                f"g={g:.3f} for QM({qm_index})–MM({mm_index}). Reconsider "
                "which bond is cut by the QM/MM boundary."
            )
        links.append(LinkAtom(qm_index, mm_index, row_of[qm_index], g))

    links.sort(key=lambda la: (la.qm_index, la.mm_index))
    return links


def link_atom_position(qm_pos, mm_pos, g):
    """Cartesian position of a link atom on the QM→MM bond.

    ``qm_pos`` / ``mm_pos`` are length-3 sequences in any single length unit;
    the returned array is in that same unit.
    """
    qm_pos = np.asarray(qm_pos, dtype=float)
    mm_pos = np.asarray(mm_pos, dtype=float)
    return qm_pos + g * (mm_pos - qm_pos)


def project_link_gradient(grad_link, g):
    """Redistribute a link-atom gradient onto its (QM host, MM host).

    Returns ``(grad_to_qm_host, grad_to_mm_host)`` such that their sum equals
    ``grad_link`` (translational invariance of the capping constraint).
    """
    grad_link = np.asarray(grad_link, dtype=float)
    return (1.0 - g) * grad_link, g * grad_link
