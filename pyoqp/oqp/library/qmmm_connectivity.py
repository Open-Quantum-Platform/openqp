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


# ---------------------------------------------------------------------------
# Frontier (M1) charge redistribution for electrostatic embedding
# ---------------------------------------------------------------------------
#
# A covalent QM/MM cut leaves an MM *host* atom M1 (the MM end of the severed
# bond) roughly a bond length (~1.5 A) from the QM density and from the capping
# link-atom hydrogen.  Embedding M1's full force-field point charge there badly
# over-polarizes the QM density.  The remedy (Lin & Truhlar, J. Phys. Chem. A
# 2005, 109, 3991) is to delete q1 from M1 and redistribute it so the embedding
# charge cloud reproduces:
#
#   * the total charge of the deleted q1                        (RC and RCD), and
#   * its dipole about M1 (which is identically zero)           (RCD only),
#
# using *virtual* point charges placed at the midpoints of the M1-M2 bonds,
# where M2 are the MM atoms bonded to M1.  Let M1 have N such neighbours and let
# d_k = (R_M2k - R_M1)/2 be the displacement of midpoint k from M1:
#
#   Z1  : delete q1.                       Conserves neither charge nor dipole.
#   RC  : + q1/N at each midpoint.         Conserves charge; not dipole.
#   RCD : + 2 q1/N at each midpoint AND    Conserves charge AND dipole:
#         - q1/N added to each M2.           sum q = N(2 q1/N) - N(q1/N) = q1;
#                                            dipole = sum_k [2q1/N d_k - q1/N (2 d_k)] = 0.
#
# Each virtual charge sits at a *fixed linear combination* of two real atom
# positions (a bond midpoint), so - exactly like the link-atom capping position -
# any electrostatic force on it redistributes onto M1 and M2 by the chain rule.
# This keeps the analytic QM/MM gradient exact, which naive charge *deletion*
# does not (deletion drops a chunk of the field with no matching force term).
#
# NOTE (ESPF): OpenQP's electrostatic embedding is ESPF (Huix-Rotllant & Ferre,
# JCTC 2021, 17, 538). There the MM potential phi_A couples to the QM *atomic
# charge operators* Q_A (H += sum_A phi_A Q_A, eq 6), NOT to the raw electron
# density via 1/|r-R_M| point-charge integrals. That coarse-grained coupling
# already suppresses the electron spill-out / over-polarization that plagues
# density-based embedding, which is why the ESPF papers use the FULL MM charges
# with "no scaling of the QM-MM electrostatic interactions" even at a covalent
# protein boundary (cryptochrome). So for ESPF the DEFAULT is full-field
# ('none'); RCD/RC/Z1 are OPTIONAL refinements (physically motivated - they keep
# the MM charge/dipole faithful near the cut - but a departure from the validated
# ESPF baseline, so opt-in).


@dataclass(frozen=True)
class VirtualCharge:
    """A redistributed (virtual) embedding point charge at an M1-M2 bond
    midpoint.

    Attributes
    ----------
    charge : float
        Point charge (elementary charge units).
    hosts : tuple[int, int]
        Absolute atom indices ``(m1_index, m2_index)`` whose positions define
        the virtual site.
    weights : tuple[float, float]
        Barycentric weights so that ``r_virtual = w0*R[m1] + w1*R[m2]`` (``0.5,
        0.5`` for a bond midpoint).  A force ``f`` on the virtual charge maps to
        ``w0*f`` on ``m1`` and ``w1*f`` on ``m2``.
    """

    charge: float
    hosts: tuple
    weights: tuple


def redistribute_frontier_charges(frontier, charge_of, scheme="rcd"):
    """Redistribute MM frontier-host (M1) charges for electrostatic embedding.

    Parameters
    ----------
    frontier : iterable of (int, sequence[int])
        One ``(m1_index, [m2_index, ...])`` per *unique* MM host atom cut by the
        QM/MM boundary, with its MM neighbours M2.
    charge_of : callable(int) -> float
        MM force-field charge (e) of an absolute atom index.
    scheme : {'rcd', 'rc', 'z1', 'none'}
        Redistribution scheme (see module notes above). ``'none'`` disables it
        (full-field embedding, the pre-RCD behaviour).

    Returns
    -------
    deleted : set[int]
        MM atoms whose embedding charge is removed (the M1 hosts).
    delta_q : dict[int, float]
        Charge added to a *real* MM atom's embedding charge (the M2 corrections).
    virtuals : list[VirtualCharge]
        Virtual midpoint charges to add to the embedding set.
    """
    scheme = (scheme or "none").lower()
    if scheme not in ("rcd", "rc", "z1", "none"):
        raise ValueError(f"Unknown frontier charge scheme: {scheme!r}")

    deleted, delta_q, virtuals = set(), {}, []
    if scheme == "none":
        return deleted, delta_q, virtuals

    for m1, m2s in frontier:
        m1 = int(m1)
        m2s = [int(m) for m in m2s]
        q1 = float(charge_of(m1))
        deleted.add(m1)
        # Nothing to redistribute onto (isolated host) -> fall back to deletion.
        if scheme == "z1" or not m2s:
            continue
        n = len(m2s)
        for m2 in m2s:
            if scheme == "rc":
                virtuals.append(VirtualCharge(q1 / n, (m1, m2), (0.5, 0.5)))
            else:  # rcd
                virtuals.append(VirtualCharge(2.0 * q1 / n, (m1, m2), (0.5, 0.5)))
                delta_q[m2] = delta_q.get(m2, 0.0) - q1 / n
    return deleted, delta_q, virtuals


def assemble_embedding_sites(mm_idx, mm_charges, mm_positions,
                             deleted, delta_q, virtuals):
    """Assemble the embedding point-charge set the QM density is embedded in.

    Combines the raw MM atoms with the frontier redistribution from
    :func:`redistribute_frontier_charges`.

    Parameters
    ----------
    mm_idx : sequence[int]
        Absolute atom indices of the real MM atoms, in array order.
    mm_charges : sequence[float]
        Force-field charges (e), aligned with ``mm_idx``.
    mm_positions : (M, 3) array
        Positions (any length unit), aligned with ``mm_idx``.
    deleted, delta_q, virtuals :
        Output of :func:`redistribute_frontier_charges`.

    Returns
    -------
    charges : (S,) ndarray
        Embedding charges.
    positions : (S, 3) ndarray
        Embedding positions (same unit as ``mm_positions``).
    scatter : list[tuple[tuple[int, float], ...]]
        For each site, the ``(absolute_atom_index, weight)`` contributions that
        build its position, so a force ``f`` on site ``s`` adds ``weight*f`` to
        each contributing real atom (identity for a real atom, 0.5/0.5 for a
        midpoint).
    """
    mm_positions = np.asarray(mm_positions, dtype=float)
    pos_of = {int(mm_idx[k]): mm_positions[k] for k in range(len(mm_idx))}

    charges, positions, scatter = [], [], []
    applied = set()
    for k in range(len(mm_idx)):
        a = int(mm_idx[k])
        if a in deleted:
            continue
        charges.append(float(mm_charges[k]) + delta_q.get(a, 0.0))
        positions.append(mm_positions[k])
        scatter.append(((a, 1.0),))
        applied.add(a)
    # A delta_q correction whose target atom is itself a deleted host (two
    # adjacent cuts / mutually-bonded frontier hosts) would be dropped by the
    # loop above. Keep it as a charge-only residual site at that atom so total
    # charge and the dipole about each host stay conserved regardless of topology.
    for a, dq in delta_q.items():
        a = int(a)
        if a not in applied and dq != 0.0:
            charges.append(dq)
            positions.append(pos_of[a])
            scatter.append(((a, 1.0),))
    for v in virtuals:
        m1, m2 = v.hosts
        w0, w1 = v.weights
        positions.append(w0 * pos_of[m1] + w1 * pos_of[m2])
        charges.append(v.charge)
        scatter.append(((m1, w0), (m2, w1)))
    return (np.asarray(charges, dtype=float),
            np.asarray(positions, dtype=float),
            scatter)
