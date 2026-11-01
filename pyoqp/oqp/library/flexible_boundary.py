"""FlexiBLE QM/MM boundary potential for diffusible solvent molecules.

Z. Shen and W. J. Glover, "Flexible boundary layer using exchange for
embedding theories. I. Theory and implementation", J. Chem. Phys. 155,
224112 (2021), doi:10.1063/5.0067855.

Each solvent molecule i of one molecular kind is represented by an anchor
(one atom, or the molecular centre of mass) at distance r_i from a boundary
centre.  For a partition of the molecules into a QM set Q and an MM set M,

    h(Q, M) = exp(-sum_{i in Q, j in M} g(r_i - r_j)),
    g(R)    = (alpha R)^3 / (1 + alpha R) for R > 0, and 0 otherwise,

so h = 1 when every QM molecule is closer to the centre than every MM
molecule.  The boundary potential is

    V = -k_B T ln[ h(actual partition) / sum_partitions h(partition) ],

where the denominator runs over the partitions with the same number of QM
molecules.  Only partitions near the distance-ordered one contribute; they
are collected as in the authors' OpenMM reference kernel
(github.com/kc5054/FlexiBLEPlugin, revision f7e5611): starting from the
ordered partition of the boundary window, a QM label moves one rank outwards
while the term stays above a threshold, and the threshold is lowered by a
fixed factor until the sum changes by no more than the initial threshold
relative to its previous value.  Unlike the reference kernel, the threshold
is lowered at least once even when the first sum is exactly 1, so a partition
just below the initial threshold is not dropped.  The force is the exact
derivative of the collected sum.

This module is an independent implementation in atomic units.  Periodic
cells are handled with the minimum-image anchor-to-centre vector, which the
reference kernel does not support; every QM molecule must then lie within
half the shortest box length of the centre.
"""

from __future__ import annotations

import math

import numpy as np

#: Largest boundary force (Hartree/bohr) allowed on an anchor that can reach
#: a periodic face of the minimum-image cell.
PERIODIC_FACE_FORCE_TOLERANCE = 1.0e-12

KB_HARTREE_PER_KELVIN = 3.166811563e-6
ANGSTROM_TO_BOHR = 1.0 / 0.52917721092


def pair_term(alpha, separation):
    """g(R) and dg/dR for an array of signed separations R = r_i - r_j."""
    separation = np.asarray(separation, dtype=np.float64)
    ar = alpha*np.where(separation > 0.0, separation, 0.0)
    value = ar**3/(1.0 + ar)
    derivative = alpha*(3.0*ar**2/(1.0 + ar) - ar**3/(1.0 + ar)**2)
    return value, derivative


class FlexiBLEBoundary:
    """One molecular kind (for example water) split into QM and MM molecules.

    ``qm_molecules`` and ``mm_molecules`` list the atom indices of each
    molecule.  ``anchor`` is the index of the anchor atom inside each molecule
    (the same for every molecule), or ``None`` for the molecular centre of
    mass.  ``center`` is ``("point", xyz_bohr)``, ``("atoms", indices)`` for
    the centre of mass of those atoms, or ``("molecules", None)`` for the
    centre of mass of all QM and MM molecules of this kind (the reference
    kernel's spherical boundary).  ``bonds`` lists atom pairs (for example
    all bonds of the topology) through which the centre atoms and, for
    ``anchor=None``, every molecule are made whole in a periodic cell;
    without it, or for an atom not bonded to the first atom of its group, the
    atom is placed at the minimum image of that first atom.
    """

    def __init__(self, qm_molecules, mm_molecules, masses, *, alpha,
                 temperature, anchor=0, center=("molecules", None),
                 threshold=0.01, scale=0.1, max_iterations=20,
                 keep_below_threshold=False, max_partitions=200000,
                 bonds=None):
        self.qm = [np.asarray(m, dtype=np.int64) for m in qm_molecules]
        self.mm = [np.asarray(m, dtype=np.int64) for m in mm_molecules]
        if not self.qm or not self.mm:
            raise ValueError(
                "FlexiBLE needs at least one QM and one MM molecule of the kind")
        sizes = {len(m) for m in self.qm + self.mm}
        if len(sizes) != 1:
            raise ValueError("FlexiBLE molecules of one kind must have equal size")
        self.molecules = self.qm + self.mm
        self.nqm = len(self.qm)
        self.nmol = len(self.molecules)
        masses = np.asarray(masses, dtype=np.float64).reshape(-1)
        self.masses = masses
        if anchor is not None:
            anchor = int(anchor)
            if not 0 <= anchor < sizes.pop():
                raise ValueError("FlexiBLE anchor is outside the molecule")
            # the integrator divides the force by the mass; a massless
            # (virtual-site) anchor would receive the whole boundary force
            anchor_masses = masses[[int(m[anchor]) for m in self.molecules]]
            if np.any(anchor_masses <= 0.0):
                raise ValueError(
                    "FlexiBLE anchor atom must carry mass; a massless virtual "
                    "site cannot be the anchor (use a real atom or com)")
        elif any(np.sum(masses[m]) <= 0.0 for m in self.molecules):
            raise ValueError("FlexiBLE molecule without mass")
        self.anchor = anchor
        self.alpha = float(alpha)
        self.temperature = float(temperature)
        self.threshold = float(threshold)
        self.scale = float(scale)
        self.max_iterations = int(max_iterations)
        self.keep_below_threshold = bool(keep_below_threshold)
        self.max_partitions = int(max_partitions)
        for name, value in (("alpha", self.alpha),
                            ("temperature", self.temperature)):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"FlexiBLE {name} must be positive and finite")
        if not 0.0 < self.threshold < 1.0 or not 0.0 < self.scale < 1.0:
            raise ValueError("FlexiBLE threshold and scale must lie in (0, 1)")
        if self.max_iterations < 2:
            raise ValueError("FlexiBLE max_iterations must be at least 2")
        kind, value = center
        if kind == "point":
            value = np.asarray(value, dtype=np.float64).reshape(3)
        elif kind == "atoms":
            value = np.asarray(value, dtype=np.int64).reshape(-1)
            if value.size == 0:
                raise ValueError("FlexiBLE centre atoms are empty")
            # the centre is mass-weighted; massless (virtual) sites alone
            # would give NaN weights
            if not np.sum(masses[value]) > 0.0:
                raise ValueError(
                    "FlexiBLE centre atoms must have a positive total mass; "
                    "massless (virtual) sites cannot define the centre")
        elif kind == "molecules":
            value = np.concatenate(self.molecules)
            # (rejected for a periodic cell in evaluate: the centre of mass
            # of a space-filling solvent is not defined there)
        else:
            raise ValueError(f"unknown FlexiBLE centre {kind!r}")
        self.center_kind = kind
        self.center_value = value
        neighbours = {}
        for i, j in (bonds if bonds is not None else ()):
            neighbours.setdefault(int(i), []).append(int(j))
            neighbours.setdefault(int(j), []).append(int(i))
        self.center_path = (self._unwrap_path(value, neighbours)
                            if kind == "atoms" else None)
        self.molecule_paths = (
            [self._unwrap_path(m, neighbours) for m in self.molecules]
            if self.anchor is None else None)
        self.last = None

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------

    @staticmethod
    def _image(vector, box):
        if box is None:
            return vector
        return vector - box*np.round(vector/box)

    @staticmethod
    def _unwrap_path(atoms, neighbours):
        """(parent, child) steps that make a group of atoms whole.

        A breadth-first tree over the bond graph ``neighbours`` from the first
        atom of the group; only the steps on the tree paths to group atoms
        are kept, in tree order.  An atom the bonds do not reach is stepped
        to directly from the first atom (its minimum image about that atom).
        """
        root = int(atoms[0])
        targets = set(int(a) for a in atoms)
        parent = {root: None}
        order = [root]
        if neighbours:
            head = 0
            while head < len(order) and not targets.issubset(parent):
                atom = order[head]
                head += 1
                for other in neighbours.get(atom, ()):
                    if other not in parent:
                        parent[other] = atom
                        order.append(other)
        needed = set()
        for atom in targets:
            while atom in parent and atom not in needed:
                needed.add(atom)
                atom = parent[atom]
        path = [(parent[a], a) for a in order if a in needed and parent[a] is not None]
        path += [(root, a) for a in sorted(targets) if a not in parent]
        return path

    def _center(self, xyz, box):
        """Centre position and, for a moving centre, (atoms, weights)."""
        if self.center_kind == "point":
            return self.center_value, None
        atoms = self.center_value
        weights = self.masses[atoms]/np.sum(self.masses[atoms])
        return weights @ self._whole(xyz, atoms, self.center_path, box), (
            atoms, weights)

    def _whole(self, xyz, atoms, path, box):
        """Positions of ``atoms`` made whole along their bond ``path``.

        A lattice translation of one atom does not change its force, so
        weights formed from these positions apply to the raw atoms.
        """
        if box is None:
            return xyz[atoms]
        whole = {int(atoms[0]): xyz[atoms[0]]}
        for parent, child in path:
            whole[child] = whole[parent] + self._image(
                xyz[child] - xyz[parent], box)
        return np.asarray([whole[int(a)] for a in atoms])

    def _anchors(self, xyz, box):
        """Anchor position of every molecule and its atom weights."""
        rows = []
        weights = []
        for index, molecule in enumerate(self.molecules):
            if self.anchor is not None:
                rows.append(xyz[molecule[self.anchor]])
                w = np.zeros(len(molecule))
                w[self.anchor] = 1.0
            else:
                m = self.masses[molecule]
                w = m/np.sum(m)
                rows.append(w @ self._whole(
                    xyz, molecule, self.molecule_paths[index], box))
            weights.append(w)
        return np.asarray(rows), weights

    # ------------------------------------------------------------------
    # partition sums
    # ------------------------------------------------------------------

    def _numerator(self, r):
        """h of the actual partition and its derivative with respect to r."""
        separation = r[:self.nqm, None] - r[None, self.nqm:]
        g, dg = pair_term(self.alpha, separation)
        value = math.exp(-float(np.sum(g)))
        derivative = np.zeros_like(r)
        derivative[:self.nqm] = -value*np.sum(dg, axis=1)
        derivative[self.nqm:] = value*np.sum(dg, axis=0)
        return value, derivative

    def _rank_weights(self, r_sorted):
        """h of each distance rank against its neighbours across the boundary.

        A QM rank p is weighed against the outer QM ranks and the first MM
        rank, an MM rank q against the last QM rank and the inner MM ranks;
        ranks whose weight is below the threshold bound the window in which
        partitions are enumerated.
        """
        nqm = self.nqm
        g = pair_term(self.alpha, r_sorted[None, :] - r_sorted[:, None])[0]
        weights = np.empty(self.nmol)
        for p in range(nqm):
            weights[p] = math.exp(-float(np.sum(g[p, p + 1:nqm + 1])))
        for q in range(nqm, self.nmol):
            weights[q] = math.exp(-float(np.sum(g[nqm - 1:q, q])))
        return weights

    def _window_bounds(self, weights, h):
        nqm = self.nqm
        lower = 0
        for p in range(nqm - 1, -1, -1):
            if weights[p] < h:
                lower = p + 1
                break
        # the reference kernel leaves this bound at 0 when no MM rank falls
        # below the threshold; every MM rank then belongs to the window
        upper = self.nmol - 1
        for q in range(nqm, self.nmol):
            if weights[q] < h:
                upper = q - 1
                break
        return lower, upper

    def _denominator(self, r, order, weights, h):
        """Sum of h over partitions reachable from the ordered window."""
        lower, upper = self._window_bounds(weights, h)
        nqm_window = self.nqm - lower
        nmm_window = upper - (self.nqm - 1)
        start = (1,)*nqm_window + (0,)*nmm_window
        window = order[lower:upper + 1]
        r_window = r[window]
        total = 0.0
        derivative = np.zeros_like(r)
        visited = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            bits = np.asarray(node, dtype=bool)
            separation = r_window[bits][:, None] - r_window[~bits][None, :]
            g, dg = pair_term(self.alpha, separation)
            value = math.exp(-float(np.sum(g)))
            if value < h and not self.keep_below_threshold:
                continue
            visited.add(node)
            if len(visited) > self.max_partitions:
                raise RuntimeError(
                    "FlexiBLE: more than "
                    f"{self.max_partitions} partitions exceed the threshold "
                    f"{h:.3e}; alpha is too small for the molecular spacing")
            total += value
            derivative[window[bits]] += -value*np.sum(dg, axis=1)
            derivative[window[~bits]] += value*np.sum(dg, axis=0)
            if value < h:
                continue
            for k in range(len(node) - 1):
                if node[k] == 1 and node[k + 1] == 0:
                    child = node[:k] + (0, 1) + node[k + 2:]
                    if child not in visited:
                        stack.append(child)
        return total, derivative

    def _check_periodic_faces(self, r, box, h, kt):
        """Refuse a cell in which an MM anchor on a face feels the boundary.

        The minimum-image radial direction of an anchor flips sign when it
        crosses a face of the cell, at r >= L/2 for the shortest edge L.  An
        MM molecule there lies on the QM side of a partition only together
        with a QM molecule on the MM side, so that partition's weight is at
        most w = exp(-g(L/2 - r_QM,max)), and its pair force at most
        kT g'(L/2 - r_QM,max) w.  Unless the threshold drops such partitions
        (w < h) or that force is negligible, the boundary force would jump
        when an MM molecule crosses a face, so the evaluation is refused.
        """
        separation = 0.5*float(np.min(box)) - float(np.max(r[:self.nqm]))
        g, dg = pair_term(self.alpha, separation)
        weight = math.exp(-float(g))
        if weight < h and not self.keep_below_threshold:
            return
        bound = kt*float(dg)*weight
        if bound > PERIODIC_FACE_FORCE_TOLERANCE:
            raise RuntimeError(
                "FlexiBLE: an MM molecule on a periodic face would still feel "
                "the boundary, so its force would jump when it crosses the "
                f"face (partition weight {weight:.3e} >= threshold {h:.3e}, "
                f"force up to {bound:.3e} Hartree/bohr); use a larger box or "
                "a larger alpha")

    # ------------------------------------------------------------------
    # potential
    # ------------------------------------------------------------------

    def evaluate(self, coordinates, box=None):
        """Energy (Hartree) and force (Hartree/bohr) on all atoms."""
        xyz = np.asarray(coordinates, dtype=np.float64)
        if box is not None:
            box = np.asarray(box, dtype=np.float64).reshape(3)
            if self.center_kind == "molecules":
                raise ValueError(
                    "FlexiBLE centre=molecules is not defined in a periodic "
                    "cell; centre the boundary on solute atoms or a point")
        center, center_weights = self._center(xyz, box)
        anchors, atom_weights = self._anchors(xyz, box)
        vectors = self._image(anchors - center, box)
        r = np.linalg.norm(vectors, axis=1)
        if box is not None and float(np.max(r[:self.nqm])) >= 0.5*float(np.min(box)):
            # A QM molecule there would make the boundary itself reach a
            # periodic image.  MM molecules there are checked below: they are
            # allowed only while their boundary force has decayed.
            far = int(np.argmax(r[:self.nqm]))
            raise RuntimeError(
                "FlexiBLE: a QM molecule is not closer to the centre than "
                "half the shortest box length; the boundary would reach its "
                f"own periodic image (QM molecule {far}, r = {r[far]:.6g} bohr)")
        # an anchor exactly on the centre has no radial direction: it gets
        # no force, but its distance stays 0 in the energy (the reference
        # kernel replaces it by 1e-8)
        zero = r <= 0.0
        units = np.where(zero[:, None], 0.0,
                         vectors/np.where(zero, 1.0, r)[:, None])

        numerator, d_numerator = self._numerator(r)
        if numerator < 1.0e-14:
            raise RuntimeError(
                "FlexiBLE: the QM/MM partition lies too far from the "
                f"distance ordering (h = {numerator:.3e})")
        order = np.argsort(r, kind="stable")
        weights = self._rank_weights(r[order])
        h = self.threshold
        gamma = self.threshold
        previous = None
        for iteration in range(1, self.max_iterations + 1):
            denominator, d_denominator = self._denominator(r, order, weights, h)
            if previous is None:
                # The reference kernel stops here when the sum is exactly 1.
                # A partition just below the initial threshold would then be
                # dropped, so the threshold is always lowered at least once.
                previous = denominator
                h *= self.scale
                continue
            if denominator - previous > gamma*previous:
                h *= self.scale
                previous = denominator
                continue
            break
        else:
            raise RuntimeError(
                "FlexiBLE: the partition sum did not converge in "
                f"{self.max_iterations} threshold reductions")

        kt = KB_HARTREE_PER_KELVIN*self.temperature
        energy = -kt*math.log(numerator/denominator)
        dlog = d_numerator/numerator - d_denominator/denominator
        anchor_force = kt*dlog[:, None]*units
        if box is not None:
            self._check_periodic_faces(r, box, h, kt)
        force = np.zeros_like(xyz)
        for molecule, weights, f in zip(self.molecules, atom_weights,
                                        anchor_force):
            force[molecule] += weights[:, None]*f[None, :]
        if center_weights is not None:
            atoms, weights = center_weights
            np.add.at(force, atoms, -weights[:, None]*np.sum(
                anchor_force, axis=0)[None, :])
        self.last = {
            "energy": energy,
            "force": force,
            "numerator": numerator,
            "denominator": denominator,
            "iterations": iteration,
            "threshold": h,
            "distances": r,
        }
        return energy, force
