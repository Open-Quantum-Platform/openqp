"""Coordinate systems for the builtin OpenQP geometry optimizer.

This module is intentionally backend-agnostic: it depends only on NumPy and can
be unit-tested without the compiled OQP electronic-structure core.  It provides

* :class:`RedundantInternalCoordinates` -- bonds, angles and dihedrals built
  from covalent connectivity, with the Wilson B-matrix, gradient/Hessian
  transforms and an iterative Cartesian back-transformation
  (Pulay/Fogarasi natural internals; Peng, Ayala, Schlegel, Frisch,
  J. Comput. Chem. 17, 49 (1996)); and
* :class:`CartesianCoordinates` -- a trivial fall-back used when the internal
  set fails to span the 3N-6 non-redundant space (e.g. linear molecules) or the
  back-transformation diverges.

Both classes share the same small interface consumed by
``builtin_engine.BuiltinEngine``::

    coords.q(x)                 -> internal/Cartesian coordinate vector
    coords.b_matrix(x)          -> dq/dx        (n_q, 3N)
    coords.grad_to_q(x, gx)     -> gradient in working coordinates
    coords.guess_hessian(x)     -> diagonal model Hessian (n_q, n_q)
    coords.q_displacement(q2,q1)-> q2 - q1 with dihedral wrapping
    coords.back_transform(x, dq)-> new Cartesian geometry for a step dq
    coords.cart_rmsd(x, x_new)  -> RMS Cartesian displacement (Bohr)

All quantities are in atomic units (Bohr, radian, Hartree/Bohr ...).
"""

from __future__ import annotations

import numpy as np

BOHR_TO_ANGSTROM = 0.52917721092
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

# Cordero et al. covalent radii (Angstrom), Dalton Trans. 2008, 2832.
# Indexed by atomic number; entry 0 is a placeholder.
_COVALENT_RADII_ANG = [
    0.00,
    0.31, 0.28,                                                  # H  He
    1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,              # Li..Ne
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,              # Na..Ar
    2.03, 1.76,                                                  # K  Ca
    1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,  # Sc..Zn
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16,                          # Ga..Kr
    2.20, 1.95,                                                  # Rb Sr
    1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,  # Y..Cd
    1.42, 1.39, 1.39, 1.38, 1.39, 1.40,                          # In..Xe
]


def covalent_radius_bohr(z: int) -> float:
    """Covalent radius of element ``z`` in Bohr (1.5 A fallback for heavy Z)."""
    if 0 < z < len(_COVALENT_RADII_ANG):
        r = _COVALENT_RADII_ANG[z]
    else:
        r = 1.5
    return r * ANGSTROM_TO_BOHR


# --------------------------------------------------------------------------- #
# Primitive internal coordinates                                              #
# --------------------------------------------------------------------------- #
# Each primitive exposes ``value(x)`` and ``derivatives(x)`` where ``x`` has
# shape (natom, 3).  ``derivatives`` returns a list of (atom_index, 3-vector)
# pairs giving dq/dr_atom; all other atoms have zero derivative.


class Bond:
    kind = "bond"

    def __init__(self, i, j):
        self.atoms = (i, j)

    def value(self, x):
        return float(np.linalg.norm(x[self.atoms[0]] - x[self.atoms[1]]))

    def derivatives(self, x):
        i, j = self.atoms
        u = x[i] - x[j]
        r = np.linalg.norm(u)
        u = u / r
        return [(i, u), (j, -u)]


class Angle:
    kind = "angle"

    def __init__(self, i, j, k):
        # j is the apex
        self.atoms = (i, j, k)

    def value(self, x):
        i, j, k = self.atoms
        u = x[i] - x[j]
        v = x[k] - x[j]
        u /= np.linalg.norm(u)
        v /= np.linalg.norm(v)
        c = float(np.clip(np.dot(u, v), -1.0, 1.0))
        return float(np.arccos(c))

    def derivatives(self, x):
        i, j, k = self.atoms
        u = x[i] - x[j]
        v = x[k] - x[j]
        lu = np.linalg.norm(u)
        lv = np.linalg.norm(v)
        u = u / lu
        v = v / lv
        c = float(np.clip(np.dot(u, v), -1.0, 1.0))
        s = np.sqrt(max(1.0 - c * c, 1.0e-12))
        di = (c * u - v) / (lu * s)
        dk = (c * v - u) / (lv * s)
        dj = -(di + dk)
        return [(i, di), (j, dj), (k, dk)]


class Dihedral:
    kind = "dihedral"

    def __init__(self, i, j, k, l):
        self.atoms = (i, j, k, l)

    def value(self, x):
        i, j, k, l = self.atoms
        f = x[i] - x[j]
        g = x[j] - x[k]
        h = x[l] - x[k]
        a = np.cross(f, g)
        b = np.cross(h, g)
        lg = np.linalg.norm(g)
        return float(np.arctan2(lg * np.dot(f, b), np.dot(a, b)))

    def derivatives(self, x, h=1.0e-6):
        # The end-atom dihedral derivatives are well known, but the central-atom
        # terms are notoriously sign-error prone.  For geometry-optimization
        # B-matrices a central finite difference of the (analytic, periodic-safe)
        # value is both exact to ~1e-9 and unambiguous, at negligible cost for
        # the small systems these optimizations target.
        i, j, k, l = self.atoms
        out = []
        for atom in (i, j, k, l):
            d = np.zeros(3)
            for c in range(3):
                xp = x.copy()
                xm = x.copy()
                xp[atom, c] += h
                xm[atom, c] -= h
                dphi = self.value(xp) - self.value(xm)
                # unwrap the +/-2*pi branch jump
                dphi = (dphi + np.pi) % (2.0 * np.pi) - np.pi
                d[c] = dphi / (2.0 * h)
            out.append((atom, d))
        return out


# --------------------------------------------------------------------------- #
# Translation-rotation coordinates (TRIC, Wang & Song, JCP 144, 214108 (2016)) #
# --------------------------------------------------------------------------- #
def _rotmat_to_expmap(R):
    """Exponential-map vector v = theta * axis of a 3x3 rotation matrix."""
    cos = float(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))
    theta = np.arccos(cos)
    skew = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    if theta < 1.0e-8:
        return 0.5 * skew                      # small-angle limit
    if abs(theta - np.pi) < 1.0e-6:
        # 180 deg: axis from the eigenvector of R with eigenvalue +1 (rare).
        w, v = np.linalg.eig(R)
        axis = np.real(v[:, int(np.argmin(np.abs(w - 1.0)))])
        axis = axis / (np.linalg.norm(axis) + 1.0e-12)
        return theta * axis
    return theta * skew / (2.0 * np.sin(theta))


def _rotation_vector(X, X0):
    """Exp-map vector of the rotation best aligning reference ``X0`` onto ``X``.

    ``X``/``X0`` are (N,3) fragment coordinates (any units).  Uses the Kabsch
    superposition; returns a 3-vector that is zero when X and X0 share an
    orientation.
    """
    Xc = X - X.mean(axis=0)
    X0c = X0 - X0.mean(axis=0)
    h = X0c.T @ Xc
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return _rotmat_to_expmap(rot)


class TranslationComponent:
    """One Cartesian centroid component of a fragment (analytic, constant B)."""

    kind = "translation"

    def __init__(self, atoms, axis):
        self.frag = list(atoms)
        self.axis = int(axis)
        self.n = len(self.frag)

    def value(self, x):
        return float(np.mean(np.asarray(x)[self.frag, self.axis]))

    def derivatives(self, x):
        d = np.zeros(3)
        d[self.axis] = 1.0 / self.n
        return [(a, d.copy()) for a in self.frag]


class RotationComponent:
    """One exp-map rotation component of a fragment relative to a reference.

    The value is evaluated analytically (Kabsch); the B-matrix row is taken by a
    central finite difference of that value -- exact to ~1e-9 and free of the
    intricate quaternion-derivative algebra, at negligible cost for the small
    systems these optimizations target.
    """

    kind = "rotation"

    def __init__(self, atoms, axis, reference_xyz):
        self.frag = list(atoms)
        self.axis = int(axis)
        self.ref = np.asarray(reference_xyz, dtype=float).reshape(-1, 3)[self.frag].copy()

    def value(self, x):
        sub = np.asarray(x, dtype=float).reshape(-1, 3)[self.frag]
        return float(_rotation_vector(sub, self.ref)[self.axis])

    def derivatives(self, x, h=1.0e-5):
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        out = []
        for a in self.frag:
            d = np.zeros(3)
            for c in range(3):
                xp = x.copy()
                xm = x.copy()
                xp[a, c] += h
                xm[a, c] -= h
                d[c] = (self.value(xp) - self.value(xm)) / (2.0 * h)
            out.append((a, d))
        return out


# Diagonal model force constants (Hartree / Bohr**2 or / rad**2), following the
# diagonal guess used by geomeTRIC (Schlegel-type values).  Translation and
# rotation use the small geomeTRIC trans/rot value.
_GUESS_FC = {"bond": 0.35, "angle": 0.20, "dihedral": 0.023,
             "translation": 0.05, "rotation": 0.05}


class RedundantInternalCoordinates:
    """Redundant internal coordinate system (bonds, angles, dihedrals)."""

    def __init__(self, primitives, natom):
        self.primitives = primitives
        self.natom = natom
        self.ndim = 3 * natom

    # -- construction ------------------------------------------------------- #
    @classmethod
    def from_geometry(cls, atoms, x, bond_factor=1.3):
        """Build internals from atomic numbers ``atoms`` and geometry ``x``.

        Returns ``None`` when fewer than two atoms are present (nothing to do)
        or when the constructed set is rank-deficient (caller falls back to
        Cartesians).
        """
        atoms = np.asarray(atoms, dtype=int).reshape(-1)
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        natom = len(atoms)
        if natom < 2:
            return None

        radii = np.array([covalent_radius_bohr(z) for z in atoms])

        # Bonds from covalent connectivity.
        neighbours = {a: set() for a in range(natom)}
        bonds = []
        for a in range(natom):
            for b in range(a + 1, natom):
                rab = np.linalg.norm(x[a] - x[b])
                if rab < bond_factor * (radii[a] + radii[b]):
                    bonds.append((a, b))
                    neighbours[a].add(b)
                    neighbours[b].add(a)

        # Connect disjoint fragments by their closest interatomic contact so the
        # internal set spans intermolecular translations/rotations.
        bonds = _connect_fragments(bonds, neighbours, x, natom)

        primitives = [Bond(a, b) for (a, b) in bonds]

        # Angles: apex j with neighbour pair (i, k); skip near-linear.
        for j in range(natom):
            nb = sorted(neighbours[j])
            for ii in range(len(nb)):
                for kk in range(ii + 1, len(nb)):
                    i, k = nb[ii], nb[kk]
                    ang = Angle(i, j, k)
                    val = ang.value(x)
                    if 0.26 < val < (np.pi - 0.26):  # ~15 deg .. 165 deg
                        primitives.append(ang)

        # Dihedrals along each bond (j-k) with terminal neighbours i and l.
        seen = set()
        for (j, k) in bonds:
            for i in sorted(neighbours[j]):
                if i == k:
                    continue
                if not _well_defined_angle(x, i, j, k):
                    continue
                for l in sorted(neighbours[k]):
                    if l == j or l == i:
                        continue
                    if not _well_defined_angle(x, j, k, l):
                        continue
                    key = (i, j, k, l) if (i, j) < (l, k) else (l, k, j, i)
                    if key in seen:
                        continue
                    seen.add(key)
                    primitives.append(Dihedral(i, j, k, l))

        coords = cls(primitives, natom)
        # Reject rank-deficient sets (e.g. linear species) -> Cartesian fallback.
        if not coords.spans_internal_space(x):
            return None
        return coords

    # -- core linear algebra ------------------------------------------------ #
    def q(self, x):
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        return np.array([p.value(x) for p in self.primitives])

    def b_matrix(self, x):
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        nq = len(self.primitives)
        b = np.zeros((nq, self.ndim))
        for row, p in enumerate(self.primitives):
            for atom, d in p.derivatives(x):
                b[row, 3 * atom:3 * atom + 3] = d
        return b

    def _g_inverse(self, b, eig_tol=1.0e-6):
        g = b @ b.T
        w, v = np.linalg.eigh(g)
        keep = w > eig_tol
        inv = (v[:, keep] / w[keep]) @ v[:, keep].T
        rank = int(np.count_nonzero(keep))
        return inv, rank

    def spans_internal_space(self, x, eig_tol=1.0e-6):
        b = self.b_matrix(x)
        _, rank = self._g_inverse(b, eig_tol)
        target = 3 * self.natom - 6 if self.natom > 2 else 3 * self.natom - 5
        return rank >= target

    def grad_to_q(self, x, gx):
        b = self.b_matrix(x)
        ginv, _ = self._g_inverse(b)
        return ginv @ (b @ np.asarray(gx, dtype=float).reshape(-1))

    def guess_hessian(self, x):
        diag = np.array([_GUESS_FC[p.kind] for p in self.primitives])
        return np.diag(diag)

    def q_displacement(self, q2, q1):
        dq = np.asarray(q2, dtype=float) - np.asarray(q1, dtype=float)
        for idx, p in enumerate(self.primitives):
            if p.kind == "dihedral":
                dq[idx] = (dq[idx] + np.pi) % (2.0 * np.pi) - np.pi
        return dq

    def back_transform(self, x, dq, tol=1.0e-6, maxiter=50):
        """Iteratively realise an internal-coordinate step ``dq`` in Cartesians.

        Returns ``(x_new, ok)``; ``ok`` is False if the iteration diverged, in
        which case the caller should shrink the step or fall back.
        """
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        q_target = self.q(x) + dq
        x_cur = x.copy()
        prev_norm = None
        best = x_cur.copy()
        for _ in range(maxiter):
            b = self.b_matrix(x_cur)
            ginv, _ = self._g_inverse(b)
            dq_remain = self.q_displacement(q_target, self.q(x_cur))
            err = np.linalg.norm(dq_remain)
            if prev_norm is not None and err > prev_norm * 1.0001 + 1.0e-10:
                # Diverging: return the best Cartesian step found so far.
                return best.reshape(-1), False
            prev_norm = err
            best = x_cur.copy()
            if err < tol:
                return x_cur.reshape(-1), True
            dx = (b.T @ (ginv @ dq_remain)).reshape(-1, 3)
            x_cur = x_cur + dx
        # Did not fully converge; accept the closest geometry but flag it.
        return best.reshape(-1), prev_norm is not None and prev_norm < 1.0e-3

    def cart_rmsd(self, x, x_new):
        a = np.asarray(x, dtype=float).reshape(-1)
        b = np.asarray(x_new, dtype=float).reshape(-1)
        return float(np.sqrt(np.mean((a - b) ** 2)))


class DelocalizedInternalCoordinates:
    """Delocalized internal coordinates (Baker, JCP 105, 192 (1996)).

    Non-redundant linear combinations of a redundant primitive set: the active
    subspace ``U`` is the eigenvectors of G = B B^T with non-zero eigenvalue,
    fixed at construction.  Coordinates are measured as displacements from the
    construction geometry, so G_dlc is well conditioned and the transforms use a
    plain solve instead of a pseudo-inverse.
    """

    def __init__(self, ric, u, q_ref):
        self.ric = ric
        self.U = u
        self.q_ref = q_ref
        self.natom = ric.natom
        self.ndim = ric.ndim
        self.primitives = ric.primitives

    @classmethod
    def from_ric(cls, ric, x, eig_tol=1.0e-6):
        b = ric.b_matrix(x)
        g = b @ b.T
        w, v = np.linalg.eigh(0.5 * (g + g.T))
        u = v[:, w > eig_tol]
        return cls(ric, u, ric.q(x))

    def _disp_prim(self, x):
        return self.ric.q_displacement(self.ric.q(x), self.q_ref)

    def q(self, x):
        return self.U.T @ self._disp_prim(x)

    def b_matrix(self, x):
        return self.U.T @ self.ric.b_matrix(x)

    def _g_inverse(self, b, eig_tol=1.0e-8):
        g = b @ b.T
        w, v = np.linalg.eigh(0.5 * (g + g.T))
        keep = w > eig_tol
        inv = (v[:, keep] / w[keep]) @ v[:, keep].T
        return inv, int(np.count_nonzero(keep))

    def grad_to_q(self, x, gx):
        b = self.b_matrix(x)
        ginv, _ = self._g_inverse(b)
        return ginv @ (b @ np.asarray(gx, dtype=float).reshape(-1))

    def guess_hessian(self, x):
        diag = np.array([_GUESS_FC[p.kind] for p in self.ric.primitives])
        return self.U.T @ (diag[:, None] * self.U)

    def q_displacement(self, q2, q1):
        return np.asarray(q2, dtype=float) - np.asarray(q1, dtype=float)

    def back_transform(self, x, dq, tol=1.0e-6, maxiter=50):
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        s_target = self.q(x) + dq
        x_cur = x.copy()
        prev = None
        best = x_cur.copy()
        for _ in range(maxiter):
            b = self.b_matrix(x_cur)
            ginv, _ = self._g_inverse(b)
            ds = s_target - self.q(x_cur)
            err = np.linalg.norm(ds)
            if prev is not None and err > prev * 1.0001 + 1.0e-10:
                return best.reshape(-1), False
            prev = err
            best = x_cur.copy()
            if err < tol:
                return x_cur.reshape(-1), True
            x_cur = x_cur + (b.T @ (ginv @ ds)).reshape(-1, 3)
        return best.reshape(-1), prev is not None and prev < 1.0e-3

    def cart_rmsd(self, x, x_new):
        a = np.asarray(x, dtype=float).reshape(-1)
        b = np.asarray(x_new, dtype=float).reshape(-1)
        return float(np.sqrt(np.mean((a - b) ** 2)))


def _build_tric(atoms, x, bond_factor=1.3):
    """TRIC: bonded internals + per-fragment translation/rotation coordinates."""
    atoms = np.asarray(atoms, dtype=int).reshape(-1)
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    natom = len(atoms)
    if natom < 2:
        return None

    bonds, neighbours, fragments = _covalent_graph(atoms, x, bond_factor)
    primitives = _make_internals(bonds, neighbours, x)

    for frag in fragments:
        for axis in range(3):
            primitives.append(TranslationComponent(frag, axis))
        if len(frag) >= 3:
            sub = x[frag] - x[frag].mean(axis=0)
            sv = np.linalg.svd(sub, compute_uv=False)
            if len(sv) >= 2 and sv[1] > 1.0e-3:        # non-collinear fragment
                for axis in range(3):
                    primitives.append(RotationComponent(frag, axis, x))

    coords = RedundantInternalCoordinates(primitives, natom)
    # TRIC must span the FULL 3N space (6 external + 3N-6 internal per the
    # connected fragments).  Requiring only 3N-6 would let the padded
    # translations mask a missing internal mode (e.g. the bend of a linear
    # molecule), so demand full rank and otherwise fall back.
    _, rank = coords._g_inverse(coords.b_matrix(x))
    if rank < 3 * natom:
        return None
    return coords


class CartesianCoordinates:
    """Trivial Cartesian coordinate system (identity B-matrix)."""

    def __init__(self, natom):
        self.natom = natom
        self.ndim = 3 * natom
        self.primitives = []

    def q(self, x):
        return np.asarray(x, dtype=float).reshape(-1)

    def b_matrix(self, x):
        return np.eye(self.ndim)

    def grad_to_q(self, x, gx):
        return np.asarray(gx, dtype=float).reshape(-1)

    def guess_hessian(self, x):
        return np.eye(self.ndim) * 0.5

    def q_displacement(self, q2, q1):
        return np.asarray(q2, dtype=float) - np.asarray(q1, dtype=float)

    def back_transform(self, x, dq, **_):
        x = np.asarray(x, dtype=float).reshape(-1)
        return x + np.asarray(dq, dtype=float).reshape(-1), True

    def cart_rmsd(self, x, x_new):
        a = np.asarray(x, dtype=float).reshape(-1)
        b = np.asarray(x_new, dtype=float).reshape(-1)
        return float(np.sqrt(np.mean((a - b) ** 2)))


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #
def _well_defined_angle(x, i, j, k, tol=0.26):
    a = Angle(i, j, k).value(x)
    return tol < a < (np.pi - tol)


def _connect_fragments(bonds, neighbours, x, natom):
    """Add auxiliary bonds so the connectivity graph is a single component."""
    comp = _components(neighbours, natom)
    while len(comp) > 1:
        # Merge the two closest fragments by their nearest atom pair.
        best = None
        for ci in range(len(comp)):
            for cj in range(ci + 1, len(comp)):
                for a in comp[ci]:
                    for b in comp[cj]:
                        r = np.linalg.norm(x[a] - x[b])
                        if best is None or r < best[0]:
                            best = (r, a, b)
        _, a, b = best
        bonds.append((a, b) if a < b else (b, a))
        neighbours[a].add(b)
        neighbours[b].add(a)
        comp = _components(neighbours, natom)
    return bonds


def _components(neighbours, natom):
    seen = set()
    comps = []
    for start in range(natom):
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            a = stack.pop()
            if a in seen:
                continue
            seen.add(a)
            comp.append(a)
            stack.extend(neighbours[a] - seen)
        comps.append(comp)
    return comps


def _covalent_graph(atoms, x, bond_factor=1.3):
    """Covalent bonds, neighbour map, and connected fragments (no auxiliaries)."""
    atoms = np.asarray(atoms, dtype=int).reshape(-1)
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    natom = len(atoms)
    radii = np.array([covalent_radius_bohr(z) for z in atoms])
    neighbours = {a: set() for a in range(natom)}
    bonds = []
    for a in range(natom):
        for b in range(a + 1, natom):
            if np.linalg.norm(x[a] - x[b]) < bond_factor * (radii[a] + radii[b]):
                bonds.append((a, b))
                neighbours[a].add(b)
                neighbours[b].add(a)
    return bonds, neighbours, _components(neighbours, natom)


def _make_internals(bonds, neighbours, x):
    """Bond/angle/dihedral primitives for a given covalent graph."""
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    primitives = [Bond(a, b) for (a, b) in bonds]
    for j in neighbours:
        nb = sorted(neighbours[j])
        for ii in range(len(nb)):
            for kk in range(ii + 1, len(nb)):
                ang = Angle(nb[ii], j, nb[kk])
                if 0.26 < ang.value(x) < (np.pi - 0.26):
                    primitives.append(ang)
    seen = set()
    for (j, k) in bonds:
        for i in sorted(neighbours[j]):
            if i == k or not _well_defined_angle(x, i, j, k):
                continue
            for l in sorted(neighbours[k]):
                if l == j or l == i or not _well_defined_angle(x, j, k, l):
                    continue
                key = (i, j, k, l) if (i, j) < (l, k) else (l, k, j, i)
                if key in seen:
                    continue
                seen.add(key)
                primitives.append(Dihedral(i, j, k, l))
    return primitives


def build_coordinates(atoms, x, coordsys="tric"):
    """Return a coordinate system for ``atoms``/``x``.

    ``coordsys`` selects the working coordinates:

    * ``tric`` (default) / ``auto`` -- translation-rotation internal coordinates;
    * ``dlc`` -- delocalized internal coordinates;
    * ``ric`` / ``internal`` -- redundant internal coordinates;
    * ``cart`` / ``cartesian`` -- Cartesians.

    All internal systems fall back to redundant internals and then Cartesians
    when their set is unavailable or rank-deficient (e.g. linear species), so the
    optimizer always has a usable coordinate system.
    """
    natom = len(np.asarray(atoms, dtype=int).reshape(-1))
    cs = (coordsys or "tric").lower()

    if cs in ("cart", "cartesian"):
        return CartesianCoordinates(natom)

    ric = RedundantInternalCoordinates.from_geometry(atoms, x)

    if cs in ("ric", "internal"):
        return ric if ric is not None else CartesianCoordinates(natom)

    if cs == "dlc":
        if ric is None:
            return CartesianCoordinates(natom)
        try:
            return DelocalizedInternalCoordinates.from_ric(ric, x)
        except Exception:
            return ric

    # tric / auto / anything else
    tric = _build_tric(atoms, x)
    if tric is not None:
        return tric
    return ric if ric is not None else CartesianCoordinates(natom)
