"""Molden orbital parsing and cube-grid generation.

Ports the parsing and evaluation logic of OpenqpView (frontend viewer) to
Python + NumPy so molecular orbitals from OpenQP Molden files can be
rendered as cube isosurfaces by Mol*. Supports Cartesian s/p/d/f/g shells;
spherical ([5D]/[7F]/[9G]) files are reported as unsupported.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903

SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
ATOMIC_NUMBER = {symbol: z for z, symbol in enumerate(SYMBOLS)}

# Molden Cartesian component ordering per shell.
CART_POWERS: dict[str, list[tuple[int, int, int]]] = {
    "s": [(0, 0, 0)],
    "p": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "d": [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    "f": [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
          (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
    "g": [(4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 1, 0), (3, 0, 1),
          (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3), (2, 2, 0),
          (2, 0, 2), (0, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)],
}


def _double_factorial(n: int) -> float:
    result = 1.0
    while n > 1:
        result *= n
        n -= 2
    return result


def _primitive_norm(alpha: float, ang: int) -> float:
    """Normalization of one primitive Cartesian Gaussian of type (l, 0, 0)."""
    return ((2.0 * alpha / math.pi) ** 0.75
            * (4.0 * alpha) ** (ang / 2.0)
            / math.sqrt(_double_factorial(2 * ang - 1)))


def _component_scale(powers: tuple[int, int, int]) -> float:
    """Cartesian component normalization relative to the pure (l,0,0) form."""
    total = sum(powers)
    if total < 2:
        return 1.0
    component = 1.0
    for p in powers:
        component *= math.factorial(p) / math.factorial(2 * p)
    pure = math.factorial(total) / math.factorial(2 * total)
    return math.sqrt(component / pure)


@dataclass
class BasisFunction:
    atom_index: int
    powers: tuple[int, int, int]
    exponents: np.ndarray
    coefficients: np.ndarray
    scale: float


@dataclass
class Orbital:
    index: int
    energy: float
    spin: str
    occupancy: float | None
    symmetry: str | None
    coefficients: np.ndarray


@dataclass
class MoldenData:
    atoms: list[tuple[str, float, float, float]]  # Angstrom
    basis: list[BasisFunction] = field(default_factory=list)
    orbitals: list[Orbital] = field(default_factory=list)
    supported: bool = True
    unsupported: list[str] = field(default_factory=list)


def parse_molden(text: str) -> MoldenData:
    lines = text.splitlines()

    # Atoms
    atoms: list[tuple[str, float, float, float]] = []
    start = next((i for i, ln in enumerate(lines) if re.match(r"^\[Atoms\]", ln.strip(), re.IGNORECASE)), -1)
    factor = 1.0
    if start >= 0:
        unit = re.search(r"\]\s*(\S+)", lines[start])
        if unit and unit.group(1).upper() == "AU":
            factor = BOHR_TO_ANGSTROM
        for ln in lines[start + 1 :]:
            if re.match(r"^\s*\[", ln):
                break
            parts = ln.split()
            if len(parts) < 6:
                continue
            atoms.append((parts[0], float(parts[3]) * factor,
                          float(parts[4]) * factor, float(parts[5]) * factor))
    data = MoldenData(atoms=atoms)

    # Spherical markers make the Cartesian evaluator invalid.
    spherical = set()
    for ln in lines:
        marker = ln.strip()
        if re.match(r"^\[[^\]]*5D", marker, re.IGNORECASE):
            spherical.add("d")
        if re.match(r"^\[[^\]]*7F", marker, re.IGNORECASE):
            spherical.add("f")
        if re.match(r"^\[[^\]]*9G", marker, re.IGNORECASE):
            spherical.add("g")

    # GTO basis
    start = next((i for i, ln in enumerate(lines) if re.match(r"^\[GTO\]", ln.strip(), re.IGNORECASE)), -1)
    if start >= 0:
        current_atom: int | None = None
        i = start + 1
        while i < len(lines):
            line = lines[i].strip()
            if re.match(r"^\[", line):
                break
            atom_match = re.match(r"^(\d+)(?:\s+0)?\s*$", line)
            shell_match = re.match(r"^([spdfgh])\s+(\d+)", line, re.IGNORECASE)
            if atom_match:
                current_atom = int(atom_match.group(1)) - 1
            elif shell_match and current_atom is not None:
                shell = shell_match.group(1).lower()
                nprim = int(shell_match.group(2))
                exps, coefs = [], []
                for p in range(nprim):
                    parts = lines[i + 1 + p].split()
                    exps.append(float(parts[0].replace("D", "E").replace("d", "e")))
                    coefs.append(float(parts[1].replace("D", "E").replace("d", "e")))
                if shell in spherical or shell not in CART_POWERS:
                    data.unsupported.append(shell.upper())
                else:
                    # Molden contraction coefficients refer to normalized
                    # primitives, so fold that normalization in here once.
                    ang = sum(CART_POWERS[shell][0])
                    coefs = [c * _primitive_norm(a, ang)
                             for a, c in zip(exps, coefs)]
                    for powers in CART_POWERS[shell]:
                        data.basis.append(BasisFunction(
                            atom_index=current_atom,
                            powers=powers,
                            exponents=np.asarray(exps),
                            coefficients=np.asarray(coefs),
                            scale=_component_scale(powers),
                        ))
                i += nprim
            i += 1
    if data.unsupported:
        data.supported = False
        return data

    # MO coefficients
    start = next((i for i, ln in enumerate(lines) if re.match(r"^\[MO\]", ln.strip(), re.IGNORECASE)), -1)
    if start >= 0:
        nao = len(data.basis)
        current: Orbital | None = None
        pending_sym: str | None = None
        for ln in lines[start + 1 :]:
            line = ln.strip()
            if re.match(r"^\[", line):
                break
            if not line:
                continue
            if m := re.match(r"^Sym=\s*(\S+)", line, re.IGNORECASE):
                pending_sym = m.group(1)
                continue
            if m := re.match(r"^Ene=\s*([-+0-9.EeDd]+)", line, re.IGNORECASE):
                current = Orbital(
                    index=len(data.orbitals) + 1,
                    energy=float(m.group(1).replace("D", "E").replace("d", "e")),
                    spin="Unknown",
                    occupancy=None,
                    symmetry=pending_sym,
                    coefficients=np.zeros(nao),
                )
                data.orbitals.append(current)
                pending_sym = None
                continue
            if current is None:
                continue
            if m := re.match(r"^Spin=\s*(\S+)", line, re.IGNORECASE):
                current.spin = m.group(1)
                continue
            if m := re.match(r"^Occup=\s*([-+0-9.EeDd]+)", line, re.IGNORECASE):
                current.occupancy = float(m.group(1).replace("D", "E").replace("d", "e"))
                continue
            if m := re.match(r"^(\d+)\s+([-+0-9.EeDd]+)", line):
                ao = int(m.group(1)) - 1
                if 0 <= ao < nao:
                    current.coefficients[ao] = float(
                        m.group(2).replace("D", "E").replace("d", "e"))
    return data


# Grid spacing in Bohr. The isosurface is built by marching cubes over these
# samples, so the spacing is what decides whether a lobe looks smooth or
# faceted; 0.12 is fine enough that the facets fall below a pixel at normal
# zoom. Large molecules relax towards MAX_SPACING to stay interactive.
TARGET_SPACING = 0.12
MAX_SPACING = 0.30


def _grid(data: MoldenData, max_points: int, margin: float = 3.4,
          target: float = TARGET_SPACING):
    """Regular grid (in Bohr) that encloses the molecule with a margin."""
    coords = np.array([[a[1], a[2], a[3]] for a in data.atoms]) / BOHR_TO_ANGSTROM
    lo = coords.min(axis=0) - margin
    hi = coords.max(axis=0) + margin
    span = hi - lo
    # As fine as TARGET_SPACING, unless that would blow the point budget.
    budget = float((span.prod() / max_points) ** (1 / 3))
    spacing = min(max(target, budget), MAX_SPACING)
    counts = np.maximum((span / spacing).astype(int) + 1, 2)
    axes = tuple(np.linspace(lo[d], lo[d] + spacing * (counts[d] - 1), counts[d])
                 for d in range(3))
    return coords, lo, spacing, counts, axes


def _point_budget(data: MoldenData, max_points: int) -> int:
    """Trade grid points against basis size to keep a redraw interactive.

    Evaluation cost is points x basis functions, so a large molecule gets a
    coarser grid — which is also when a coarser grid is least visible, since
    the whole molecule is drawn smaller.
    """
    WORK = 8e7
    return int(max(120_000, min(max_points, WORK / max(1, len(data.basis)))))


def _ao_norms(data: MoldenData) -> np.ndarray:
    """1/sqrt(<mu|mu>) per basis function.

    Molden files do not agree on whether the contraction itself is
    normalized, so every AO is scaled to unit norm here. That is the
    convention the MO coefficients follow, which is what makes the density
    integrate to the electron count.
    """
    norms = np.ones(len(data.basis))
    for index, bf in enumerate(data.basis):
        total = 0.0
        for ei, ci in zip(bf.exponents, bf.coefficients):
            for ej, cj in zip(bf.exponents, bf.coefficients):
                gamma = ei + ej
                term = 1.0
                for d in range(3):
                    term *= _overlap_1d(bf.powers[d], bf.powers[d], 0.0, 0.0, gamma)
                total += ci * cj * term
        total *= bf.scale * bf.scale
        norms[index] = 1.0 / math.sqrt(total) if total > 1e-30 else 1.0
    return norms


def _ao_slab(data: MoldenData, coords: np.ndarray, norms: np.ndarray,
             xs: np.ndarray, ys: np.ndarray, zs: np.ndarray) -> np.ndarray:
    """Every basis function on one slab, shaped (n_basis, n_points).

    A Cartesian Gaussian is separable — exp(-a r^2) is the product of three
    one-dimensional factors — so the exponentials are evaluated along each
    axis once and combined by outer products. That replaces one exp() per
    grid point per primitive with three cheap multiplies, which is what
    makes a grid fine enough to look smooth affordable.
    """
    out = np.empty((len(data.basis), xs.size * ys.size * zs.size), dtype=np.float32)
    for index, bf in enumerate(data.basis):
        ax, ay, az = coords[bf.atom_index]
        dx, dy, dz = xs - ax, ys - ay, zs - az
        px, py, pz = bf.powers
        block = np.zeros((xs.size, ys.size, zs.size))
        for exp_a, c in zip(bf.exponents, bf.coefficients):
            fx = (dx**px) * np.exp(-exp_a * dx * dx)
            fy = (dy**py) * np.exp(-exp_a * dy * dy)
            fz = (dz**pz) * np.exp(-exp_a * dz * dz)
            block += c * (fx[:, None, None] * fy[None, :, None] * fz[None, None, :])
        out[index] = (bf.scale * norms[index] * block).ravel()
    return out


def _slabs(count: int, basis_size: int, points_per_plane: int):
    """Split the slow axis so one slab of AO values stays a sane size."""
    per_slab = max(1, int(3e7 / max(1, basis_size * max(1, points_per_plane))))
    for start in range(0, count, per_slab):
        yield start, min(start + per_slab, count)


def _evaluate(data: MoldenData, coords, axes, weights: np.ndarray) -> np.ndarray:
    """Evaluate sum_mu w[k, mu] phi_mu(r) on the grid, for every row k.

    One matrix product per slab keeps this in BLAS rather than in Python,
    which matters once a molecule has a few hundred basis functions.
    """
    xs, ys, zs = axes
    norms = _ao_norms(data)
    result = np.zeros((weights.shape[0], xs.size, ys.size, zs.size), dtype=np.float32)
    single = weights.astype(np.float32)
    plane = ys.size * zs.size
    for start, stop in _slabs(xs.size, len(data.basis), plane):
        aos = _ao_slab(data, coords, norms, xs[start:stop], ys, zs)
        result[:, start:stop] = (single @ aos).reshape(
            weights.shape[0], stop - start, ys.size, zs.size)
    return result


def _cube_text(data: MoldenData, title: str, subtitle: str, coords, lo,
               spacing, counts, values: np.ndarray) -> str:
    out = [
        title,
        subtitle,
        f"{len(data.atoms):5d} {lo[0]:12.6f} {lo[1]:12.6f} {lo[2]:12.6f}",
        f"{counts[0]:5d} {spacing:12.6f} {0.0:12.6f} {0.0:12.6f}",
        f"{counts[1]:5d} {0.0:12.6f} {spacing:12.6f} {0.0:12.6f}",
        f"{counts[2]:5d} {0.0:12.6f} {0.0:12.6f} {spacing:12.6f}",
    ]
    for (symbol, *_), (x, y, z) in zip(data.atoms, coords):
        z_num = ATOMIC_NUMBER.get(symbol, 0)
        out.append(f"{z_num:5d} {float(z_num):12.6f} {x:12.6f} {y:12.6f} {z:12.6f}")
    # Gaussian cube format: z runs fastest.
    flat = values.reshape(counts[0] * counts[1], counts[2])
    for row in flat:
        for i in range(0, len(row), 6):
            out.append(" ".join(f"{v:13.5E}" for v in row[i : i + 6]))
    return "\n".join(out) + "\n"


def orbital_cube(data: MoldenData, mo_index: int, max_points: int = 1_100_000) -> str:
    """Gaussian cube text of one MO evaluated on a regular grid (Bohr)."""
    orbital = data.orbitals[mo_index - 1]
    coords, lo, spacing, counts, axes = _grid(data, _point_budget(data, max_points))
    values = _evaluate(data, coords, axes,
                       np.asarray(orbital.coefficients, dtype=float)[None, :])[0]
    return _cube_text(
        data,
        f"OQP Studio molecular orbital {orbital.index}",
        f"MO {orbital.index}  E={orbital.energy:.6f}  occ={orbital.occupancy}",
        coords, lo, spacing, counts, values,
    )


def orbital_coeff_cube(data: MoldenData, coefficients: np.ndarray,
                       title: str, subtitle: str,
                       max_points: int = 1_100_000) -> str:
    """Gaussian cube for an AO-coefficient vector not stored as a Molden MO."""
    vector = np.asarray(coefficients, dtype=float).ravel()
    if vector.size != len(data.basis):
        raise ValueError("orbital coefficient count does not match the Molden basis")
    coords, lo, spacing, counts, axes = _grid(data, _point_budget(data, max_points))
    values = _evaluate(data, coords, axes, vector[None, :])[0]
    return _cube_text(data, title, subtitle, coords, lo, spacing, counts, values)


def matrix_density_cube(data: MoldenData, density_ao: np.ndarray,
                        title: str, subtitle: str,
                        max_points: int = 700_000) -> str:
    """Gaussian cube for a real AO density matrix, including signed fields."""
    matrix = np.asarray(density_ao, dtype=float)
    nbf = len(data.basis)
    if matrix.shape != (nbf, nbf):
        raise ValueError("density matrix dimensions do not match the Molden basis")
    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues, vectors = np.linalg.eigh(matrix)
    largest = float(np.max(np.abs(eigenvalues))) if eigenvalues.size else 0.0
    keep = np.abs(eigenvalues) > max(largest * 1e-12, 1e-14)
    coords, lo, spacing, counts, axes = _grid(data, _point_budget(data, max_points))
    if np.any(keep):
        orbitals = _evaluate(data, coords, axes, vectors[:, keep].T)
        values = np.einsum("i,ixyz->xyz", eigenvalues[keep], orbitals * orbitals)
    else:
        values = np.zeros(tuple(int(value) for value in counts))
    return _cube_text(data, title, subtitle, coords, lo, spacing, counts, values)


def _occupancies(data: MoldenData) -> list[float]:
    """Occupation numbers, defaulted for files that omit them."""
    spins = {o.spin for o in data.orbitals}
    default = 1.0 if len(spins) > 1 else 2.0
    return [default if o.occupancy is None else float(o.occupancy)
            for o in data.orbitals]


def density_grid(data: MoldenData, coords, axes, spin: str = "total") -> np.ndarray:
    """Electron density (or alpha-beta spin density) on the grid.

    rho(r) = sum_i n_i |psi_i(r)|^2 over the occupied orbitals, which is exact
    for a single-determinant reference and is what a Molden file describes.
    """
    unrestricted = len({o.spin.lower() for o in data.orbitals}) > 1
    rows, weights = [], []
    for orbital, occupancy in zip(data.orbitals, _occupancies(data)):
        if occupancy <= 1e-8:
            continue
        alpha_orbital = orbital.spin.lower().startswith("a")
        weight = occupancy
        if spin == "spin":
            if unrestricted:
                weight = occupancy if alpha_orbital else -occupancy
            else:
                # A restricted file carries one orbital per spatial function,
                # so its spin density is n_alpha - n_beta = min(n,1)-max(n-1,0):
                # zero for a doubly occupied orbital, one for a singly one.
                weight = min(occupancy, 1.0) - max(occupancy - 1.0, 0.0)
        elif spin == "alpha":
            weight = occupancy if (unrestricted and alpha_orbital) else (
                min(occupancy, 1.0) if not unrestricted else 0.0)
        elif spin == "beta":
            weight = occupancy if (unrestricted and not alpha_orbital) else (
                max(occupancy - 1.0, 0.0) if not unrestricted else 0.0)
        if abs(weight) < 1e-10:
            continue
        rows.append(np.asarray(orbital.coefficients, dtype=float))
        weights.append(weight)
    if not rows:
        return np.zeros((axes[0].size, axes[1].size, axes[2].size))
    orbitals = _evaluate(data, coords, axes, np.vstack(rows))
    density = np.zeros(orbitals.shape[1:])
    for psi, weight in zip(orbitals, weights):
        density += weight * psi * psi
    return density


def density_cube(data: MoldenData, spin: str = "total",
                 max_points: int = 700_000) -> str:
    coords, lo, spacing, counts, axes = _grid(data, _point_budget(data, max_points))
    values = density_grid(data, coords, axes, spin)
    label = {"total": "electron density", "spin": "spin density",
             "alpha": "alpha density", "beta": "beta density"}[spin]
    return _cube_text(data, f"OQP Studio {label}",
                      f"{label} (e/bohr^3)", coords, lo, spacing, counts, values)


def mulliken_charges(data: MoldenData) -> list[float]:
    """Mulliken atomic charges from the Molden orbitals.

    The AO overlap is built analytically from the Gaussian basis, so this
    works from a Molden file alone, without the run's density matrix.
    """
    overlap = ao_overlap(data)
    natoms = len(data.atoms)
    populations = np.zeros(natoms)
    for orbital, occupancy in zip(data.orbitals, _occupancies(data)):
        if occupancy <= 1e-8:
            continue
        c = np.asarray(orbital.coefficients, dtype=float)
        # Mulliken partition: gross population of AO mu is n * c_mu (S c)_mu.
        gross = occupancy * c * (overlap @ c)
        for index, bf in enumerate(data.basis):
            populations[bf.atom_index] += gross[index]
    charges = []
    for index, (symbol, *_) in enumerate(data.atoms):
        charges.append(float(ATOMIC_NUMBER.get(symbol, 0) - populations[index]))
    return charges


def _binomial_prefactor(s: int, l1: int, l2: int, pa: float, pb: float) -> float:
    total = 0.0
    for t in range(s + 1):
        u = s - t
        if t > l1 or u > l2:
            continue
        total += (math.comb(l1, t) * math.comb(l2, u)
                  * (pa ** (l1 - t)) * (pb ** (l2 - u)))
    return total


def _overlap_1d(l1: int, l2: int, pa: float, pb: float, gamma: float) -> float:
    total = 0.0
    for i in range(1 + (l1 + l2) // 2):
        total += (_binomial_prefactor(2 * i, l1, l2, pa, pb)
                  * _double_factorial(2 * i - 1) / (2 * gamma) ** i)
    return total * math.sqrt(math.pi / gamma)


def ao_overlap(data: MoldenData) -> np.ndarray:
    """Analytic overlap matrix of the Cartesian Gaussian basis.

    Molden files disagree about whether contractions are normalized, so the
    matrix is scaled to a unit diagonal — the convention Mulliken analysis
    assumes anyway.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in data.atoms]) / BOHR_TO_ANGSTROM
    n = len(data.basis)
    matrix = np.zeros((n, n))
    for i, bi in enumerate(data.basis):
        ai = coords[bi.atom_index]
        for j in range(i, n):
            bj = data.basis[j]
            aj = coords[bj.atom_index]
            ab2 = float(np.sum((ai - aj) ** 2))
            total = 0.0
            for ei, ci in zip(bi.exponents, bi.coefficients):
                for ej, cj in zip(bj.exponents, bj.coefficients):
                    gamma = ei + ej
                    p = (ei * ai + ej * aj) / gamma
                    pre = math.exp(-ei * ej * ab2 / gamma)
                    term = pre
                    for d in range(3):
                        term *= _overlap_1d(bi.powers[d], bj.powers[d],
                                            p[d] - ai[d], p[d] - aj[d], gamma)
                    total += ci * cj * term
            matrix[i, j] = matrix[j, i] = total * bi.scale * bj.scale
    norms = _ao_norms(data)
    return matrix * np.outer(norms, norms)


def esp_cube(data: MoldenData, charges: list[float] | None = None,
             max_points: int = 700_000) -> str:
    """Molecular electrostatic potential from atomic point charges.

    V(r) = sum_A q_A / |r - R_A|, with q_A taken from the run's own charges
    when the job exported them and computed by Mulliken analysis otherwise.
    This is the charge-model MEP that most viewers draw; it is not the exact
    density integral, and the cube header says so.
    """
    # The potential falls off as 1/r, so its isosurfaces reach much further
    # from the nuclei than an orbital's do; a tight box would cut them and
    # leave a ring where the surface meets the box face. It is also a cheap
    # field to evaluate, so the larger box costs almost nothing.
    coords, lo, spacing, counts, axes = _grid(data, max_points, margin=7.0,
                                              target=0.18)
    if charges is None or len(charges) != len(data.atoms):
        charges = mulliken_charges(data)
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    values = np.zeros(gx.shape)
    for (x, y, z), q in zip(coords, charges):
        r = np.sqrt((gx - x) ** 2 + (gy - y) ** 2 + (gz - z) ** 2)
        np.clip(r, 0.3, None, out=r)     # keep the nuclear cusp finite
        values += q / r
    return _cube_text(data, "OQP Studio electrostatic potential",
                      "MEP from atomic point charges (Hartree/e)",
                      coords, lo, spacing, counts, values)


def atoms_to_xyz(data: MoldenData, title: str = "geometry") -> str:
    lines = [str(len(data.atoms)), title]
    for symbol, x, y, z in data.atoms:
        lines.append(f"{symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


@dataclass
class NormalMode:
    index: int
    frequency: float  # cm^-1
    displacements: np.ndarray  # (natoms, 3) in Bohr


@dataclass
class VibrationData:
    atoms: list[tuple[str, float, float, float]]  # Angstrom
    modes: list[NormalMode] = field(default_factory=list)
    intensities: list[float] = field(default_factory=list)


def parse_vibrations(text: str) -> VibrationData:
    """Frequencies and normal modes from [FREQ]/[FR-COORD]/[FR-NORM-COORD]."""
    lines = text.splitlines()

    def section(header: str) -> list[str]:
        start = next(
            (i for i, ln in enumerate(lines) if ln.strip().upper() == header), -1)
        if start < 0:
            return []
        out = []
        for ln in lines[start + 1 :]:
            if re.match(r"^\s*\[", ln):
                break
            out.append(ln)
        return out

    frequencies = [float(ln) for ln in section("[FREQ]") if ln.strip()]
    intensities = [float(ln) for ln in section("[INT]") if ln.strip()]

    atoms: list[tuple[str, float, float, float]] = []
    for ln in section("[FR-COORD]"):
        parts = ln.split()
        if len(parts) >= 4:
            atoms.append((parts[0], *(float(v) * BOHR_TO_ANGSTROM for v in parts[1:4])))

    modes: list[NormalMode] = []
    current: list[list[float]] | None = None
    for ln in section("[FR-NORM-COORD]"):
        if re.match(r"^\s*vibration", ln, re.IGNORECASE):
            if current:
                modes.append(NormalMode(len(modes) + 1, 0.0, np.asarray(current)))
            current = []
            continue
        parts = ln.split()
        if current is not None and len(parts) >= 3:
            current.append([float(v) for v in parts[:3]])
    if current:
        modes.append(NormalMode(len(modes) + 1, 0.0, np.asarray(current)))
    for mode, freq in zip(modes, frequencies):
        mode.frequency = freq
    return VibrationData(atoms=atoms, modes=modes, intensities=intensities)


def mode_trajectory(vib: VibrationData, mode_index: int,
                    frames: int = 24, amplitude: float = 0.6) -> str:
    """Multi-frame XYZ oscillating along one normal mode (for animation)."""
    mode = vib.modes[mode_index - 1]
    base = np.array([[a[1], a[2], a[3]] for a in vib.atoms])
    disp = mode.displacements * BOHR_TO_ANGSTROM
    out: list[str] = []
    for k in range(frames):
        factor = amplitude * math.sin(2 * math.pi * k / frames)
        coords = base + factor * disp
        out.append(str(len(vib.atoms)))
        out.append(f"mode {mode.index}  {mode.frequency:.2f} cm-1  frame {k + 1}/{frames}")
        for (symbol, *_), (x, y, z) in zip(vib.atoms, coords):
            out.append(f"{symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(out) + "\n"
