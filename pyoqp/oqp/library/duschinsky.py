"""Duschinsky transformation from two Cartesian geometries and Hessians.

The implementation is deliberately independent of the compiled OpenQP library.  It
uses the mass-weighted Cartesian convention

``x = M**(-1/2) L Q``

where the columns of ``L`` are orthonormal vibrational vectors.  After a
mass-weighted proper Kabsch rotation of the excited-state geometry into the
ground-state frame, the returned transformation is

``Q_excited = J @ Q_ground + K``

with ``J = L_excited.T @ L_ground`` and
``K = L_excited.T @ M**(1/2) (R_ground - R_excited)``.

Translations and rotations are removed independently at each equilibrium
geometry.  Consequently, a genuine geometry change can make the two vibrational
subspaces slightly different.  The code reports that physically meaningful
departure through :class:`OrthogonalityDiagnostics`; it never silently replaces
``J`` by a nearest orthogonal matrix.

The module stops at the harmonic-coordinate transformation.  It does not
evaluate Franck--Condon factors, Herzberg--Teller terms, transition moments, or
spectral intensities.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class NormalModeSet:
    """Projected normal modes at one equilibrium geometry.

    ``mass_weighted_modes`` and ``cartesian_modes`` store modes as columns.
    The former obey ``L.T @ L = I`` and the latter obey
    ``C.T @ M @ C = I``.  ``force_constants`` are the eigenvalues of the
    projected mass-weighted Hessian; negative values are retained.
    """

    geometry: FloatArray
    masses: FloatArray
    force_constants: FloatArray
    mass_weighted_modes: FloatArray
    cartesian_modes: FloatArray
    external_basis: FloatArray
    inertia_moments: FloatArray
    external_rank: int
    is_linear: bool
    hessian_symmetry_residual: float
    external_contamination: float


@dataclass(frozen=True)
class ModeAlignment:
    """Gauge transformation applied to the excited-state normal modes."""

    modes: FloatArray
    transformation: FloatArray
    overlap_with_ground: FloatArray
    aligned_degenerate_blocks: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class OrthogonalityDiagnostics:
    """Left and right orthogonality residuals of a Duschinsky matrix."""

    left_max_abs: float
    right_max_abs: float
    left_frobenius: float
    right_frobenius: float
    singular_values: FloatArray
    determinant: float | None

    @property
    def max_abs(self) -> float:
        """Largest elementwise residual from either orthogonality identity."""

        return max(self.left_max_abs, self.right_max_abs)


@dataclass(frozen=True)
class DuschinskyResult:
    """Duschinsky rotation, displacement, and alignment diagnostics."""

    J: FloatArray
    K: FloatArray
    ground_modes: NormalModeSet
    excited_modes: NormalModeSet
    mode_alignment: ModeAlignment
    rotation_excited_to_ground: FloatArray
    translation_excited_to_ground: FloatArray
    aligned_excited_geometry: FloatArray
    mass_weighted_rmsd: float
    eckart_residual: FloatArray
    orthogonality: OrthogonalityDiagnostics


def _finite_array(name: str, values: ArrayLike) -> FloatArray:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _validate_geometry(name: str, geometry: ArrayLike) -> FloatArray:
    array = _finite_array(name, geometry)
    if array.ndim != 2 or array.shape[1] != 3 or array.shape[0] == 0:
        raise ValueError(f"{name} must have shape (N, 3), got {array.shape}")
    return array.copy()


def _validate_masses(
    masses: ArrayLike,
    natom: int,
    *,
    name: str = "masses",
) -> FloatArray:
    array = _finite_array(name, masses)
    if array.shape != (natom,):
        raise ValueError(f"{name} must have shape ({natom},), got {array.shape}")
    if np.any(array <= 0.0):
        raise ValueError(f"all {name} must be positive")
    return array.copy()


def _validate_atoms(name: str, atoms: ArrayLike, natom: int) -> np.ndarray:
    array = np.asarray(atoms)
    if array.shape != (natom,):
        raise ValueError(f"{name} must have shape ({natom},), got {array.shape}")
    if array.dtype.kind in "fc" and not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array.copy()


def _validate_state_correspondence(
    ground_masses: ArrayLike,
    excited_masses: ArrayLike,
    natom: int,
    *,
    ground_atoms: ArrayLike | None,
    excited_atoms: ArrayLike | None,
    mass_relative_tolerance: float,
    mass_absolute_tolerance: float,
) -> FloatArray:
    if mass_relative_tolerance < 0.0 or mass_absolute_tolerance < 0.0:
        raise ValueError("mass tolerances must be non-negative")

    ground_mass = _validate_masses(
        ground_masses, natom, name="ground_masses"
    )
    excited_mass = _validate_masses(
        excited_masses, natom, name="excited_masses"
    )
    if not np.allclose(
        ground_mass,
        excited_mass,
        rtol=mass_relative_tolerance,
        atol=mass_absolute_tolerance,
    ):
        difference = float(np.max(np.abs(ground_mass - excited_mass)))
        raise ValueError(
            "atomic mass mismatch between ground and excited states: "
            f"maximum absolute difference is {difference:.3e}"
        )

    if (ground_atoms is None) != (excited_atoms is None):
        raise ValueError(
            "ground_atoms and excited_atoms must either both be given or both "
            "be omitted"
        )
    if ground_atoms is not None:
        ground_atom_array = _validate_atoms(
            "ground_atoms", ground_atoms, natom
        )
        excited_atom_array = _validate_atoms(
            "excited_atoms", excited_atoms, natom
        )
        if not np.array_equal(ground_atom_array, excited_atom_array):
            raise ValueError(
                "atom identity or order mismatch between ground and excited states"
            )

    return ground_mass


def _validate_hessian(
    name: str,
    hessian: ArrayLike,
    ncoord: int,
    symmetry_tolerance: float,
) -> tuple[FloatArray, float]:
    array = _finite_array(name, hessian)
    if array.shape != (ncoord, ncoord):
        raise ValueError(
            f"{name} must have shape ({ncoord}, {ncoord}), got {array.shape}"
        )
    if symmetry_tolerance < 0.0:
        raise ValueError("hessian_symmetry_tolerance must be non-negative")
    absolute_residual = float(np.max(np.abs(array - array.T), initial=0.0))
    scale = max(1.0, float(np.max(np.abs(array), initial=0.0)))
    relative_residual = absolute_residual / scale
    if relative_residual > symmetry_tolerance:
        raise ValueError(
            f"{name} is not symmetric: relative maximum residual "
            f"{relative_residual:.3e} exceeds {symmetry_tolerance:.3e}"
        )
    return 0.5 * (array + array.T), relative_residual


def _center_of_mass(geometry: FloatArray, masses: FloatArray) -> FloatArray:
    return np.sum(geometry * masses[:, None], axis=0) / np.sum(masses)


def _inertia_tensor(centered: FloatArray, masses: FloatArray) -> FloatArray:
    tensor = np.zeros((3, 3), dtype=float)
    identity = np.eye(3)
    for mass, position in zip(masses, centered, strict=True):
        tensor += mass * (
            np.dot(position, position) * identity - np.outer(position, position)
        )
    return tensor


def _external_and_vibrational_bases(
    geometry: FloatArray,
    masses: FloatArray,
    rank_tolerance: float,
) -> tuple[FloatArray, FloatArray, FloatArray, int, bool]:
    if rank_tolerance <= 0.0:
        raise ValueError("external_rank_tolerance must be positive")

    natom = geometry.shape[0]
    ncoord = 3 * natom
    center = _center_of_mass(geometry, masses)
    centered = geometry - center
    sqrt_mass = np.sqrt(masses)

    candidates = np.zeros((ncoord, 6), dtype=float)
    for axis in range(3):
        candidates[axis::3, axis] = sqrt_mass

    axes = np.eye(3)
    for atom, (position, root_mass) in enumerate(
        zip(centered, sqrt_mass, strict=True)
    ):
        for axis in range(3):
            candidates[3 * atom : 3 * atom + 3, 3 + axis] = (
                root_mass * np.cross(axes[axis], position)
            )

    # Column normalization prevents molecular size from deciding the numerical
    # rank when translations and rotations have very different dimensions.
    norms = np.linalg.norm(candidates, axis=0)
    nonzero = np.zeros(6, dtype=bool)
    nonzero[:3] = True
    rotational_scale = float(np.max(norms[3:], initial=0.0))
    if rotational_scale > 0.0:
        nonzero[3:] = norms[3:] > rank_tolerance * rotational_scale
    normalized = candidates[:, nonzero] / norms[nonzero]
    u, singular_values, _ = np.linalg.svd(normalized, full_matrices=True)
    if singular_values.size:
        rank = int(
            np.count_nonzero(
                singular_values > rank_tolerance * float(singular_values[0])
            )
        )
    else:
        rank = 0

    if natom == 1:
        expected_ranks = (3,)
    else:
        expected_ranks = (5, 6)
    if rank not in expected_ranks:
        raise ValueError(
            "geometry does not define a noncoincident linear or nonlinear molecule: "
            f"external rank is {rank}"
        )

    external = u[:, :rank]
    vibrational = u[:, rank:]
    inertia = np.linalg.eigvalsh(_inertia_tensor(centered, masses))
    is_linear = natom > 1 and rank == 5
    return external, vibrational, inertia, rank, is_linear


def projected_normal_modes(
    geometry: ArrayLike,
    masses: ArrayLike,
    hessian: ArrayLike,
    *,
    hessian_symmetry_tolerance: float = 1.0e-10,
    external_rank_tolerance: float = 1.0e-10,
) -> NormalModeSet:
    """Return physical normal modes after translation/rotation projection.

    Parameters use any self-consistent Cartesian length, energy, and mass units.
    The force constants and the displacement vector returned later inherit those
    choices.  For OpenQP output the usual convention is bohr, hartree/bohr**2,
    and atomic mass units.
    """

    xyz = _validate_geometry("geometry", geometry)
    mass = _validate_masses(masses, xyz.shape[0])
    ncoord = 3 * xyz.shape[0]
    hess, symmetry_residual = _validate_hessian(
        "hessian", hessian, ncoord, hessian_symmetry_tolerance
    )

    external, vibrational, inertia, rank, is_linear = (
        _external_and_vibrational_bases(xyz, mass, external_rank_tolerance)
    )
    coordinate_masses = np.repeat(mass, 3)
    inverse_root_mass = 1.0 / np.sqrt(coordinate_masses)
    mass_weighted_hessian = (
        inverse_root_mass[:, None] * hess * inverse_root_mass[None, :]
    )

    reduced_hessian = vibrational.T @ mass_weighted_hessian @ vibrational
    force_constants, reduced_modes = np.linalg.eigh(
        0.5 * (reduced_hessian + reduced_hessian.T)
    )
    modes = vibrational @ reduced_modes
    cartesian_modes = inverse_root_mass[:, None] * modes
    contamination = float(
        np.max(np.abs(external.T @ modes), initial=0.0)
    )

    return NormalModeSet(
        geometry=xyz,
        masses=mass,
        force_constants=force_constants,
        mass_weighted_modes=modes,
        cartesian_modes=cartesian_modes,
        external_basis=external,
        inertia_moments=inertia,
        external_rank=rank,
        is_linear=is_linear,
        hessian_symmetry_residual=symmetry_residual,
        external_contamination=contamination,
    )


def mass_weighted_kabsch(
    reference_geometry: ArrayLike,
    moving_geometry: ArrayLike,
    masses: ArrayLike,
) -> tuple[FloatArray, FloatArray, FloatArray, float, FloatArray]:
    """Align ``moving_geometry`` to ``reference_geometry`` with a proper rotation.

    The returned rotation follows a row-vector convention:
    ``aligned = moving @ rotation + translation``.  The final vector is the
    rotational Eckart residual ``sum_i m_i (r_moving_i x r_reference_i)`` in
    the common center-of-mass frame.
    """

    reference = _validate_geometry("reference_geometry", reference_geometry)
    moving = _validate_geometry("moving_geometry", moving_geometry)
    if moving.shape != reference.shape:
        raise ValueError(
            "reference_geometry and moving_geometry must have the same shape"
        )
    mass = _validate_masses(masses, reference.shape[0])

    reference_com = _center_of_mass(reference, mass)
    moving_com = _center_of_mass(moving, mass)
    reference_centered = reference - reference_com
    moving_centered = moving - moving_com

    covariance = moving_centered.T @ (mass[:, None] * reference_centered)
    left, _, right_transpose = np.linalg.svd(covariance)
    handedness = np.linalg.det(left @ right_transpose)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if handedness >= 0.0 else -1.0
    rotation = left @ correction @ right_transpose
    translation = reference_com - moving_com @ rotation
    aligned = moving @ rotation + translation

    difference = aligned - reference
    rmsd = float(
        np.sqrt(np.sum(mass[:, None] * difference**2) / np.sum(mass))
    )
    aligned_centered = aligned - reference_com
    eckart = np.sum(
        mass[:, None] * np.cross(aligned_centered, reference_centered), axis=0
    )
    return rotation, translation, aligned, rmsd, eckart


def rotate_cartesian_hessian(
    hessian: ArrayLike,
    rotation_excited_to_ground: ArrayLike,
    natom: int,
) -> FloatArray:
    """Rotate a Cartesian Hessian into the aligned coordinate frame."""

    if natom <= 0:
        raise ValueError("natom must be positive")
    rotation = _finite_array(
        "rotation_excited_to_ground", rotation_excited_to_ground
    )
    if rotation.shape != (3, 3):
        raise ValueError(
            "rotation_excited_to_ground must have shape (3, 3), "
            f"got {rotation.shape}"
        )
    rotation_residual = np.max(np.abs(rotation.T @ rotation - np.eye(3)))
    if rotation_residual > 1.0e-10 or np.linalg.det(rotation) <= 0.0:
        raise ValueError("rotation_excited_to_ground must be a proper rotation")

    matrix = _finite_array("hessian", hessian)
    ncoord = 3 * natom
    if matrix.shape != (ncoord, ncoord):
        raise ValueError(
            f"hessian must have shape ({ncoord}, {ncoord}), got {matrix.shape}"
        )
    cartesian_rotation = np.kron(np.eye(natom), rotation.T)
    return cartesian_rotation @ matrix @ cartesian_rotation.T


def _degenerate_blocks(
    values: FloatArray,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> tuple[tuple[int, int], ...]:
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("degeneracy tolerances must be non-negative")
    if values.size == 0:
        return ()
    blocks: list[tuple[int, int]] = []
    start = 0
    for index in range(1, values.size):
        scale = max(abs(float(values[index - 1])), abs(float(values[index])))
        threshold = absolute_tolerance + relative_tolerance * scale
        if abs(float(values[index] - values[index - 1])) > threshold:
            blocks.append((start, index))
            start = index
    blocks.append((start, values.size))
    return tuple(blocks)


def align_mode_gauge(
    ground_modes: ArrayLike,
    excited_modes: ArrayLike,
    ground_force_constants: ArrayLike,
    excited_force_constants: ArrayLike,
    *,
    degeneracy_relative_tolerance: float = 1.0e-8,
    degeneracy_absolute_tolerance: float = 1.0e-12,
) -> ModeAlignment:
    """Fix mode phases and align common degenerate subspaces.

    Nondegenerate modes receive only a sign change.  A Procrustes rotation is
    applied only when the same contiguous index block is degenerate in both
    states.  General nondegenerate Duschinsky mixing is therefore preserved.
    """

    ground = _finite_array("ground_modes", ground_modes)
    excited = _finite_array("excited_modes", excited_modes)
    if ground.ndim != 2 or excited.shape != ground.shape:
        raise ValueError(
            "ground_modes and excited_modes must have the same two-dimensional shape"
        )
    nmode = ground.shape[1]
    ground_values = _finite_array(
        "ground_force_constants", ground_force_constants
    )
    excited_values = _finite_array(
        "excited_force_constants", excited_force_constants
    )
    if ground_values.shape != (nmode,) or excited_values.shape != (nmode,):
        raise ValueError(f"force constants must each have shape ({nmode},)")

    identity = np.eye(nmode)
    if np.max(np.abs(ground.T @ ground - identity), initial=0.0) > 1.0e-9:
        raise ValueError("ground_modes must have orthonormal columns")
    if np.max(np.abs(excited.T @ excited - identity), initial=0.0) > 1.0e-9:
        raise ValueError("excited_modes must have orthonormal columns")

    ground_blocks = set(
        _degenerate_blocks(
            ground_values,
            degeneracy_relative_tolerance,
            degeneracy_absolute_tolerance,
        )
    )
    excited_blocks = set(
        _degenerate_blocks(
            excited_values,
            degeneracy_relative_tolerance,
            degeneracy_absolute_tolerance,
        )
    )
    common_degenerate = tuple(
        sorted(
            block
            for block in ground_blocks & excited_blocks
            if block[1] - block[0] > 1
        )
    )

    transformation = np.eye(nmode)
    covered: set[int] = set()
    for start, stop in common_degenerate:
        cross_overlap = excited[:, start:stop].T @ ground[:, start:stop]
        left, _, right_transpose = np.linalg.svd(cross_overlap)
        transformation[start:stop, start:stop] = left @ right_transpose
        covered.update(range(start, stop))

    aligned = excited @ transformation
    for index in range(nmode):
        if index in covered:
            continue
        overlaps = aligned[:, index].T @ ground
        anchor = int(np.argmax(np.abs(overlaps)))
        if float(overlaps[anchor]) < 0.0:
            transformation[index, index] *= -1.0
            aligned[:, index] *= -1.0

    return ModeAlignment(
        modes=aligned,
        transformation=transformation,
        overlap_with_ground=aligned.T @ ground,
        aligned_degenerate_blocks=common_degenerate,
    )


def orthogonality_diagnostics(matrix: ArrayLike) -> OrthogonalityDiagnostics:
    """Return two-sided orthogonality diagnostics for a rectangular matrix."""

    array = _finite_array("matrix", matrix)
    if array.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    left = array @ array.T - np.eye(array.shape[0])
    right = array.T @ array - np.eye(array.shape[1])
    determinant = (
        float(np.linalg.det(array))
        if array.shape[0] == array.shape[1]
        else None
    )
    return OrthogonalityDiagnostics(
        left_max_abs=float(np.max(np.abs(left), initial=0.0)),
        right_max_abs=float(np.max(np.abs(right), initial=0.0)),
        left_frobenius=float(np.linalg.norm(left)),
        right_frobenius=float(np.linalg.norm(right)),
        singular_values=np.linalg.svd(array, compute_uv=False),
        determinant=determinant,
    )


def compute_duschinsky(
    ground_geometry: ArrayLike,
    excited_geometry: ArrayLike,
    masses: ArrayLike,
    ground_hessian: ArrayLike,
    excited_hessian: ArrayLike,
    *,
    hessian_symmetry_tolerance: float = 1.0e-10,
    external_rank_tolerance: float = 1.0e-10,
    degeneracy_relative_tolerance: float = 1.0e-8,
    degeneracy_absolute_tolerance: float = 1.0e-12,
    orthogonality_tolerance: float | None = None,
) -> DuschinskyResult:
    """Compute the physical-mode Duschinsky matrix and displacement vector.

    Atom order and isotopic masses must be identical in the two states.  A
    linear-to-nonlinear change is rejected because the two harmonic coordinate
    spaces then have different dimensions.  ``K`` has units of
    ``sqrt(mass) * length``.
    """

    ground_xyz = _validate_geometry("ground_geometry", ground_geometry)
    excited_xyz = _validate_geometry("excited_geometry", excited_geometry)
    if excited_xyz.shape != ground_xyz.shape:
        raise ValueError("ground and excited geometries must have the same shape")
    mass = _validate_masses(masses, ground_xyz.shape[0])

    rotation, translation, aligned_excited, rmsd, eckart = mass_weighted_kabsch(
        ground_xyz, excited_xyz, mass
    )
    aligned_excited_hessian = rotate_cartesian_hessian(
        excited_hessian, rotation, ground_xyz.shape[0]
    )

    normal_mode_options = dict(
        hessian_symmetry_tolerance=hessian_symmetry_tolerance,
        external_rank_tolerance=external_rank_tolerance,
    )
    ground = projected_normal_modes(
        ground_xyz, mass, ground_hessian, **normal_mode_options
    )
    excited = projected_normal_modes(
        aligned_excited, mass, aligned_excited_hessian, **normal_mode_options
    )
    if ground.mass_weighted_modes.shape[1] != excited.mass_weighted_modes.shape[1]:
        raise ValueError(
            "ground and excited geometries have different physical-mode counts; "
            "a linear-to-nonlinear Duschinsky transformation is not supported"
        )

    alignment = align_mode_gauge(
        ground.mass_weighted_modes,
        excited.mass_weighted_modes,
        ground.force_constants,
        excited.force_constants,
        degeneracy_relative_tolerance=degeneracy_relative_tolerance,
        degeneracy_absolute_tolerance=degeneracy_absolute_tolerance,
    )
    coordinate_masses = np.repeat(mass, 3)
    aligned_cartesian_modes = (
        alignment.modes / np.sqrt(coordinate_masses)[:, None]
    )
    excited = replace(
        excited,
        mass_weighted_modes=alignment.modes,
        cartesian_modes=aligned_cartesian_modes,
        external_contamination=float(
            np.max(
                np.abs(excited.external_basis.T @ alignment.modes), initial=0.0
            )
        ),
    )

    J = excited.mass_weighted_modes.T @ ground.mass_weighted_modes
    mass_weighted_shift = np.sqrt(coordinate_masses) * (
        ground_xyz - aligned_excited
    ).reshape(-1)
    K = excited.mass_weighted_modes.T @ mass_weighted_shift
    diagnostics = orthogonality_diagnostics(J)
    if orthogonality_tolerance is not None:
        if orthogonality_tolerance < 0.0:
            raise ValueError("orthogonality_tolerance must be non-negative")
        if diagnostics.max_abs > orthogonality_tolerance:
            raise ValueError(
                "Duschinsky matrix orthogonality residual "
                f"{diagnostics.max_abs:.3e} exceeds "
                f"{orthogonality_tolerance:.3e}"
            )

    return DuschinskyResult(
        J=J,
        K=K,
        ground_modes=ground,
        excited_modes=excited,
        mode_alignment=alignment,
        rotation_excited_to_ground=rotation,
        translation_excited_to_ground=translation,
        aligned_excited_geometry=aligned_excited,
        mass_weighted_rmsd=rmsd,
        eckart_residual=eckart,
        orthogonality=diagnostics,
    )


def compute_duschinsky_from_arrays(
    ground_geometry: ArrayLike,
    ground_masses: ArrayLike,
    ground_hessian: ArrayLike,
    excited_geometry: ArrayLike,
    excited_masses: ArrayLike,
    excited_hessian: ArrayLike,
    *,
    ground_atoms: ArrayLike | None = None,
    excited_atoms: ArrayLike | None = None,
    mass_relative_tolerance: float = 1.0e-12,
    mass_absolute_tolerance: float = 1.0e-12,
    hessian_symmetry_tolerance: float = 1.0e-10,
    external_rank_tolerance: float = 1.0e-10,
    degeneracy_relative_tolerance: float = 1.0e-8,
    degeneracy_absolute_tolerance: float = 1.0e-12,
    orthogonality_tolerance: float | None = None,
) -> DuschinskyResult:
    """Compute a Duschinsky transformation from two explicit state records.

    Each state supplies its own geometry, isotopic masses, and Cartesian
    Hessian.  The mass arrays must agree within the requested tolerances.  If
    atomic identities are supplied, both arrays are required and their order
    must agree exactly.  Coordinates and Hessians must use one consistent unit
    system; ``K`` then has units of ``sqrt(mass) * length``.

    This helper returns normal modes, ``J``, ``K``, and diagnostics only.  It
    does not calculate Franck--Condon or Herzberg--Teller intensities.
    """

    ground_xyz = _validate_geometry("ground_geometry", ground_geometry)
    excited_xyz = _validate_geometry("excited_geometry", excited_geometry)
    if excited_xyz.shape != ground_xyz.shape:
        raise ValueError("ground and excited geometries must have the same shape")

    masses = _validate_state_correspondence(
        ground_masses,
        excited_masses,
        ground_xyz.shape[0],
        ground_atoms=ground_atoms,
        excited_atoms=excited_atoms,
        mass_relative_tolerance=mass_relative_tolerance,
        mass_absolute_tolerance=mass_absolute_tolerance,
    )
    return compute_duschinsky(
        ground_xyz,
        excited_xyz,
        masses,
        ground_hessian,
        excited_hessian,
        hessian_symmetry_tolerance=hessian_symmetry_tolerance,
        external_rank_tolerance=external_rank_tolerance,
        degeneracy_relative_tolerance=degeneracy_relative_tolerance,
        degeneracy_absolute_tolerance=degeneracy_absolute_tolerance,
        orthogonality_tolerance=orthogonality_tolerance,
    )


def _molecule_geometry(molecule: object, state_name: str) -> FloatArray:
    getter = getattr(molecule, "get_system", None)
    if not callable(getter):
        raise TypeError(f"{state_name}_molecule must provide get_system()")
    geometry = _finite_array(f"{state_name}_geometry", getter())
    if geometry.ndim == 1:
        if geometry.size == 0 or geometry.size % 3:
            raise ValueError(
                f"{state_name}_molecule.get_system() must contain 3N coordinates"
            )
        geometry = geometry.reshape((-1, 3))
    return _validate_geometry(f"{state_name}_geometry", geometry)


def _molecule_array(
    molecule: object,
    accessor: str,
    state_name: str,
) -> object:
    getter = getattr(molecule, accessor, None)
    if not callable(getter):
        raise TypeError(f"{state_name}_molecule must provide {accessor}()")
    return getter()


def _molecule_hessian(
    molecule: object,
    override: ArrayLike | None,
    state_name: str,
) -> ArrayLike:
    if override is not None:
        return override
    hessian = _molecule_array(molecule, "get_hess", state_name)
    if np.asarray(hessian).size == 0:
        raise ValueError(
            f"{state_name}_molecule contains no Hessian; pass "
            f"{state_name}_hessian explicitly"
        )
    return hessian


def compute_duschinsky_from_molecules(
    ground_molecule: object,
    excited_molecule: object,
    *,
    ground_hessian: ArrayLike | None = None,
    excited_hessian: ArrayLike | None = None,
    mass_relative_tolerance: float = 1.0e-12,
    mass_absolute_tolerance: float = 1.0e-12,
    hessian_symmetry_tolerance: float = 1.0e-10,
    external_rank_tolerance: float = 1.0e-10,
    degeneracy_relative_tolerance: float = 1.0e-8,
    degeneracy_absolute_tolerance: float = 1.0e-12,
    orthogonality_tolerance: float | None = None,
) -> DuschinskyResult:
    """Compute a Duschinsky transformation from two OpenQP molecules.

    Geometry, masses, and atomic identities are read through ``get_system()``,
    ``get_mass()``, and ``get_atoms()``.  By default, each Cartesian Hessian is
    read through ``get_hess()``; an explicit Hessian keyword can replace either
    stored value.  The two molecules must describe the same atoms, isotope
    masses, and atom order.

    A returned :class:`DuschinskyResult` is sufficient for subsequent
    vibrational-coordinate analysis but is not a Franck--Condon or
    Herzberg--Teller intensity calculation.
    """

    ground_geometry = _molecule_geometry(ground_molecule, "ground")
    excited_geometry = _molecule_geometry(excited_molecule, "excited")
    ground_masses = _molecule_array(ground_molecule, "get_mass", "ground")
    excited_masses = _molecule_array(excited_molecule, "get_mass", "excited")
    ground_atoms = _molecule_array(ground_molecule, "get_atoms", "ground")
    excited_atoms = _molecule_array(excited_molecule, "get_atoms", "excited")

    return compute_duschinsky_from_arrays(
        ground_geometry,
        ground_masses,
        _molecule_hessian(ground_molecule, ground_hessian, "ground"),
        excited_geometry,
        excited_masses,
        _molecule_hessian(excited_molecule, excited_hessian, "excited"),
        ground_atoms=ground_atoms,
        excited_atoms=excited_atoms,
        mass_relative_tolerance=mass_relative_tolerance,
        mass_absolute_tolerance=mass_absolute_tolerance,
        hessian_symmetry_tolerance=hessian_symmetry_tolerance,
        external_rank_tolerance=external_rank_tolerance,
        degeneracy_relative_tolerance=degeneracy_relative_tolerance,
        degeneracy_absolute_tolerance=degeneracy_absolute_tolerance,
        orthogonality_tolerance=orthogonality_tolerance,
    )


__all__ = (
    "DuschinskyResult",
    "ModeAlignment",
    "NormalModeSet",
    "OrthogonalityDiagnostics",
    "align_mode_gauge",
    "compute_duschinsky",
    "compute_duschinsky_from_arrays",
    "compute_duschinsky_from_molecules",
    "mass_weighted_kabsch",
    "orthogonality_diagnostics",
    "projected_normal_modes",
    "rotate_cartesian_hessian",
)
