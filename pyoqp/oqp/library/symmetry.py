"""Backend-free symmetry utilities for one-electron diagnostics.

This module currently provides metadata-only diagnostics used by the symmetry
planning gates. It does not change SCF/integral/response execution behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import numpy as np


@dataclass(frozen=True)
class OneElectronBlockLeakSummary:
    """Normalized summary payload for one-electron block leakage checks."""

    max_off_block_abs: float
    max_off_block_indices: list[tuple[int, int]]
    off_block_element_count: int
    within_tolerance: bool
    orthogonality_ok: bool
    orthogonality_max_deviation: float
    status: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "max_off_block_abs": float(self.max_off_block_abs),
            "max_off_block_indices": [list(pair) for pair in self.max_off_block_indices],
            "off_block_element_count": int(self.off_block_element_count),
            "within_tolerance": bool(self.within_tolerance),
            "transform_orthogonality": {
                "ok": bool(self.orthogonality_ok),
                "max_deviation": float(self.orthogonality_max_deviation),
            },
            "status": self.status,
            "use_integral_symmetry": False,
            "use_response_symmetry": False,
        }


def _to_float_array(values: Any) -> np.ndarray:
    return np.asarray(values, dtype=float)


def _as_list(values: Any) -> list[Any]:
    if values is None:
        return []
    if isinstance(values, list):
        return values
    return list(values)


def one_electron_block_diagnostics(
    matrix: Any,
    symmetry_adapted_transform: Any,
    basis_labels: Iterable[str],
    tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Compute off-block leakage diagnostics in a symmetry-adapted basis.

    Parameters
    ----------
    matrix:
        One-electron matrix (overlap or core-Hamiltonian like).
    symmetry_adapted_transform:
        AO-space transformation matrix to symmetry-adapted AO basis.
    basis_labels:
        Per-AO symmetry labels in the transformed basis.
    tolerance:
        Off-block absolute-coupling tolerance.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    m = _to_float_array(matrix)
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError("matrix must be a square 2D array")

    label_list = [str(lbl) for lbl in _as_list(basis_labels)]
    if len(label_list) != m.shape[0]:
        raise ValueError("basis_labels length must match matrix dimension")

    u = _to_float_array(symmetry_adapted_transform)
    if u.shape != (m.shape[0], m.shape[0]):
        raise ValueError("symmetry_adapted_transform must be square with matching dimension")

    # Transform into symmetry-adapted AO basis.
    transformed = u.T @ m @ u

    # Identify off-block entries by label mismatch.
    labels = np.array(label_list)
    mismatch = labels[:, None] != labels[None, :]
    off_block_values = np.abs(transformed[mismatch])
    if off_block_values.size == 0:
        max_off_block = 0.0
        max_indices: list[tuple[int, int]] = []
        off_count = 0
    else:
        max_off_block = float(np.max(off_block_values))
        flat_index = int(np.argmax(off_block_values))
        # Recover one flattened index from compressed mismatch view to AO pair indices.
        mismatch_positions = np.argwhere(mismatch)
        pair = mismatch_positions[flat_index]
        max_indices = [(int(pair[0]), int(pair[1]))]
        off_count = int(mismatch.sum())

    orthogonality = u.T @ u
    orthogonality_max_deviation = float(np.max(np.abs(orthogonality - np.eye(u.shape[0]))))
    orthogonality_ok = orthogonality_max_deviation <= tolerance

    payload = OneElectronBlockLeakSummary(
        max_off_block_abs=max_off_block,
        max_off_block_indices=max_indices,
        off_block_element_count=off_count,
        within_tolerance=(max_off_block <= tolerance),
        orthogonality_ok=orthogonality_ok,
        orthogonality_max_deviation=orthogonality_max_deviation,
        status="diagnostic_only_no_reductions",
    ).as_dict()
    payload["transform_shape"] = tuple(int(x) for x in u.shape)

    return payload


def update_one_electron_block_diagnostics(
    symmetry_metadata: Mapping[str, Any] | None,
    matrices: Mapping[str, Any],
    symmetry_adapted_transform: Any,
    basis_labels: Iterable[str],
    tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Attach named one-electron diagnostics under
    ``symmetry_metadata['one_electron_block_diagnostics']``.

    This is metadata-only and intentionally does not enable any symmetry
    acceleration behavior.
    """

    if symmetry_metadata is None:
        metadata: dict[str, Any] = {}
    else:
        metadata = dict(symmetry_metadata)

    diagnostics: dict[str, Any] = {}
    for key, matrix in matrices.items():
        diagnostics[str(key)] = one_electron_block_diagnostics(
            matrix, symmetry_adapted_transform, basis_labels,
            tolerance=tolerance,
        )

    metadata["one_electron_block_diagnostics"] = diagnostics
    metadata["use_integral_symmetry"] = False
    metadata["use_response_symmetry"] = False
    metadata["one_electron_block_diagnostics_status"] = "diagnostic_only_no_reductions"

    return metadata


def _cartesian_shell_size(l: int) -> int:
    """Return Cartesian basis-function count for a shell with angular momentum l."""

    if l < 0:
        raise ValueError("angular momentum must be non-negative")
    if l == 0:
        return 1
    if l == 1:
        return 3
    if l == 2:
        return 6
    raise ValueError("only s/p/d Cartesian shells are supported in this metadata-only scaffold")


def _component_signs(l: int, signs: tuple[int, int, int]) -> list[int]:
    """Per-component characters of a Cartesian shell under diag(sx, sy, sz).

    Cartesian d components are ordered xx, yy, zz, xy, xz, yz.
    """

    sx, sy, sz = signs
    if l == 0:
        return [1]
    if l == 1:
        return [sx, sy, sz]
    if l == 2:
        return [1, 1, 1, sx * sy, sx * sz, sy * sz]
    raise ValueError("only s/p/d Cartesian shells are supported in this metadata-only scaffold")


def _normalize_shells(shells: Iterable[Any]) -> list[tuple[int, int]]:
    """Normalize shell specs to (atom_index, l) tuples."""

    normalized: list[tuple[int, int]] = []
    for shell in shells:
        if isinstance(shell, Mapping):
            atom = int(shell["atom"])
            l = int(shell["l"])
        else:
            atom, l = (int(x) for x in shell)
        if atom < 0:
            raise ValueError("shell atom index must be non-negative")
        normalized.append((atom, l))
    return normalized


def _ao_operator_maps(
    shells: list[tuple[int, int]],
    operations: Iterable[Mapping[str, Any]],
) -> tuple[int, list[tuple[np.ndarray, np.ndarray]]]:
    """Build per-operation AO maps (target_index, sign) for sign-matrix ops.

    Operation payloads follow ``symmetry_detect``: ``matrix`` must be a
    diagonal +-1 matrix and ``permutation[a]`` is the atom that atom ``a``
    is mapped onto.
    """

    # AO offsets per shell, and per-atom shell lists for cross-atom mapping.
    offsets: list[int] = []
    n_ao = 0
    atom_shells: dict[int, list[int]] = {}
    for idx, (atom, l) in enumerate(shells):
        offsets.append(n_ao)
        n_ao += _cartesian_shell_size(l)
        atom_shells.setdefault(atom, []).append(idx)

    signatures = {
        atom: tuple(shells[idx][1] for idx in idx_list)
        for atom, idx_list in atom_shells.items()
    }

    maps: list[tuple[np.ndarray, np.ndarray]] = []
    for op in operations:
        matrix = _to_float_array(op["matrix"])
        diag = np.diagonal(matrix)
        if not np.allclose(matrix, np.diag(diag)) or not np.allclose(np.abs(diag), 1.0):
            raise ValueError("operations must be diagonal sign matrices (D2h family)")
        signs = tuple(int(round(s)) for s in diag)

        permutation = [int(a) for a in op["permutation"]]
        target = np.zeros(n_ao, dtype=int)
        sign = np.zeros(n_ao, dtype=float)
        for idx, (atom, l) in enumerate(shells):
            mapped_atom = permutation[atom]
            if signatures[atom] != signatures.get(mapped_atom):
                raise ValueError("symmetry-equivalent atoms must carry identical shell lists")
            shell_rank = atom_shells[atom].index(idx)
            mapped_shell = atom_shells[mapped_atom][shell_rank]
            comp_signs = _component_signs(l, signs)
            for comp, comp_sign in enumerate(comp_signs):
                target[offsets[idx] + comp] = offsets[mapped_shell] + comp
                sign[offsets[idx] + comp] = float(comp_sign)
        maps.append((target, sign))

    return n_ao, maps


def build_symmetry_adapted_transform(
    shells: Iterable[Any],
    operations: Iterable[Mapping[str, Any]],
    character_table: Mapping[str, Iterable[int]],
) -> tuple[np.ndarray, list[str]]:
    """Build an orthogonal symmetry-adapted AO transform with irrep labels.

    Parameters
    ----------
    shells:
        Shell specs as (atom_index, l) pairs or {'atom','l'} mappings,
        Cartesian s/p/d only, in AO order.
    operations:
        Abelian-group operations from ``symmetry_detect.detect_point_group``
        (``operations`` payload), in character-table column order.
    character_table:
        Irrep -> characters mapping matching the operation order.

    Returns
    -------
    (U, labels):
        Orthogonal (n_ao, n_ao) transform whose columns are SALCs, and the
        per-column irrep labels.
    """

    shell_list = _normalize_shells(shells)
    op_list = list(operations)
    n_ao, maps = _ao_operator_maps(shell_list, op_list)

    table = {str(irrep): [int(c) for c in row] for irrep, row in character_table.items()}
    for irrep, row in table.items():
        if len(row) != len(op_list):
            raise ValueError(f"character row for {irrep} does not match operation count")

    columns: list[np.ndarray] = []
    labels: list[str] = []
    order = len(op_list)
    for seed in range(n_ao):
        for irrep, chars in table.items():
            vec = np.zeros(n_ao)
            for (target, sign), chi in zip(maps, chars):
                vec[target[seed]] += chi * sign[seed]
            vec /= float(order)
            # Project out previously accepted SALCs.
            for col in columns:
                vec -= np.dot(col, vec) * col
            norm = float(np.linalg.norm(vec))
            if norm > 1.0e-8:
                columns.append(vec / norm)
                labels.append(irrep)

    if len(columns) != n_ao:
        raise ValueError(
            f"projection produced {len(columns)} SALCs for {n_ao} AOs; "
            "check operations/permutations consistency"
        )

    u = np.array(columns).T
    return u, labels


# Relative normalization of Cartesian d components (xx, yy, zz, xy, xz, yz):
# the cross monomials carry sqrt(3) relative to the squares.
_D_COMPONENT_NORMS = np.array([1.0, 1.0, 1.0, np.sqrt(3.0), np.sqrt(3.0), np.sqrt(3.0)])
_D_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def _shell_block(l: int, op_matrix: np.ndarray) -> np.ndarray:
    """Representation of an orthogonal 3x3 operation on one Cartesian shell.

    Acts on coefficient columns: c' = block @ c. For sign-diagonal operations
    this reduces to the diagonal of component signs.
    """

    if l == 0:
        return np.ones((1, 1))
    if l == 1:
        return op_matrix.copy()
    if l == 2:
        mono = np.zeros((6, 6))
        for col, (a, b) in enumerate(_D_PAIRS):
            for row, (cc, dd) in enumerate(_D_PAIRS):
                value = op_matrix[cc, a] * op_matrix[dd, b]
                if cc != dd:
                    value += op_matrix[dd, a] * op_matrix[cc, b]
                mono[row, col] = value
        # Rescale to normalized Cartesian components.
        return mono * (_D_COMPONENT_NORMS[None, :] / _D_COMPONENT_NORMS[:, None])
    raise ValueError("only s/p/d Cartesian shells are supported in this metadata-only scaffold")


def _ao_operator_matrix(
    shells: list[tuple[int, int]],
    op: Mapping[str, Any],
    matrix_key: str = "matrix",
) -> np.ndarray:
    """Dense AO representation of one (possibly input-frame) operation."""

    matrix = _to_float_array(op[matrix_key])
    if matrix.shape != (3, 3) or not np.allclose(matrix @ matrix.T, np.eye(3), atol=1.0e-10):
        raise ValueError("operation matrices must be orthogonal 3x3 arrays")

    offsets: list[int] = []
    n_ao = 0
    atom_shells: dict[int, list[int]] = {}
    for idx, (atom, l) in enumerate(shells):
        offsets.append(n_ao)
        n_ao += _cartesian_shell_size(l)
        atom_shells.setdefault(atom, []).append(idx)

    signatures = {
        atom: tuple(shells[idx][1] for idx in idx_list)
        for atom, idx_list in atom_shells.items()
    }

    permutation = [int(a) for a in op["permutation"]]
    t = np.zeros((n_ao, n_ao))
    for idx, (atom, l) in enumerate(shells):
        mapped_atom = permutation[atom]
        if signatures[atom] != signatures.get(mapped_atom):
            raise ValueError("symmetry-equivalent atoms must carry identical shell lists")
        shell_rank = atom_shells[atom].index(idx)
        mapped_shell = atom_shells[mapped_atom][shell_rank]
        block = _shell_block(l, matrix)
        size = _cartesian_shell_size(l)
        t[offsets[mapped_shell]:offsets[mapped_shell] + size,
          offsets[idx]:offsets[idx] + size] = block
    return t


def assign_mo_irreps(
    mo_coefficients: Any,
    overlap: Any,
    shells: Iterable[Any],
    operations: Iterable[Mapping[str, Any]],
    character_table: Mapping[str, Iterable[int]],
    tolerance: float = 1.0e-4,
    matrix_key: str = "matrix",
) -> dict[str, Any]:
    """Assign abelian irrep labels to molecular orbitals (metadata only).

    For each MO ``m`` the character under operation ``O`` is
    ``chi(O) = <m|S O|m> / <m|S|m>``; the MO gets the irrep whose character
    row matches within ``tolerance``, otherwise the label 'mixed'.

    ``matrix_key`` selects which operation matrix to use; pass
    ``'matrix_input_frame'`` when the MO coefficients live in the original
    input coordinates rather than the standard orientation.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    c = _to_float_array(mo_coefficients)
    if c.ndim != 2:
        raise ValueError("mo_coefficients must be a 2D (n_ao, n_mo) array")
    s = _to_float_array(overlap)
    if s.shape != (c.shape[0], c.shape[0]):
        raise ValueError("overlap must be square and match the AO dimension")

    shell_list = _normalize_shells(shells)
    op_list = list(operations)
    n_ao = sum(_cartesian_shell_size(l) for _, l in shell_list)
    if n_ao != c.shape[0]:
        raise ValueError("shells do not match the AO dimension of mo_coefficients")

    table = {str(irrep): [int(ch) for ch in row] for irrep, row in character_table.items()}

    sc = s @ c
    denominator = np.einsum('im,im->m', sc, c)
    characters = np.zeros((c.shape[1], len(op_list)))
    for iop, op in enumerate(op_list):
        t = _ao_operator_matrix(shell_list, op, matrix_key=matrix_key)
        transformed = t @ c
        numerator = np.einsum('im,im->m', sc, transformed)
        characters[:, iop] = numerator / denominator

    labels: list[str] = []
    deviations: list[float] = []
    for imo in range(c.shape[1]):
        best_label = 'mixed'
        best_dev = np.inf
        for irrep, row in table.items():
            dev = float(np.max(np.abs(characters[imo] - np.asarray(row, dtype=float))))
            if dev < best_dev:
                best_dev = dev
                best_label = irrep
        labels.append(best_label if best_dev <= tolerance else 'mixed')
        deviations.append(best_dev)

    return {
        'labels': labels,
        'characters': characters.tolist(),
        'max_deviation': float(np.max(deviations)) if deviations else 0.0,
        'status': 'label_only_no_reductions',
    }


def product_irrep(labels: Iterable[str], character_table: Mapping[str, Iterable[int]]) -> str:
    """Direct product of abelian irreps, e.g. b1 x b2 -> a2 in C2v.

    Returns 'mixed' if any input label is not in the table (e.g. a
    symmetry-broken 'mixed' orbital).
    """

    table = {str(irrep): np.asarray(list(row), dtype=int) for irrep, row in character_table.items()}
    if not table:
        raise ValueError("character table must not be empty")
    n_ops = len(next(iter(table.values())))

    product = np.ones(n_ops, dtype=int)
    for label in labels:
        row = table.get(str(label))
        if row is None:
            return 'mixed'
        product = product * row

    for irrep, row in table.items():
        if np.array_equal(row, product):
            return irrep
    return 'mixed'


def assign_state_irreps(
    amplitudes: Any,
    occ_coefficients: Any,
    vir_coefficients: Any,
    overlap: Any,
    shells: Iterable[Any],
    operations: Iterable[Mapping[str, Any]],
    character_table: Mapping[str, Iterable[int]],
    reference_labels: Iterable[str] = (),
    tolerance: float = 1.0e-3,
    matrix_key: str = "matrix",
) -> dict[str, Any]:
    """Assign abelian irrep labels to excitation amplitudes (metadata only).

    ``amplitudes`` holds X_ia per state, shape (n_states, n_occ, n_vir),
    where i indexes the occupied set described by ``occ_coefficients``
    (n_ao, n_occ) and a the virtual set of ``vir_coefficients`` (n_ao,
    n_vir). The transition character per operation is
    ``<X, U X V^T> / <X, X>`` with ``U/V`` the MO representations of the
    operation in the two sets. The total state irrep is the direct product
    of the transition irrep with ``reference_labels`` (e.g. the SOMOs of a
    spin-flip reference); for a closed-shell reference leave it empty.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    x = _to_float_array(amplitudes)
    if x.ndim != 3:
        raise ValueError("amplitudes must have shape (n_states, n_occ, n_vir)")
    c_occ = _to_float_array(occ_coefficients)
    c_vir = _to_float_array(vir_coefficients)
    s = _to_float_array(overlap)
    if c_occ.shape != (s.shape[0], x.shape[1]) or c_vir.shape != (s.shape[0], x.shape[2]):
        raise ValueError("coefficient shapes do not match amplitudes/overlap")

    shell_list = _normalize_shells(shells)
    op_list = list(operations)
    table = {str(irrep): [int(ch) for ch in row] for irrep, row in character_table.items()}

    norms = np.einsum('sia,sia->s', x, x)
    characters = np.zeros((x.shape[0], len(op_list)))
    for iop, op in enumerate(op_list):
        t = _ao_operator_matrix(shell_list, op, matrix_key=matrix_key)
        u = c_occ.T @ s @ t @ c_occ
        v = c_vir.T @ s @ t @ c_vir
        transformed = np.einsum('ji,sia,ba->sjb', u, x, v)
        characters[:, iop] = np.einsum('sia,sia->s', x, transformed) / norms

    reference = list(reference_labels)
    labels: list[str] = []
    transition_labels: list[str] = []
    deviations: list[float] = []
    for istate in range(x.shape[0]):
        best_label = 'mixed'
        best_dev = np.inf
        for irrep, row in table.items():
            dev = float(np.max(np.abs(characters[istate] - np.asarray(row, dtype=float))))
            if dev < best_dev:
                best_dev = dev
                best_label = irrep
        transition = best_label if best_dev <= tolerance else 'mixed'
        transition_labels.append(transition)
        labels.append(product_irrep([transition] + reference, table))
        deviations.append(best_dev)

    return {
        'labels': labels,
        'transition_labels': transition_labels,
        'reference_labels': reference,
        'characters': characters.tolist(),
        'max_deviation': float(np.max(deviations)) if deviations else 0.0,
        'status': 'label_only_no_reductions',
    }


def assign_mode_irreps(
    modes: Any,
    operations: Iterable[Mapping[str, Any]],
    character_table: Mapping[str, Iterable[int]],
    tolerance: float = 1.0e-3,
    matrix_key: str = "matrix",
) -> dict[str, Any]:
    """Assign abelian irrep labels to Cartesian normal modes (metadata only).

    ``modes`` holds one mode per row with 3*natom displacement components.
    Each operation acts by permuting atoms and rotating the per-atom
    displacement vectors; the character is ``<v|T v> / <v|v>``. Use
    ``matrix_key='matrix_input_frame'`` for modes in input coordinates.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    v = _to_float_array(modes)
    if v.ndim != 2 or v.shape[1] % 3 != 0:
        raise ValueError("modes must be a 2D (n_modes, 3*natom) array")
    natom = v.shape[1] // 3

    op_list = list(operations)
    table = {str(irrep): [int(ch) for ch in row] for irrep, row in character_table.items()}

    displacements = v.reshape(v.shape[0], natom, 3)
    norms = np.einsum('mai,mai->m', displacements, displacements)
    characters = np.zeros((v.shape[0], len(op_list)))
    for iop, op in enumerate(op_list):
        matrix = _to_float_array(op[matrix_key])
        if matrix.shape != (3, 3):
            raise ValueError("operation matrices must be 3x3 arrays")
        permutation = [int(a) for a in op["permutation"]]
        if len(permutation) != natom:
            raise ValueError("operation permutation does not match the atom count")
        transformed = np.zeros_like(displacements)
        for atom in range(natom):
            transformed[:, permutation[atom], :] = displacements[:, atom, :] @ matrix.T
        characters[:, iop] = np.einsum('mai,mai->m', displacements, transformed) / norms

    labels: list[str] = []
    deviations: list[float] = []
    for imode in range(v.shape[0]):
        best_label = 'mixed'
        best_dev = np.inf
        for irrep, row in table.items():
            dev = float(np.max(np.abs(characters[imode] - np.asarray(row, dtype=float))))
            if dev < best_dev:
                best_dev = dev
                best_label = irrep
        labels.append(best_label if best_dev <= tolerance else 'mixed')
        deviations.append(best_dev)

    return {
        'labels': labels,
        'characters': characters.tolist(),
        'max_deviation': float(np.max(deviations)) if deviations else 0.0,
        'status': 'label_only_no_reductions',
    }
