"""Small dense Full Configuration Interaction energy driver."""

from __future__ import annotations

import os
from dataclasses import dataclass
from itertools import combinations
import math
from math import comb

import numpy as np

import oqp
from oqp.library.cas_orbitals import active_core_lists
from oqp.utils.file_utils import dump_log, print_module_banner


FCI_INTEGRAL_BACKENDS = {"native"}
FCI_SOLVERS = {"auto", "dense", "davidson"}
CAS_ORBITAL_SOURCES = {"rhf", "json"}
CAS_LOCALIZE_OPTIONS = {"none"}
CAS_SORT_OPTIONS = {"energy", "none"}
CI_ROOT_TRACKING_OPTIONS = {"energy"}
STATE_AVERAGE_SPIN_BLOCKS = {"diagnostic"}
STATE_AVERAGE_ROOT_TRACKING = {"overlap"}
TARGET_SPIN_MULTIPLICITIES = {
    "singlet": 1,
    "doublet": 2,
    "triplet": 3,
    "quartet": 4,
    "quintet": 5,
    "sextet": 6,
    "septet": 7,
    "octet": 8,
    "nonet": 9,
}


def _bit_count(value: int) -> int:
    return bin(int(value)).count("1")


@dataclass(frozen=True)
class FCISettings:
    nroot: int = 1
    active_electrons: int | tuple[int, int] | list[int] = 0
    active_orbitals: int = 0
    frozen_core: int = 0
    max_det: int = 5000
    max_memory: int = 2048
    eig_tol: float = 1.0e-10
    integral_backend: str = "native"
    integral_cutoff: float = 5.0e-11
    solver: str = "auto"
    davidson_maxiter: int = 100
    davidson_subspace: int = 0
    print_ci_vectors: bool = False
    ci_print_threshold: float = 5.0e-2
    save_ci_vectors: bool = False
    save_rdm: bool = False
    state_average_enabled: bool = False
    state_average_weights: tuple[float, ...] = ()
    state_average_nstate: int = 0
    state_average_target_roots: tuple[int, ...] = ()
    state_average_equal_weights: bool = True
    active_orbital_indices: tuple[int, ...] = ()
    core_orbital_indices: tuple[int, ...] = ()
    orbital_source: str = "rhf"
    orbital_file: str = ""
    target_spin: str = "any"


def _annihilate(det: int, orb: int) -> tuple[int, int] | None:
    bit = 1 << int(orb)
    det = int(det)
    if not det & bit:
        return None
    phase = -1 if _bit_count(det & (bit - 1)) % 2 else 1
    return det ^ bit, phase


def _create(det: int, orb: int) -> tuple[int, int] | None:
    bit = 1 << int(orb)
    det = int(det)
    if det & bit:
        return None
    phase = -1 if _bit_count(det & (bit - 1)) % 2 else 1
    return det | bit, phase


def _occupied(det: int, norb: int) -> list[int]:
    return [orb for orb in range(norb) if det & (1 << orb)]


def _string_basis(norb: int, nelec: int) -> list[int]:
    return [sum(1 << orb for orb in occ) for occ in combinations(range(norb), nelec)]


def _determinants(norb: int, nelec: tuple[int, int]) -> list[int]:
    alpha = _string_basis(norb, nelec[0])
    beta = _string_basis(norb, nelec[1])
    return [a | (b << norb) for a in alpha for b in beta]


def _ci_vector_log_entries(
    ci_vectors: np.ndarray,
    determinants: list[int],
    norb: int,
    *,
    threshold: float,
    active_orbital_indices: str = "",
    root_indices: np.ndarray | list[int] | tuple[int, ...] | None = None,
) -> dict[str, object]:
    """Return determinant-resolved CI-vector rows for log output."""
    coeffs = _real_array(ci_vectors, "CI vectors")
    if coeffs.ndim == 1:
        coeffs = coeffs[:, None]
    if coeffs.ndim != 2 or coeffs.shape[0] != len(determinants):
        raise ValueError("CI vectors must have shape (ndet, nroot)")
    threshold = _direct_float_setting(threshold, "CI print threshold")
    if threshold < 0:
        raise ValueError("CI print threshold cannot be negative")

    labels = [token.strip() for token in str(active_orbital_indices or "").split(",") if token.strip()]
    if len(labels) != int(norb):
        labels = [str(index + 1) for index in range(int(norb))]
    if root_indices is None:
        roots = tuple(range(coeffs.shape[1]))
    else:
        root_values = _integer_vector(root_indices, "CI root indices")
        if root_values.shape != (coeffs.shape[1],):
            raise ValueError("CI root indices length must match CI vector roots")
        roots = tuple(int(index) for index in root_values)

    entries = []
    for det_index, det in enumerate(determinants):
        row = coeffs[det_index]
        if not np.any(np.abs(row) >= threshold):
            continue
        alpha_occ = _occupied(det, norb)
        beta_occ = _occupied(det >> norb, norb)
        entries.append(
            {
                "index": int(det_index),
                "alpha": ",".join(labels[orb] for orb in alpha_occ) or "-",
                "beta": ",".join(labels[orb] for orb in beta_occ) or "-",
                "coefficients": [float(value) for value in row],
            }
        )
    return {"threshold": threshold, "root_indices": roots, "entries": entries}


def _bool_setting(value, label: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "t", "yes", "y", "on", ".true."}:
            return True
        if normalized in {"0", "false", "f", "no", "n", "off", ".false.", ""}:
            return False
        raise ValueError(f"{label} must be a boolean")
    if isinstance(value, (int, np.integer)) and int(value) in (0, 1):
        return bool(value)
    raise ValueError(f"{label} must be a boolean")


def _int_setting(value, label: str, *, minimum: int | None = None) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be an integer")
    if isinstance(value, (int, np.integer)):
        parsed = int(value)
    elif isinstance(value, (float, np.floating)):
        number = float(value)
        if not math.isfinite(number) or not number.is_integer():
            raise ValueError(f"{label} must be an integer")
        parsed = int(number)
    elif isinstance(value, str):
        raw = value.strip()
        if not raw:
            raise ValueError(f"{label} must be an integer")
        try:
            parsed = int(raw, 10)
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer") from exc
    else:
        raise ValueError(f"{label} must be an integer")
    if minimum is not None and parsed < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    return parsed


def _float_setting(value, label: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be a finite real number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite real number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _direct_float_setting(value, label: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value,
        (int, float, np.integer, np.floating),
    ):
        raise ValueError(f"{label} must be a finite real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _real_array(values, label: str) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must contain finite real values") from exc
    if raw.dtype.kind not in {"i", "u", "f"}:
        raise ValueError(f"{label} must contain finite real values")
    array = np.ascontiguousarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must be finite")
    return array


def _integer_vector(values, label: str) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must contain integer root indices") from exc
    if raw.ndim == 0:
        raw = raw.reshape(1)
    if raw.size == 0:
        return np.ascontiguousarray([], dtype=np.int64)
    if raw.dtype.kind not in {"i", "u"}:
        raise ValueError(f"{label} must contain integer root indices")
    return np.ascontiguousarray(raw.reshape(-1), dtype=np.int64)


def _tuple_setting(value, cast, label: str) -> tuple:
    if value is None:
        return ()
    if isinstance(value, str):
        if not value.strip():
            return ()
        values = [item.strip() for item in value.split(",") if item.strip()]
    elif isinstance(value, (list, tuple)):
        values = value
    else:
        values = (value,)
    return tuple(cast(item, label) for item in values)


def _active_electrons_setting(value, label: str) -> int | tuple[int, int]:
    def parse_count(item, item_label):
        return _int_setting(item, item_label, minimum=0)

    if isinstance(value, str) and "," in value:
        values = _tuple_setting(value, parse_count, label)
    elif isinstance(value, (list, tuple)):
        values = _tuple_setting(value, parse_count, label)
    else:
        return _int_setting(value, label, minimum=0)
    if len(values) != 2:
        raise ValueError(f"{label} must be a non-negative integer or (nalpha, nbeta) pair")
    return (int(values[0]), int(values[1]))


def _choice_setting(value, label: str, choices: set[str]) -> str:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise ValueError(f"{label} must be one of {', '.join(sorted(choices))}")
    normalized = value.strip().lower()
    if normalized not in choices:
        raise ValueError(f"{label} must be one of {', '.join(sorted(choices))}")
    return normalized


def _string_setting(value, label: str) -> str:
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise ValueError(f"{label} must be a string")
    return value.strip()


def _target_spin_setting(value, label: str) -> str:
    if value is None:
        return "any"
    if isinstance(value, str):
        normalized = value.strip().lower()
        if not normalized:
            return "any"
        if normalized == "any" or normalized in TARGET_SPIN_MULTIPLICITIES:
            return normalized
        return str(_int_setting(normalized, label, minimum=1))
    return str(_int_setting(value, label, minimum=1))


def _save_npz_artifact(mol, suffix: str, **arrays) -> str:
    log_dir = str(getattr(mol, "log_path", "") or ".")
    project = os.path.basename(str(getattr(mol, "project_name", "") or "oqp"))
    if not project:
        project = "oqp"
    os.makedirs(log_dir, exist_ok=True)
    path = os.path.join(log_dir, f"{project}.{suffix}.npz")
    np.savez(path, **arrays)
    return path


def compute_s2(
    ci_vector: np.ndarray,
    determinants: list[int],
    norb: int,
    nelec: int | tuple[int, int] | list[int],
    *,
    det_index: dict[int, int] | None = None,
) -> float:
    """Return ``<S^2>`` for one CI vector in the determinant basis."""
    dets = list(determinants)
    coeff = _real_array(ci_vector, "CI vector").reshape(-1)
    if coeff.size != len(dets):
        raise ValueError("CI vector length must match the determinant count")

    norm = np.linalg.norm(coeff)
    if norm <= 0.0:
        raise ValueError("CI vector norm must be non-zero")
    coeff = coeff / norm

    nalpha, nbeta = _as_nelec_pair(nelec)
    ms = 0.5 * (nalpha - nbeta)
    s2 = ms * (ms + 1.0)

    if det_index is None:
        det_index = {det: idx for idx, det in enumerate(dets)}
    spin_flip = 0.0
    for col, det in enumerate(dets):
        c_col = coeff[col]
        if c_col == 0.0:
            continue
        for q in range(norb):
            ann_beta = _annihilate(det, norb + q)
            if ann_beta is None:
                continue
            det_q, phase_q = ann_beta
            cre_alpha_q = _create(det_q, q)
            if cre_alpha_q is None:
                continue
            det_alpha_q, phase_alpha_q = cre_alpha_q
            for p in range(norb):
                ann_alpha = _annihilate(det_alpha_q, p)
                if ann_alpha is None:
                    continue
                det_p, phase_p = ann_alpha
                cre_beta_p = _create(det_p, norb + p)
                if cre_beta_p is None:
                    continue
                det_beta_p, phase_beta_p = cre_beta_p
                row = det_index.get(det_beta_p)
                if row is not None:
                    spin_flip += (
                        coeff[row]
                        * phase_q
                        * phase_alpha_q
                        * phase_p
                        * phase_beta_p
                        * c_col
                    )

    return float(s2 + spin_flip)


def fci_spin_diagnostics(
    ci_vectors: np.ndarray,
    determinants: list[int],
    norb: int,
    nelec: int | tuple[int, int] | list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-root ``(<S^2>, multiplicity)`` diagnostics."""
    dets = list(determinants)
    coeffs = _real_array(ci_vectors, "CI vectors")
    if coeffs.ndim == 1:
        coeffs = coeffs[:, None]
    if coeffs.ndim != 2 or coeffs.shape[0] != len(dets):
        raise ValueError("CI vectors must have shape (ndet, nroot)")

    det_index = {det: idx for idx, det in enumerate(dets)}
    s2 = np.array(
        [
            compute_s2(coeffs[:, root], dets, norb, nelec, det_index=det_index)
            for root in range(coeffs.shape[1])
        ],
        dtype=np.float64,
    )
    multiplicity = np.rint(np.sqrt(np.maximum(0.0, 1.0 + 4.0 * s2))).astype(np.int64)
    multiplicity = np.maximum(multiplicity, 1)
    return np.ascontiguousarray(s2, dtype=np.float64), np.ascontiguousarray(
        multiplicity,
        dtype=np.int64,
    )


def _target_spin_multiplicity(target_spin: str) -> int | None:
    target = _target_spin_setting(target_spin, "target_spin")
    if target == "any":
        return None
    if target in TARGET_SPIN_MULTIPLICITIES:
        return TARGET_SPIN_MULTIPLICITIES[target]
    multiplicity = int(target)
    if multiplicity < 1:
        raise ValueError("target_spin multiplicity must be positive")
    return multiplicity


def _state_average_payload(
    settings: FCISettings,
    energies: np.ndarray,
    root_indices: np.ndarray | list[int] | tuple[int, ...] | None,
) -> dict[str, np.ndarray | float] | None:
    """Return normalized fixed-orbital state-average data for solved roots."""
    if not settings.state_average_enabled:
        return None

    energies = _state_average_real_vector(energies, "state-average energies")
    root_count = int(energies.shape[0])
    if root_count < 1:
        raise ValueError("State-average bookkeeping requires at least one solved CI root")

    configured_nstate = _int_setting(
        settings.state_average_nstate,
        "state_average.nstate",
        minimum=0,
    )
    requested_nroot = _int_setting(settings.nroot, "nroot", minimum=1)
    configured_target_roots = _state_average_integer_vector(
        settings.state_average_target_roots,
        "state_average.target_roots",
    )

    if configured_target_roots.size:
        target_roots = tuple(int(root) for root in configured_target_roots)
        if configured_nstate and configured_nstate != len(target_roots):
            raise ValueError("state_average.nstate must match state_average.target_roots length")
    else:
        nstate = int(configured_nstate or requested_nroot)
        if nstate < 1:
            raise ValueError("state_average.nstate must be positive when state averaging is enabled")
        target_roots = tuple(range(nstate))

    if not target_roots:
        raise ValueError("State-average target root list cannot be empty")
    if len(set(target_roots)) != len(target_roots):
        raise ValueError("state_average.target_roots must contain unique roots")
    if any(root < 0 or root >= root_count for root in target_roots):
        raise ValueError(
            f"State-average target roots must be between 0 and the solved root count minus one ({root_count - 1})"
        )

    if settings.state_average_equal_weights:
        weights = np.full(len(target_roots), 1.0 / len(target_roots), dtype=np.float64)
    else:
        weights = _state_average_real_vector(
            settings.state_average_weights,
            "state_average.weights",
        )
        if weights.shape != (len(target_roots),):
            raise ValueError(
                "state_average.weights length must match the number of state-average target roots"
            )
        if np.any(weights < 0.0):
            raise ValueError("state_average.weights cannot be negative")
        total = float(np.sum(weights))
        if total <= 0.0:
            raise ValueError("state_average.weights must contain a positive total weight")
        weights = weights / total

    roots = np.asarray(target_roots, dtype=np.int64)
    if root_indices is None:
        selected_root_indices = roots.copy()
    else:
        root_index_values = _state_average_integer_vector(
            root_indices,
            "state_average.root_indices",
        )
        if root_index_values.shape != (root_count,):
            raise ValueError("state_average.root_indices length must match solved CI roots")
        selected_root_indices = root_index_values[roots]
    state_average_energy = float(np.dot(weights, energies[roots]))

    return {
        "energy": state_average_energy,
        "weights": np.ascontiguousarray(weights, dtype=np.float64),
        "roots": np.ascontiguousarray(roots, dtype=np.int64),
        "root_indices": np.ascontiguousarray(selected_root_indices, dtype=np.int64),
    }


def _state_average_real_vector(values, label: str) -> np.ndarray:
    array = _real_array(values, label)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional array")
    return np.ascontiguousarray(array, dtype=np.float64)


def _state_average_integer_vector(values, label: str) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must contain integer root indices") from exc
    if raw.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional integer array")
    if raw.size == 0:
        return np.ascontiguousarray([], dtype=np.int64)
    if raw.dtype.kind not in {"i", "u"}:
        raise ValueError(f"{label} must contain integer root indices")
    return np.ascontiguousarray(raw, dtype=np.int64)


def _filter_roots_by_target_spin(
    energies: np.ndarray,
    coeffs: np.ndarray,
    s2: np.ndarray,
    multiplicity: np.ndarray,
    *,
    target_spin: str,
    requested_nroot: int,
    ci_label: str,
    ci_section: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return root arrays filtered to the requested spin multiplicity."""
    target_multiplicity = _target_spin_multiplicity(target_spin)
    root_indices = np.arange(np.asarray(energies).shape[0], dtype=np.int64)
    if target_multiplicity is None:
        return energies, coeffs, s2, multiplicity, root_indices

    keep = np.flatnonzero(np.asarray(multiplicity, dtype=np.int64) == target_multiplicity)
    if keep.size == 0:
        raise ValueError(
            f"{ci_label} target_spin={target_spin} found no matching roots among "
            f"{len(root_indices)} solved roots; choose another target_spin for this active space."
        )
    if keep.size < int(requested_nroot):
        raise ValueError(
            f"{ci_label} target_spin={target_spin} found only {keep.size} matching roots among "
            f"{len(root_indices)} solved roots; lower {ci_section} nroot or choose another target_spin."
        )
    keep = keep[: int(requested_nroot)]
    return (
        np.asarray(energies)[keep],
        np.asarray(coeffs)[:, keep],
        np.asarray(s2)[keep],
        np.asarray(multiplicity)[keep],
        root_indices[keep],
    )


def _electron_count_component(value, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must contain integer electron counts")
    count = int(value)
    if count < 0:
        raise ValueError(f"{label} cannot contain negative electron counts")
    return count


def _as_nelec_pair(nelec: int | tuple[int, int] | list[int]) -> tuple[int, int]:
    if isinstance(nelec, (tuple, list)):
        if len(nelec) != 2:
            raise ValueError("nelec must be an integer or an (nalpha, nbeta) pair")
        return (
            _electron_count_component(nelec[0], "nelec"),
            _electron_count_component(nelec[1], "nelec"),
        )
    total = _electron_count_component(nelec, "nelec")
    if total % 2:
        raise ValueError("integer nelec implies a closed-shell singlet and must be even")
    return total // 2, total // 2


def _spin_orbital_integrals(
    h1e: np.ndarray,
    eri: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    norb = h1e.shape[0]
    nspin = 2 * norb
    hspin = np.zeros((nspin, nspin), dtype=float)
    gspin = np.zeros((nspin, nspin, nspin, nspin), dtype=float)

    hspin[:norb, :norb] = h1e
    hspin[norb:, norb:] = h1e

    for p in range(nspin):
        ps = 0 if p < norb else 1
        pp = p % norb
        for q in range(nspin):
            qs = 0 if q < norb else 1
            qq = q % norb
            for r in range(nspin):
                if ps != (0 if r < norb else 1):
                    continue
                rr = r % norb
                for s in range(nspin):
                    if qs == (0 if s < norb else 1):
                        gspin[p, q, r, s] = eri[pp, rr, qq, s % norb]

    return hspin, gspin


def _iter_hamiltonian_elements(dets, det_index, hspin, gspin, nspin, cutoff):
    """Yield ``(row, col, value)`` contributions to the FCI Hamiltonian.

    This is the single second-quantized enumeration shared by the dense and the
    sparse/iterative paths: the ``a_p^+ a_q`` one-electron terms and the
    ``0.5 (pq|rs) a_p^+ a_q^+ a_s a_r`` two-electron terms in the spin-orbital
    basis. Several contributions to the same ``(row, col)`` are summed by the
    caller.
    """
    for col, det in enumerate(dets):
        occ = _occupied(det, nspin)

        for q in occ:
            ann_q = _annihilate(det, q)
            if ann_q is None:
                continue
            det_q, phase_q = ann_q
            for p in range(nspin):
                hval = hspin[p, q]
                if abs(hval) <= cutoff:
                    continue
                cre_p = _create(det_q, p)
                if cre_p is None:
                    continue
                det_pq, phase_p = cre_p
                row = det_index.get(det_pq)
                if row is not None:
                    yield row, col, hval * phase_q * phase_p

        for r in occ:
            ann_r = _annihilate(det, r)
            if ann_r is None:
                continue
            det_r, phase_r = ann_r
            for s in _occupied(det_r, nspin):
                ann_s = _annihilate(det_r, s)
                if ann_s is None:
                    continue
                det_rs, phase_s = ann_s
                for q in range(nspin):
                    cre_q = _create(det_rs, q)
                    if cre_q is None:
                        continue
                    det_qrs, phase_q = cre_q
                    for p in range(nspin):
                        gval = gspin[p, q, r, s]
                        if abs(gval) <= cutoff:
                            continue
                        cre_p = _create(det_qrs, p)
                        if cre_p is None:
                            continue
                        det_pqrs, phase_p = cre_p
                        row = det_index.get(det_pqrs)
                        if row is not None:
                            yield row, col, 0.5 * gval * phase_r * phase_s * phase_q * phase_p


def _fci_lib_threads() -> int:
    import os
    return max(1, min(8, os.cpu_count() or 1))


def _lib_dense_hamiltonian(dets, hspin, gspin, nspin, cutoff):
    """Dense determinant Hamiltonian through the liboqp OpenMP engine, or None
    when the symbol is missing (Python enumeration remains the fallback)."""
    backend = _lib_backend()
    if backend is None or nspin > 62:
        return None
    lib, ffi = backend
    if not hasattr(lib, "fci_dense_hamiltonian"):
        return None
    det_arr = np.ascontiguousarray(np.asarray(dets, dtype=np.int64))
    ndet = int(det_arr.size)
    h = np.zeros((ndet, ndet), dtype=np.float64)
    hs = np.ascontiguousarray(hspin, dtype=np.float64)
    gs = np.ascontiguousarray(gspin, dtype=np.float64)
    lib.fci_dense_hamiltonian(
        int(nspin), ndet, ffi.cast("int64_t *", det_arr.ctypes.data),
        ffi.cast("double *", hs.ctypes.data), ffi.cast("double *", gs.ctypes.data),
        float(cutoff), ffi.cast("double *", h.ctypes.data), _fci_lib_threads())
    # the engine symmetrizes in place; the C-order view of the symmetric
    # column-major result is the matrix itself
    return h


class _LibHamiltonianOperator:
    """Matrix-free determinant Hamiltonian (liboqp engine) with the sparse-
    matrix duck type _davidson consumes: .shape, .dot, .diagonal, .toarray.
    Applications are exactly 0.5*(H_raw + H_raw^T) x, matching the dense
    symmetrization."""

    def __init__(self, dets, hspin, gspin, nspin, cutoff):
        self._dets = np.ascontiguousarray(np.asarray(dets, dtype=np.int64))
        self._hs = np.ascontiguousarray(hspin, dtype=np.float64)
        self._gs = np.ascontiguousarray(gspin, dtype=np.float64)
        self._nspin = int(nspin)
        self._cutoff = float(cutoff)
        ndet = int(self._dets.size)
        self.shape = (ndet, ndet)

    def diagonal(self):
        lib, ffi = _lib_backend()
        diag = np.zeros(self.shape[0], dtype=np.float64)
        lib.fci_hamiltonian_diag(
            self._nspin, self.shape[0],
            ffi.cast("int64_t *", self._dets.ctypes.data),
            ffi.cast("double *", self._hs.ctypes.data),
            ffi.cast("double *", self._gs.ctypes.data),
            ffi.cast("double *", diag.ctypes.data))
        return diag

    def dot(self, x):
        lib, ffi = _lib_backend()
        x = np.asarray(x, dtype=np.float64)
        squeeze = x.ndim == 1
        xb = np.ascontiguousarray(x[:, None] if squeeze else x)
        y = np.zeros_like(xb)
        status = int(lib.fci_hamiltonian_matvec(
            self._nspin, self.shape[0],
            ffi.cast("int64_t *", self._dets.ctypes.data),
            ffi.cast("double *", self._hs.ctypes.data),
            ffi.cast("double *", self._gs.ctypes.data),
            self._cutoff, int(xb.shape[1]),
            ffi.cast("double *", xb.ctypes.data),
            ffi.cast("double *", y.ctypes.data), _fci_lib_threads()))
        if status != 0:
            raise RuntimeError("liboqp fci_hamiltonian_matvec failed")
        return y[:, 0] if squeeze else y

    def toarray(self):
        return _lib_dense_hamiltonian(self._dets, self._hs, self._gs,
                                      self._nspin, self._cutoff)


def _lib_matvec_available(nspin) -> bool:
    backend = _lib_backend()
    if backend is None or nspin > 62:
        return False
    lib, _ffi = backend
    return hasattr(lib, "fci_hamiltonian_matvec") and hasattr(lib, "fci_hamiltonian_diag")


def _build_dense_hamiltonian(dets, det_index, hspin, gspin, nspin, cutoff):
    """Assemble the explicit dense FCI Hamiltonian (8*ndet**2 bytes)."""
    lib_h = _lib_dense_hamiltonian(dets, hspin, gspin, nspin, cutoff)
    if lib_h is not None:
        return lib_h
    ndet = len(dets)
    hamiltonian = np.zeros((ndet, ndet), dtype=float)
    for row, col, value in _iter_hamiltonian_elements(
        dets, det_index, hspin, gspin, nspin, cutoff
    ):
        hamiltonian[row, col] += value
    return 0.5 * (hamiltonian + hamiltonian.T)


def _build_sparse_hamiltonian(dets, det_index, hspin, gspin, nspin, cutoff):
    """Assemble the FCI Hamiltonian as a sparse CSR matrix (no ndet**2 storage)."""
    import scipy.sparse as sp

    ndet = len(dets)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for row, col, value in _iter_hamiltonian_elements(
        dets, det_index, hspin, gspin, nspin, cutoff
    ):
        rows.append(row)
        cols.append(col)
        vals.append(value)
    hamiltonian = sp.coo_matrix(
        (vals, (rows, cols)), shape=(ndet, ndet), dtype=float
    ).tocsr()
    # The enumeration already produces a symmetric matrix; symmetrize defensively
    # (matching the dense path) so the iterative solver sees an exactly symmetric
    # operator.
    return (0.5 * (hamiltonian + hamiltonian.transpose())).tocsr()


def _davidson_subspace_limit(ndet: int, nroot: int, davidson_subspace: int = 0) -> int:
    requested = int(davidson_subspace)
    if requested < 0:
        raise ValueError("FCI Davidson subspace override cannot be negative")
    if requested == 0:
        requested = max(4 * int(nroot), 2 * int(nroot) + 20)
    if requested < int(nroot):
        raise ValueError("FCI Davidson subspace must be 0 or at least nroot")
    return min(int(ndet), requested)


def _symmetric_eigh_jacobi(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Diagonalize a small real symmetric matrix without external LAPACK.

    OpenQP is built against ILP64 BLAS/LAPACK, while common Python NumPy builds
    may be LP64. Calling ``np.linalg.eigh`` after loading the native backend can
    therefore hit an incompatible LAPACK ABI. The active-space dense CI path is
    intentionally small, so a Jacobi solver is a robust local diagonalizer here.
    """
    a = np.array(matrix, dtype=np.float64, copy=True)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("symmetric eigensolver requires a square matrix")
    n = int(a.shape[0])
    if n == 0:
        raise ValueError("symmetric eigensolver requires a non-empty matrix")
    if n == 1:
        return np.array([float(a[0, 0])], dtype=np.float64), np.eye(1, dtype=np.float64)

    # Start from the explicitly symmetrized matrix used by the old LAPACK path.
    for row in range(n):
        for col in range(row + 1, n):
            value = 0.5 * (float(a[row, col]) + float(a[col, row]))
            a[row, col] = value
            a[col, row] = value

    eigvecs = np.eye(n, dtype=np.float64)
    scale = max(float(np.max(np.abs(a))), 1.0)
    tol = max(1.0e-14 * scale, np.finfo(np.float64).eps * n * scale)
    max_rotations = max(100, 25 * n * n * n)

    for _ in range(max_rotations):
        pivot_value = 0.0
        pivot = (0, 1)
        for row in range(n - 1):
            for col in range(row + 1, n):
                value = abs(float(a[row, col]))
                if value > pivot_value:
                    pivot_value = value
                    pivot = (row, col)
        if pivot_value <= tol:
            break

        p, q = pivot
        apq = float(a[p, q])
        if apq == 0.0:
            continue
        app = float(a[p, p])
        aqq = float(a[q, q])
        tau = (aqq - app) / (2.0 * apq)
        if tau == 0.0:
            tangent = 1.0
        else:
            tangent = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine

        a[p, p] = app - tangent * apq
        a[q, q] = aqq + tangent * apq
        a[p, q] = 0.0
        a[q, p] = 0.0
        for row in range(n):
            if row == p or row == q:
                continue
            arp = float(a[row, p])
            arq = float(a[row, q])
            new_rp = cosine * arp - sine * arq
            new_rq = sine * arp + cosine * arq
            a[row, p] = new_rp
            a[p, row] = new_rp
            a[row, q] = new_rq
            a[q, row] = new_rq

        for row in range(n):
            vrp = float(eigvecs[row, p])
            vrq = float(eigvecs[row, q])
            eigvecs[row, p] = cosine * vrp - sine * vrq
            eigvecs[row, q] = sine * vrp + cosine * vrq
    else:
        raise ValueError(
            f"FCI dense Jacobi diagonalization did not converge after {max_rotations} rotations"
        )

    eigvals = np.diag(a).copy()
    order = np.argsort(eigvals, kind="mergesort")
    return (
        np.ascontiguousarray(eigvals[order], dtype=np.float64),
        np.ascontiguousarray(eigvecs[:, order], dtype=np.float64),
    )


def _eigen_residual_norm(
    matrix: np.ndarray,
    eigenvalue: float,
    eigenvector: np.ndarray,
) -> float:
    # One BLAS matvec.  The original pure-Python double loop made this O(n^2)
    # Python work PER ROOT -- verifying every root of a 4900-determinant CASCI
    # cost ~20 minutes while the eigh itself took seconds.
    vec = np.asarray(eigenvector, dtype=np.float64)
    residual = np.asarray(matrix, dtype=np.float64) @ vec - float(eigenvalue) * vec
    return float(np.linalg.norm(residual))


def _lib_backend():
    """(lib, ffi) of the loaded native core, or None when unavailable."""
    try:
        from oqp import lib, ffi
    except Exception:
        return None
    return lib, ffi


def _lib_dsyevd(sym: np.ndarray):
    """Diagonalize through liboqp's own ILP64 LAPACK (eigen::diag_symm_full).

    ABI-proof by construction (same integer width as the linked LAPACK).
    Returns (eigvals, eigvecs) or None when the symbol is unavailable or
    LAPACK reports failure.  LAPACK writes eigenvectors column-major into the
    buffer, so the C-order view is transposed back."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "oqp_dsyevd"):
        return None
    n = int(sym.shape[0])
    a = np.ascontiguousarray(sym, dtype=np.float64).copy()
    w = np.zeros(n, dtype=np.float64)
    info = int(lib.oqp_dsyevd(n, ffi.cast("double *", a.ctypes.data),
                              ffi.cast("double *", w.ctypes.data)))
    if info != 0:
        return None
    return w, np.ascontiguousarray(a.T)


def _symmetric_eigh(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Diagonalize a real symmetric matrix with a LAPACK fallback guard."""
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("symmetric eigensolver requires a square matrix")
    sym = np.ascontiguousarray(0.5 * (matrix + matrix.T), dtype=np.float64)
    scale = max(float(np.max(np.abs(sym))), 1.0) if sym.size else 1.0
    residual_tol = max(1.0e-8 * scale, 100.0 * np.finfo(np.float64).eps * sym.shape[0] * scale)
    def _verified(eigvals, eigvecs):
        if not (np.all(np.isfinite(eigvals)) and np.all(np.isfinite(eigvecs))):
            return None
        # batched residual: max_k || A v_k - w_k v_k ||, one BLAS matmul
        worst = float(np.max(np.linalg.norm(
            sym @ eigvecs - eigvecs * eigvals[None, :], axis=0))) if sym.size else 0.0
        if worst > residual_tol:
            return None
        return (np.ascontiguousarray(eigvals, dtype=np.float64),
                np.ascontiguousarray(eigvecs, dtype=np.float64))

    try:
        result = _verified(*np.linalg.eigh(sym))
        if result is not None:
            return result
    except Exception:
        pass
    # numpy's LAPACK failed the residual gate (e.g. an LP64/ILP64 ABI clash):
    # retry through liboqp's own LAPACK before falling back to Jacobi.
    lib_result = _lib_dsyevd(sym)
    if lib_result is not None:
        result = _verified(*lib_result)
        if result is not None:
            return result
    return _symmetric_eigh_jacobi(sym)


def _davidson(hamiltonian, diag, nroot, *, tol=1.0e-10, max_iter=100, max_subspace=0):
    """Lowest ``nroot`` eigenpairs of a symmetric operator via block Davidson.

    ``hamiltonian`` only needs a ``.dot`` matrix-vector product; ``diag`` is the
    operator diagonal used as the Davidson preconditioner.
    """
    import scipy.sparse as sp

    ndet = hamiltonian.shape[0]
    max_subspace = _davidson_subspace_limit(ndet, nroot, max_subspace)
    # Tiny spaces: a direct dense diagonalization is cheaper and more robust than
    # iterating, and it also covers nroot close to ndet.
    if ndet <= max(2 * nroot + 20, 40):
        dense = (hamiltonian.toarray() if hasattr(hamiltonian, "toarray")
                 else np.asarray(hamiltonian))
        eigvals, eigvecs = _symmetric_eigh(0.5 * (dense + dense.T))
        return eigvals[:nroot], eigvecs[:, :nroot]

    diag = np.asarray(diag, dtype=float)
    matvec = hamiltonian.dot

    # Initial guess: unit vectors on the nroot lowest diagonal elements.
    order = np.argsort(diag)
    basis = np.zeros((ndet, nroot), dtype=float)
    for i in range(nroot):
        basis[order[i], i] = 1.0
    basis, _ = np.linalg.qr(basis)
    sigma = matvec(basis)

    residual_norms = np.full(nroot, np.inf)
    for _ in range(max_iter):
        sub = basis.T @ sigma
        sub = 0.5 * (sub + sub.T)
        theta, vecs = _symmetric_eigh(sub)
        theta = theta[:nroot]
        vecs = vecs[:, :nroot]
        ritz = basis @ vecs
        residual = sigma @ vecs - ritz * theta
        residual_norms = np.linalg.norm(residual, axis=0)
        if np.max(residual_norms) < tol:
            return theta, ritz

        # Preconditioned corrections for the unconverged roots.
        additions = []
        for i in range(nroot):
            if residual_norms[i] < tol:
                continue
            denom = theta[i] - diag
            denom = np.where(np.abs(denom) < 1.0e-12, 1.0e-12, denom)
            additions.append(residual[:, i] / denom)

        # Restart from the current Ritz vectors when the subspace would overflow.
        if not additions or basis.shape[1] + len(additions) > max_subspace:
            basis, _ = np.linalg.qr(ritz)
            sigma = matvec(basis)
            continue

        new_cols = []
        current = basis
        for vec in additions:
            for _ in range(2):  # orthogonalize twice for numerical stability
                vec = vec - current @ (current.T @ vec)
            norm = np.linalg.norm(vec)
            if norm > 1.0e-8:
                vec = vec / norm
                current = np.column_stack([current, vec])
                new_cols.append(vec)
        if not new_cols:
            return theta, ritz
        sigma = np.column_stack([sigma, matvec(np.column_stack(new_cols))])
        basis = current

    raise ValueError(
        f"FCI Davidson did not converge in {max_iter} iterations "
        f"(max residual {np.max(residual_norms):.3e} > eig_tol={tol:.3e})"
    )


def solve_fci(
    h1e: np.ndarray,
    eri: np.ndarray,
    nelec: int | tuple[int, int] | list[int],
    *,
    ecore: float = 0.0,
    nroot: int = 1,
    max_det: int = 5000,
    max_memory: int = 2048,
    eig_tol: float = 1.0e-10,
    integral_cutoff: float = 0.0,
    solver: str = "auto",
    davidson_maxiter: int = 100,
    davidson_subspace: int = 0,
    target_spin: str = "any",
    active_section: str = "[fci]",
    ci_section: str = "[fci]",
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the FCI eigenproblem in a determinant basis.

    The two-electron tensor uses chemists' notation ``(pq|rs)`` over spatial
    orbitals. Two solvers are available:

    * ``"dense"`` builds the explicit ``ndet x ndet`` Hamiltonian and calls
      :func:`numpy.linalg.eigh` (exact, but ``8*ndet**2`` bytes of storage);
    * ``"davidson"`` builds a sparse Hamiltonian and runs a block Davidson
      iteration for the lowest ``nroot`` states, avoiding the dense matrix.

    ``"auto"`` (the default) uses the dense solver while its Hamiltonian fits the
    ``max_memory`` budget and switches to Davidson once it would not. Both paths
    share the same verified second-quantized matrix-element enumeration, so they
    agree to numerical precision. When ``target_spin`` is not ``"any"``, the
    solver expands the internal root window as needed and returns the lowest
    roots matching the requested spin multiplicity. This is still meant for
    small active spaces.
    """

    h1e = _real_array(h1e, "h1e")
    eri = _real_array(eri, "eri")
    if h1e.ndim != 2 or h1e.shape[0] != h1e.shape[1]:
        raise ValueError("h1e must be a square matrix")
    norb = h1e.shape[0]
    if eri.shape != (norb, norb, norb, norb):
        raise ValueError("eri must have shape (norb, norb, norb, norb)")

    nroot = _int_setting(nroot, "nroot", minimum=1)
    max_det = _int_setting(max_det, "max_det", minimum=1)
    max_memory = _int_setting(max_memory, "max_memory", minimum=1)
    davidson_maxiter = _int_setting(
        davidson_maxiter,
        f"{ci_section} davidson_maxiter",
        minimum=1,
    )
    davidson_subspace = _int_setting(
        davidson_subspace,
        f"{ci_section} davidson_subspace",
        minimum=0,
    )
    ecore = _float_setting(ecore, "ecore")
    eig_tol = _float_setting(eig_tol, "eig_tol")
    integral_cutoff = _float_setting(integral_cutoff, "integral_cutoff")
    target_spin = _target_spin_setting(target_spin, "target_spin")
    target_multiplicity = _target_spin_multiplicity(target_spin)
    if eig_tol <= 0.0:
        raise ValueError("eig_tol must be positive")
    if integral_cutoff < 0.0:
        raise ValueError("integral_cutoff cannot be negative")

    nalpha, nbeta = _as_nelec_pair(nelec)
    if min(nalpha, nbeta) < 0 or max(nalpha, nbeta) > norb:
        raise ValueError("electron count is incompatible with the orbital count")

    ndet = comb(norb, nalpha) * comb(norb, nbeta)
    if ndet < 1:
        raise ValueError("FCI determinant space is empty")
    if ndet > max_det:
        raise ValueError(
            f"FCI determinant space has {ndet} determinants, exceeding max_det={max_det}. "
            f"Reduce the active space ({active_section} frozen_core/active_orbitals) "
            f"or raise {active_section} max_det."
        )
    if nroot < 1 or nroot > ndet:
        raise ValueError(f"nroot must be between 1 and the determinant count ({ndet})")

    solver = _choice_setting(solver, "solver", FCI_SOLVERS)
    if davidson_subspace and davidson_subspace < nroot:
        raise ValueError(f"{ci_section} davidson_subspace must be 0 or at least nroot ({nroot})")

    # The explicit dense Hamiltonian costs 8*ndet**2 bytes; the Davidson path only
    # needs a sparse Hamiltonian plus a handful of length-ndet vectors.
    dense_bytes = 8 * ndet * ndet
    budget_bytes = max_memory * 1024 * 1024
    if solver == "auto":
        solver = "dense" if dense_bytes <= budget_bytes else "davidson"
    if solver == "dense" and dense_bytes > budget_bytes:
        raise ValueError(
            f"FCI dense Hamiltonian for {ndet} determinants needs "
            f"~{dense_bytes / 1024 ** 3:.2f} GiB, exceeding the {active_section} max_memory budget "
            f"of {max_memory} MiB. Reduce the active space, raise {active_section} max_memory, "
            f"or use the iterative solver ({ci_section} solver=davidson)."
        )

    dets = _determinants(norb, (nalpha, nbeta))
    det_index = {det: idx for idx, det in enumerate(dets)}
    nspin = 2 * norb
    hspin, gspin = _spin_orbital_integrals(h1e, eri)

    if solver == "dense":
        hamiltonian = _build_dense_hamiltonian(
            dets, det_index, hspin, gspin, nspin, integral_cutoff
        )
        eigvals, eigvecs = _symmetric_eigh(hamiltonian)
        if target_multiplicity is None:
            selected_indices = np.arange(nroot, dtype=np.int64)
            selected_eigvals = eigvals[:nroot]
            selected_eigvecs = eigvecs[:, :nroot]
        else:
            # Spin-classify only a low-root PREFIX and grow it on demand: the
            # diagnostics cost is per-vector, and classifying all ndet
            # eigenvectors of a large dense solve (e.g. 4900 for CAS(8,8))
            # dominated the entire calculation while the filter only ever
            # consumes the lowest few roots.
            limit = min(ndet, max(4 * nroot, nroot + 12))
            while True:
                s2, multiplicity = fci_spin_diagnostics(
                    eigvecs[:, :limit], dets, norb, (nalpha, nbeta))
                try:
                    (selected_eigvals, selected_eigvecs, _s2, _multiplicity,
                     selected_indices) = _filter_roots_by_target_spin(
                        eigvals[:limit],
                        eigvecs[:, :limit],
                        s2,
                        multiplicity,
                        target_spin=target_spin,
                        requested_nroot=nroot,
                        ci_label="FCI",
                        ci_section=ci_section,
                    )
                    break
                except ValueError:
                    if limit >= ndet:
                        raise
                    limit = min(ndet, 2 * limit)

        for root in selected_indices:
            residual = _eigen_residual_norm(hamiltonian, eigvals[root], eigvecs[:, root])
            if residual > eig_tol:
                raise ValueError(
                    f"FCI diagonalization residual {residual:.3e} exceeds eig_tol={eig_tol:.3e}"
                )
    else:
        # Iterative path: never forms the dense ndet x ndet matrix.  The
        # liboqp matrix-free operator (Fortran/OpenMP enumeration per
        # application) is preferred; the Python-built scipy sparse matrix
        # remains the fallback and the numerical pin.
        if _lib_matvec_available(nspin):
            hamiltonian = _LibHamiltonianOperator(
                dets, hspin, gspin, nspin, integral_cutoff
            )
        else:
            hamiltonian = _build_sparse_hamiltonian(
                dets, det_index, hspin, gspin, nspin, integral_cutoff
            )
        solve_nroot = nroot
        while True:
            effective_subspace = (
                max(davidson_subspace, solve_nroot) if davidson_subspace else 0
            )
            max_subspace = _davidson_subspace_limit(
                ndet,
                solve_nroot,
                effective_subspace,
            )
            work_bytes = 8 * ndet * (2 * max_subspace + 4 * solve_nroot)
            if work_bytes > budget_bytes:
                raise ValueError(
                    f"FCI Davidson working set for {ndet} determinants needs "
                    f"~{work_bytes / 1024 ** 3:.2f} GiB, exceeding the {active_section} max_memory "
                    f"budget of {max_memory} MiB. Reduce the active space or raise {active_section} max_memory."
                )
            eigvals, eigvecs = _davidson(
                hamiltonian,
                hamiltonian.diagonal(),
                solve_nroot,
                tol=eig_tol,
                max_iter=davidson_maxiter,
                max_subspace=effective_subspace,
            )
            if target_multiplicity is None:
                selected_eigvals = eigvals[:nroot]
                selected_eigvecs = eigvecs[:, :nroot]
                break

            s2, multiplicity = fci_spin_diagnostics(eigvecs, dets, norb, (nalpha, nbeta))
            try:
                selected_eigvals, selected_eigvecs, _s2, _multiplicity, _root_indices = (
                    _filter_roots_by_target_spin(
                        eigvals,
                        eigvecs,
                        s2,
                        multiplicity,
                        target_spin=target_spin,
                        requested_nroot=nroot,
                        ci_label="FCI",
                        ci_section=ci_section,
                    )
                )
                break
            except ValueError:
                if solve_nroot >= ndet:
                    raise
                solve_nroot = min(ndet, max(solve_nroot + 1, 2 * solve_nroot))

    return selected_eigvals + ecore, selected_eigvecs


def _unpack_lower_triangle(packed: np.ndarray, n: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            matrix[i, j] = packed[idx]
            matrix[j, i] = packed[idx]
            idx += 1
    return matrix


def _transform_integrals(
    hcore_ao: np.ndarray,
    eri_ao: np.ndarray,
    coeff: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    h1e = coeff.T @ hcore_ao @ coeff
    eri = np.einsum(
        "up,vq,wr,xs,uvwx->pqrs",
        coeff,
        coeff,
        coeff,
        coeff,
        eri_ao,
        optimize=True,
    )
    return h1e, eri


def _validate_ci_settings(
    settings: FCISettings,
    *,
    ci_section: str,
) -> FCISettings:
    if settings.eig_tol <= 0.0:
        raise ValueError(f"{ci_section}.eig_tol must be positive")
    if settings.integral_cutoff < 0.0:
        raise ValueError(f"{ci_section}.integral_cutoff cannot be negative")
    if settings.davidson_subspace and settings.davidson_subspace < settings.nroot:
        raise ValueError(
            f"{ci_section}.davidson_subspace must be 0 or at least {ci_section}.nroot"
        )
    if settings.ci_print_threshold < 0.0:
        raise ValueError(f"{ci_section}.ci_print_threshold cannot be negative")
    return settings


def _validate_casci_settings(settings: FCISettings) -> FCISettings:
    active = tuple(settings.active_orbital_indices or ())
    core = tuple(settings.core_orbital_indices or ())

    if active and len(set(active)) != len(active):
        raise ValueError("cas.active_orbital_indices must be unique")
    if core and len(set(core)) != len(core):
        raise ValueError("cas.core_orbital_indices must be unique")
    if core and not active:
        raise ValueError(
            "cas.active_orbital_indices is required when cas.core_orbital_indices is set"
        )
    if active and core and set(active) & set(core):
        raise ValueError(
            "cas.core_orbital_indices must be disjoint from cas.active_orbital_indices"
        )
    if active and settings.active_orbitals and settings.active_orbitals != len(active):
        raise ValueError(
            "cas.active_orbitals must match the length of cas.active_orbital_indices"
        )

    if settings.state_average_enabled:
        target_roots = tuple(settings.state_average_target_roots or ())
        if len(set(target_roots)) != len(target_roots):
            raise ValueError("state_average.target_roots must be unique")
        if target_roots and max(target_roots) >= settings.nroot:
            raise ValueError("state_average.target_roots must refer to solved CI root slots")
        if target_roots and settings.state_average_nstate not in (0, len(target_roots)):
            raise ValueError(
                "state_average.nstate must match state_average.target_roots length"
            )

        averaged_root_count = (
            len(target_roots)
            if target_roots
            else int(settings.state_average_nstate or settings.nroot)
        )
        if averaged_root_count < 1:
            raise ValueError("state_average.nstate must request at least one root")
        if averaged_root_count > settings.nroot:
            raise ValueError("state_average.nstate cannot exceed ci.nroot")

        if not settings.state_average_equal_weights:
            weights = tuple(settings.state_average_weights or ())
            if len(weights) != averaged_root_count:
                raise ValueError(
                    "state_average.weights must match the averaged root count"
                )
            if any(weight < 0.0 for weight in weights):
                raise ValueError("state_average.weights cannot be negative")
            if sum(weights) <= 0.0:
                raise ValueError("state_average.weights must have a positive total weight")
    return settings


def _config_section(config: dict, section: str, label: str) -> dict:
    if not isinstance(config, dict):
        raise ValueError(f"{label} config must be a mapping")
    values = config.get(section, {})
    if values is None:
        return {}
    if not isinstance(values, dict):
        raise ValueError(f"{label} [{section}] config section must be a mapping")
    return values


def _settings_from_config(config: dict) -> FCISettings:
    raw = _config_section(config, "fci", "FCI")
    settings = FCISettings(
        nroot=_int_setting(raw.get("nroot", 1), "fci.nroot", minimum=1),
        active_electrons=_active_electrons_setting(
            raw.get("active_electrons", 0),
            "fci.active_electrons",
        ),
        active_orbitals=_int_setting(
            raw.get("active_orbitals", 0),
            "fci.active_orbitals",
            minimum=0,
        ),
        frozen_core=_int_setting(raw.get("frozen_core", 0), "fci.frozen_core", minimum=0),
        max_det=_int_setting(raw.get("max_det", 5000), "fci.max_det", minimum=1),
        max_memory=_int_setting(raw.get("max_memory", 2048), "fci.max_memory", minimum=1),
        eig_tol=_float_setting(raw.get("eig_tol", 1.0e-10), "fci.eig_tol"),
        integral_backend=_choice_setting(
            raw.get("integral_backend", "native"),
            "fci.integral_backend",
            FCI_INTEGRAL_BACKENDS,
        ),
        integral_cutoff=_float_setting(
            raw.get("integral_cutoff", 5.0e-11),
            "fci.integral_cutoff",
        ),
        solver=_choice_setting(raw.get("solver", "auto"), "fci.solver", FCI_SOLVERS),
        davidson_maxiter=_int_setting(
            raw.get("davidson_maxiter", 100),
            "fci.davidson_maxiter",
            minimum=1,
        ),
        davidson_subspace=_int_setting(
            raw.get("davidson_subspace", 0),
            "fci.davidson_subspace",
            minimum=0,
        ),
        print_ci_vectors=_bool_setting(
            raw.get("print_ci_vectors", False),
            "fci.print_ci_vectors",
        ),
        ci_print_threshold=_float_setting(
            raw.get("ci_print_threshold", 5.0e-2),
            "fci.ci_print_threshold",
        ),
        save_ci_vectors=_bool_setting(
            raw.get("save_ci_vectors", False),
            "fci.save_ci_vectors",
        ),
        save_rdm=_bool_setting(raw.get("save_rdm", False), "fci.save_rdm"),
        target_spin=_target_spin_setting(raw.get("target_spin", "any"), "fci.target_spin"),
    )
    return _validate_ci_settings(settings, ci_section="fci")


def settings_from_casci_config(config: dict) -> FCISettings:
    """Map first-class CAS/CI input sections onto the verified CI backend.

    CASCI intentionally reads only the first-class ``[cas]`` and ``[ci]``
    sections.  The legacy ``[fci]`` section remains reserved for
    ``method=fci`` so validation and runtime behavior cannot diverge.
    """
    cas = _config_section(config, "cas", "CASCI")
    ci = _config_section(config, "ci", "CASCI")
    state_average = _config_section(config, "state_average", "CASCI")
    _choice_setting(cas.get("localize", "none"), "cas.localize", CAS_LOCALIZE_OPTIONS)
    _choice_setting(
        cas.get("sort_orbitals", "energy"),
        "cas.sort_orbitals",
        CAS_SORT_OPTIONS,
    )
    _choice_setting(
        ci.get("root_tracking", "energy"),
        "ci.root_tracking",
        CI_ROOT_TRACKING_OPTIONS,
    )
    _choice_setting(
        state_average.get("spin_blocks", "diagnostic"),
        "state_average.spin_blocks",
        STATE_AVERAGE_SPIN_BLOCKS,
    )
    _choice_setting(
        state_average.get("root_tracking", "overlap"),
        "state_average.root_tracking",
        STATE_AVERAGE_ROOT_TRACKING,
    )
    spin_adapted = _bool_setting(ci.get("spin_adapted", False), "ci.spin_adapted")
    if spin_adapted:
        raise ValueError("ci.spin_adapted is not implemented; set ci.spin_adapted=false")
    settings = FCISettings(
        nroot=_int_setting(ci.get("nroot", 1), "ci.nroot", minimum=1),
        active_electrons=_active_electrons_setting(
            cas.get("active_electrons", 0),
            "cas.active_electrons",
        ),
        active_orbitals=_int_setting(
            cas.get("active_orbitals", 0),
            "cas.active_orbitals",
            minimum=0,
        ),
        frozen_core=_int_setting(cas.get("frozen_core", 0), "cas.frozen_core", minimum=0),
        max_det=_int_setting(cas.get("max_det", 5000), "cas.max_det", minimum=1),
        max_memory=_int_setting(cas.get("max_memory", 2048), "cas.max_memory", minimum=1),
        eig_tol=_float_setting(ci.get("eig_tol", 1.0e-10), "ci.eig_tol"),
        integral_backend=_choice_setting(
            ci.get("integral_backend", "native"),
            "ci.integral_backend",
            FCI_INTEGRAL_BACKENDS,
        ),
        integral_cutoff=_float_setting(ci.get("integral_cutoff", 5.0e-11), "ci.integral_cutoff"),
        solver=_choice_setting(ci.get("solver", "auto"), "ci.solver", FCI_SOLVERS),
        davidson_maxiter=_int_setting(
            ci.get("davidson_maxiter", 100),
            "ci.davidson_maxiter",
            minimum=1,
        ),
        davidson_subspace=_int_setting(
            ci.get("davidson_subspace", 0),
            "ci.davidson_subspace",
            minimum=0,
        ),
        print_ci_vectors=_bool_setting(
            ci.get("print_ci_vectors", False),
            "ci.print_ci_vectors",
        ),
        ci_print_threshold=_float_setting(
            ci.get("ci_print_threshold", 5.0e-2),
            "ci.ci_print_threshold",
        ),
        save_ci_vectors=_bool_setting(
            ci.get("save_ci_vectors", False),
            "ci.save_ci_vectors",
        ),
        save_rdm=_bool_setting(ci.get("save_rdm", False), "ci.save_rdm"),
        state_average_enabled=_bool_setting(
            state_average.get("enabled", False),
            "state_average.enabled",
        ),
        state_average_weights=_tuple_setting(
            state_average.get("weights", ()),
            _float_setting,
            "state_average.weights",
        ),
        state_average_nstate=_int_setting(
            state_average.get("nstate", 0),
            "state_average.nstate",
            minimum=0,
        ),
        state_average_target_roots=_tuple_setting(
            state_average.get("target_roots", ()),
            lambda item, label: _int_setting(item, label, minimum=0),
            "state_average.target_roots",
        ),
        state_average_equal_weights=_bool_setting(
            state_average.get("equal_weights", True),
            "state_average.equal_weights",
        ),
        active_orbital_indices=_tuple_setting(
            cas.get("active_orbital_indices", ()),
            lambda item, label: _int_setting(item, label, minimum=1),
            "cas.active_orbital_indices",
        ),
        core_orbital_indices=_tuple_setting(
            cas.get("core_orbital_indices", ()),
            lambda item, label: _int_setting(item, label, minimum=1),
            "cas.core_orbital_indices",
        ),
        orbital_source=_choice_setting(
            cas.get("orbital_source", "rhf"),
            "cas.orbital_source",
            CAS_ORBITAL_SOURCES,
        ),
        orbital_file=_string_setting(cas.get("orbital_file", ""), "cas.orbital_file"),
        target_spin=_target_spin_setting(ci.get("target_spin", "any"), "ci.target_spin"),
    )
    if settings.orbital_source == "json" and not settings.orbital_file:
        raise ValueError("cas.orbital_file is required when cas.orbital_source=json")
    return _validate_ci_settings(_validate_casci_settings(settings), ci_section="ci")


def settings_from_casscf_scaffold_config(config: dict) -> FCISettings:
    """Return CASCI backend settings for the explicit zero-update CASSCF scaffold."""
    casscf = _config_section(config, "casscf", "CASSCF")
    raw_max_macro_iterations = casscf.get("max_macro_iterations", -1)
    if isinstance(raw_max_macro_iterations, bool):
        raise ValueError(
            "method=casscf requires [casscf] max_macro_iterations as integer 0, not a boolean."
        )
    max_macro_iterations = _int_setting(
        raw_max_macro_iterations,
        "casscf.max_macro_iterations",
    )
    if max_macro_iterations != 0:
        raise ValueError(
            "method=casscf is only available as a zero-orbital-update scaffold; "
            "set [casscf] max_macro_iterations=0 or use method=casci."
        )
    return settings_from_casci_config(config)


def _active_nelec_from_setting(
    active_electrons: int | tuple[int, int] | list[int],
    nalpha: int,
    nbeta: int,
    frozen_core: int,
) -> tuple[tuple[int, int], int]:
    if isinstance(active_electrons, (list, tuple)):
        active_nelec = _as_nelec_pair(active_electrons)
        inactive_alpha = int(nalpha) - active_nelec[0]
        inactive_beta = int(nbeta) - active_nelec[1]
        if inactive_alpha < 0 or inactive_beta < 0 or inactive_alpha != inactive_beta:
            raise ValueError("active_electrons must leave a non-negative closed-shell core")
        inferred_core = inactive_alpha
        if frozen_core not in (0, inferred_core):
            raise ValueError("frozen_core is inconsistent with active_electrons")
        return active_nelec, inferred_core

    active_electron_total = int(active_electrons)
    if active_electron_total:
        inactive_electrons = int(nalpha) + int(nbeta) - active_electron_total
        if inactive_electrons < 0 or inactive_electrons % 2:
            raise ValueError("active_electrons must leave a non-negative closed-shell core")
        inferred_core = inactive_electrons // 2
        if frozen_core not in (0, inferred_core):
            raise ValueError("frozen_core is inconsistent with active_electrons")
        frozen_core = inferred_core

    active_nelec = (int(nalpha) - int(frozen_core), int(nbeta) - int(frozen_core))
    if active_nelec[0] < 0 or active_nelec[1] < 0:
        raise ValueError("active_electrons must leave a non-negative closed-shell core")
    return active_nelec, int(frozen_core)


def _active_space(
    h1e: np.ndarray,
    eri: np.ndarray,
    nelec: tuple[int, int],
    ecore: float,
    settings: FCISettings,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int], float, dict[str, int | str]]:
    h1e = _real_array(h1e, "h1e")
    eri = _real_array(eri, "eri")
    if h1e.ndim != 2 or h1e.shape[0] != h1e.shape[1]:
        raise ValueError("h1e must be a square matrix")
    norb = h1e.shape[0]
    if eri.shape != (norb, norb, norb, norb):
        raise ValueError("eri must have shape (norb, norb, norb, norb)")

    nalpha, nbeta = _as_nelec_pair(nelec)
    ecore = _float_setting(ecore, "ecore")
    explicit_active = tuple(settings.active_orbital_indices or ())
    explicit_core = tuple(settings.core_orbital_indices or ())

    if explicit_active or explicit_core:
        active_list, core_list = active_core_lists(
            {"cas": {"active_orbital_indices": explicit_active, "core_orbital_indices": explicit_core}},
            norb,
        )
        if not active_list:
            raise ValueError("Explicit CAS selection requires cas.active_orbital_indices")
        frozen_core = len(core_list)
        if settings.frozen_core and settings.frozen_core != frozen_core:
            raise ValueError("frozen_core must match the length of core_orbital_indices")
        active_norb = len(active_list)
        active_nelec, inferred_core = _active_nelec_from_setting(
            settings.active_electrons,
            nalpha,
            nbeta,
            frozen_core,
        )
        if inferred_core != frozen_core:
            raise ValueError("active_electrons is inconsistent with core_orbital_indices")
        if settings.active_orbitals and settings.active_orbitals != active_norb:
            raise ValueError("active_orbitals must match the length of active_orbital_indices")
        if active_norb < max(active_nelec):
            raise ValueError("active_orbital_indices selects too few orbitals for the active electron count")

        active = np.array(active_list, dtype=int)
        core = np.array(core_list, dtype=int)
        h_active = h1e[np.ix_(active, active)].copy()
        eri_active = eri[np.ix_(active, active, active, active)].copy()
        ecore_active = float(ecore)
        if frozen_core:
            for i in core:
                ecore_active += 2.0 * h1e[i, i]
            for i in core:
                for j in core:
                    ecore_active += 2.0 * eri[i, i, j, j] - eri[i, j, j, i]
            for p, pp in enumerate(active):
                for q, qq in enumerate(active):
                    for i in core:
                        h_active[p, q] += 2.0 * eri[pp, qq, i, i] - eri[pp, i, i, qq]
        selection = "explicit"
    else:
        frozen_core = settings.frozen_core
        active_nelec, frozen_core = _active_nelec_from_setting(
            settings.active_electrons,
            nalpha,
            nbeta,
            frozen_core,
        )

        if frozen_core < 0 or frozen_core > min(nalpha, nbeta):
            raise ValueError("frozen_core is incompatible with the electron count")

        active_norb = settings.active_orbitals or (norb - frozen_core)
        if active_norb < max(active_nelec):
            raise ValueError("active_orbitals is too small for the active electron count")
        if frozen_core + active_norb > norb:
            raise ValueError("active_orbitals extends beyond the available MO space")

        core = range(frozen_core)
        active = slice(frozen_core, frozen_core + active_norb)
        h_active = h1e[active, active].copy()
        eri_active = eri[active, active, active, active].copy()
        ecore_active = float(ecore)

        if frozen_core:
            for i in core:
                ecore_active += 2.0 * h1e[i, i]
            for i in core:
                for j in core:
                    ecore_active += 2.0 * eri[i, i, j, j] - eri[i, j, j, i]
            for p in range(active_norb):
                pp = p + frozen_core
                for q in range(active_norb):
                    qq = q + frozen_core
                    for i in core:
                        h_active[p, q] += 2.0 * eri[pp, qq, i, i] - eri[pp, i, i, qq]
        active_list = list(range(frozen_core, frozen_core + active_norb))
        core_list = list(range(frozen_core))
        selection = "sequential"

    metadata = {
        "norb": norb,
        "active_orbitals": active_norb,
        "active_electrons": sum(active_nelec),
        "active_alpha_electrons": active_nelec[0],
        "active_beta_electrons": active_nelec[1],
        "frozen_core": frozen_core,
        "determinants": comb(active_norb, active_nelec[0]) * comb(active_norb, active_nelec[1]),
        "active_orbital_indices": ",".join(str(i + 1) for i in active_list),
        "core_orbital_indices": ",".join(str(i + 1) for i in core_list),
        "orbital_selection": selection,
        "orbital_source": settings.orbital_source,
    }
    return h_active, eri_active, active_nelec, ecore_active, metadata


class FCI:
    """FCI energy calculation using OpenQP-native RHF orbitals and integrals."""

    method_label = "fci"
    data_prefix = "FCI"
    log_title = "PyOQP: Full Configuration Interaction"
    active_section = "[fci]"
    ci_section = "[fci]"

    def __init__(self, mol):
        self.mol = mol
        self.settings = _settings_from_config(mol.config)

    def energy(self, ref_energy=None):
        if self.settings.integral_backend != "native":
            raise ValueError(f"Only {self.method_label}.integral_backend=native is implemented")
        if self.mol.config["input"]["functional"]:
            raise ValueError(f"{self.data_prefix} requires HF integrals; unset input.functional")
        if self.mol.config["scf"]["type"] != "rhf":
            raise ValueError(f"{self.data_prefix} MVP supports only closed-shell RHF references")
        if self.mol.data["nelec_A"] != self.mol.data["nelec_B"]:
            raise ValueError(f"{self.data_prefix} MVP supports only closed-shell singlets")

        h1e, eri, nelec, ecore, metadata = self._native_mo_integrals()
        determinants = _determinants(h1e.shape[0], nelec)
        target_multiplicity = _target_spin_multiplicity(self.settings.target_spin)
        determinant_count = int(metadata["determinants"])

        def solve_and_diagnose(window_nroot: int):
            davidson_subspace = self.settings.davidson_subspace
            if davidson_subspace:
                davidson_subspace = max(davidson_subspace, int(window_nroot))
            window_energies, window_coeffs = solve_fci(
                h1e,
                eri,
                nelec,
                ecore=ecore,
                nroot=window_nroot,
                max_det=self.settings.max_det,
                max_memory=self.settings.max_memory,
                eig_tol=self.settings.eig_tol,
                integral_cutoff=self.settings.integral_cutoff,
                solver=self.settings.solver,
                davidson_maxiter=self.settings.davidson_maxiter,
                davidson_subspace=davidson_subspace,
                active_section=self.active_section,
                ci_section=self.ci_section,
            )
            window_s2, window_multiplicity = fci_spin_diagnostics(
                window_coeffs,
                determinants,
                h1e.shape[0],
                nelec,
            )
            return window_energies, window_coeffs, window_s2, window_multiplicity

        if target_multiplicity is None:
            energies, coeffs, s2, multiplicity = solve_and_diagnose(self.settings.nroot)
            root_indices = np.arange(len(energies), dtype=np.int64)
        else:
            solve_nroot = min(determinant_count, max(1, int(self.settings.nroot)))
            while True:
                energies, coeffs, s2, multiplicity = solve_and_diagnose(solve_nroot)
                try:
                    energies, coeffs, s2, multiplicity, root_indices = _filter_roots_by_target_spin(
                        energies,
                        coeffs,
                        s2,
                        multiplicity,
                        target_spin=self.settings.target_spin,
                        requested_nroot=self.settings.nroot,
                        ci_label=self.data_prefix,
                        ci_section=self.ci_section,
                    )
                    break
                except ValueError:
                    if solve_nroot >= determinant_count:
                        raise
                    solve_nroot = min(determinant_count, max(solve_nroot + 1, 2 * solve_nroot))

        self.mol.energies = energies.tolist()
        self.mol.data[f"OQP::{self.data_prefix}_ENERGIES"] = np.ascontiguousarray(energies, dtype=np.float64)
        self.mol.data[f"OQP::{self.data_prefix}_CI_VECTORS"] = np.ascontiguousarray(coeffs, dtype=np.float64)
        self.mol.data[f"OQP::{self.data_prefix}_S2"] = np.ascontiguousarray(s2, dtype=np.float64)
        self.mol.data[f"OQP::{self.data_prefix}_MULTIPLICITY"] = np.ascontiguousarray(
            multiplicity,
            dtype=np.int64,
        )
        self.mol.data[f"OQP::{self.data_prefix}_DET_COUNT"] = np.ascontiguousarray(
            [metadata["determinants"]],
            dtype=np.int64,
        )
        self.mol.data[f"OQP::{self.data_prefix}_ROOT_INDICES"] = np.ascontiguousarray(
            root_indices,
            dtype=np.int64,
        )
        state_average_info = _state_average_payload(self.settings, energies, root_indices)
        if state_average_info is not None:
            prefix = self.method_label.upper()
            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ENERGY"] = np.ascontiguousarray(
                [state_average_info["energy"]],
                dtype=np.float64,
            )
            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_WEIGHTS"] = state_average_info["weights"]
            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ROOTS"] = state_average_info["roots"]
            self.mol.data[f"OQP::{prefix}_STATE_AVERAGE_ROOT_INDICES"] = state_average_info[
                "root_indices"
            ]
        if self.settings.save_ci_vectors:
            active_alpha_electrons = int(metadata.get("active_alpha_electrons", nelec[0]))
            active_beta_electrons = int(metadata.get("active_beta_electrons", nelec[1]))
            artifact_arrays = {
                "energies": np.ascontiguousarray(energies, dtype=np.float64),
                "ci_vectors": np.ascontiguousarray(coeffs, dtype=np.float64),
                "s2": np.ascontiguousarray(s2, dtype=np.float64),
                "multiplicity": np.ascontiguousarray(multiplicity, dtype=np.int64),
                "root_indices": np.ascontiguousarray(root_indices, dtype=np.int64),
                "determinant_count": np.asarray([metadata["determinants"]], dtype=np.int64),
                "active_electrons": np.asarray(
                    [
                        active_alpha_electrons,
                        active_beta_electrons,
                    ],
                    dtype=np.int64,
                ),
                "active_orbital_indices": np.asarray(str(metadata.get("active_orbital_indices", ""))),
                "core_orbital_indices": np.asarray(str(metadata.get("core_orbital_indices", ""))),
                "target_spin": np.asarray(str(self.settings.target_spin)),
            }
            if state_average_info is not None:
                artifact_arrays.update(
                    {
                        "state_average_energy": np.asarray(
                            [state_average_info["energy"]],
                            dtype=np.float64,
                        ),
                        "state_average_weights": state_average_info["weights"],
                        "state_average_roots": state_average_info["roots"],
                        "state_average_root_indices": state_average_info["root_indices"],
                    }
                )
            _save_npz_artifact(self.mol, f"{self.method_label}_ci_vectors", **artifact_arrays)
        if self.settings.save_rdm and self.method_label == "fci":
            rdm1_roots, rdm2_roots, natural_occupation_roots = self._store_native_rdm_tensors(
                coeffs,
                determinants,
                h1e.shape[0],
            )
            self._save_native_rdm_artifact(
                metadata,
                rdm1_roots,
                rdm2_roots,
                natural_occupation_roots,
                energies,
                root_indices,
                nelec,
            )
        self.mol.mol_energy.energy = float(
            state_average_info["energy"] if state_average_info is not None else energies[0]
        )
        ci_vector_log = None
        if self.settings.print_ci_vectors:
            ci_vector_log = _ci_vector_log_entries(
                coeffs,
                determinants,
                h1e.shape[0],
                threshold=self.settings.ci_print_threshold,
                active_orbital_indices=str(metadata.get("active_orbital_indices", "")),
                root_indices=root_indices,
            )

        print_module_banner(
            self.mol, self.data_prefix,
            str(self.log_title).replace("PyOQP:", "").strip())
        dump_log(
            self.mol,
            title=self.log_title,
            section="fci",
            info={
                "method": self.method_label,
                "ci_label": self.data_prefix,
                "hf_energy": ref_energy[0] if ref_energy else None,
                "energies": self.mol.energies,
                "s2": s2.tolist(),
                "multiplicity": multiplicity.tolist(),
                "root_indices": root_indices.tolist(),
                "target_spin": self.settings.target_spin,
                "state_average": None
                if state_average_info is None
                else {
                    "energy": state_average_info["energy"],
                    "weights": state_average_info["weights"].tolist(),
                    "roots": state_average_info["roots"].tolist(),
                    "root_indices": state_average_info["root_indices"].tolist(),
                },
                "ci_vector_log": ci_vector_log,
                **metadata,
            },
        )
        return self.mol.energies

    def _store_native_rdm_tensors(self, ci_vectors, determinants, norb: int):
        from oqp.library.rdm import (
            make_rdm1_spatial,
            make_rdm2_spatial,
            natural_orbital_occupations,
        )

        coeffs = _real_array(ci_vectors, f"{self.data_prefix} CI vectors")
        if coeffs.ndim == 1:
            coeffs = coeffs[:, None]
        if coeffs.ndim != 2 or coeffs.shape[0] != len(determinants):
            raise ValueError("CI vectors must have shape (ndet, nroot)")
        rdm1_roots = np.stack(
            [
                make_rdm1_spatial(coeffs[:, root], determinants, norb)
                for root in range(coeffs.shape[1])
            ],
            axis=0,
        )
        rdm2_roots = np.stack(
            [
                make_rdm2_spatial(coeffs[:, root], determinants, norb)
                for root in range(coeffs.shape[1])
            ],
            axis=0,
        )
        natural_occupation_roots = np.stack(
            [natural_orbital_occupations(rdm1_roots[root]) for root in range(coeffs.shape[1])],
            axis=0,
        )
        prefix = self.method_label.upper()
        self.mol.data[f"OQP::{prefix}_RDM1"] = np.ascontiguousarray(
            rdm1_roots[0],
            dtype=np.float64,
        )
        self.mol.data[f"OQP::{prefix}_RDM2"] = np.ascontiguousarray(
            rdm2_roots[0],
            dtype=np.float64,
        )
        self.mol.data[f"OQP::{prefix}_RDM1_ROOTS"] = np.ascontiguousarray(
            rdm1_roots,
            dtype=np.float64,
        )
        self.mol.data[f"OQP::{prefix}_RDM2_ROOTS"] = np.ascontiguousarray(
            rdm2_roots,
            dtype=np.float64,
        )
        self.mol.data[f"OQP::{prefix}_NATURAL_OCCUPATIONS"] = np.ascontiguousarray(
            natural_occupation_roots[0],
            dtype=np.float64,
        )
        self.mol.data[f"OQP::{prefix}_NATURAL_OCCUPATIONS_ROOTS"] = np.ascontiguousarray(
            natural_occupation_roots,
            dtype=np.float64,
        )
        return rdm1_roots, rdm2_roots, natural_occupation_roots

    def _save_native_rdm_artifact(
        self,
        metadata,
        rdm1_roots,
        rdm2_roots,
        natural_occupation_roots,
        energies,
        root_indices,
        active_nelec,
    ) -> None:
        active_alpha_electrons = int(metadata.get("active_alpha_electrons", active_nelec[0]))
        active_beta_electrons = int(metadata.get("active_beta_electrons", active_nelec[1]))
        _save_npz_artifact(
            self.mol,
            f"{self.method_label}_rdm",
            rdm1=np.ascontiguousarray(rdm1_roots[0], dtype=np.float64),
            rdm2=np.ascontiguousarray(rdm2_roots[0], dtype=np.float64),
            rdm1_roots=np.ascontiguousarray(rdm1_roots, dtype=np.float64),
            rdm2_roots=np.ascontiguousarray(rdm2_roots, dtype=np.float64),
            natural_occupations=np.ascontiguousarray(
                natural_occupation_roots[0],
                dtype=np.float64,
            ),
            natural_occupation_roots=np.ascontiguousarray(
                natural_occupation_roots,
                dtype=np.float64,
            ),
            energies=np.ascontiguousarray(energies, dtype=np.float64),
            root_indices=np.ascontiguousarray(root_indices, dtype=np.int64),
            active_electrons=np.asarray(
                [
                    active_alpha_electrons,
                    active_beta_electrons,
                ],
                dtype=np.int64,
            ),
            active_orbital_indices=np.asarray(str(metadata.get("active_orbital_indices", ""))),
            core_orbital_indices=np.asarray(str(metadata.get("core_orbital_indices", ""))),
            target_spin=np.asarray(str(self.settings.target_spin)),
        )

    def _check_ao_eri_budget(self, nbf: int) -> None:
        """Guard the dense AO ERI allocation (nbf**4 doubles) before building it."""
        dense_bytes = 8 * int(nbf) ** 4
        budget_bytes = max(1, int(self.settings.max_memory)) * 1024 * 1024
        if dense_bytes > budget_bytes:
            raise ValueError(
                f"Dense AO ERI for nbf={nbf} needs ~{dense_bytes / 1024 ** 3:.2f} GiB, "
                f"exceeding the {self.active_section} max_memory budget of {self.settings.max_memory} MiB. "
                "Dense FCI is only intended for small basis sets."
            )

    def _native_mo_integrals(self):
        nbf = int(self.mol.data.get_basis()["nbf"])
        self._check_ao_eri_budget(nbf)
        oqp.fci_ao_integrals(self.mol)

        hcore = _unpack_lower_triangle(np.asarray(self.mol.data["OQP::Hcore"], dtype=float), nbf)
        coeff = np.asarray(self.mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
        eri_ao = np.asarray(self.mol.data["OQP::AO_ERI"], dtype=float).reshape(
            (nbf, nbf, nbf, nbf),
            order="F",
        )

        h1e, eri = _transform_integrals(hcore, eri_ao, coeff)
        nelec = (int(self.mol.data["nelec_A"]), int(self.mol.data["nelec_B"]))
        ecore = float(self.mol.mol_energy.nenergy)
        result = _active_space(h1e, eri, nelec, ecore, self.settings)
        result[4]["orbital_source"] = "rhf"
        return result
