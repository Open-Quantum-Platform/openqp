"""Small dense Full Configuration Interaction energy driver."""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from itertools import combinations
import math
import warnings
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
    #: Spatial-symmetry filter on the returned roots. "any" (the default)
    #: leaves root selection exactly as it was; a point-group irrep name keeps
    #: only roots whose dominant irrep is that one.
    target_irrep: str = "any"
    #: Weight a root must carry in its dominant irrep to count as a symmetry
    #: eigenstate at all.
    irrep_min_purity: float = 0.5


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
    array = _as_f64c(raw)
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
    """Write one .npz artifact, from the MPI world root only.

    Every rank reaches this with the same destination path, so an unguarded
    write means concurrent truncation of the same ZIP -- a nondeterministic or
    corrupt artifact rather than a clean one.  `world_rank`, not `rank`: inside
    an MPIManager.task_groups split every group root has rank == 0 and they all
    share this path, exactly as `mpi_dump` guards the shared log.
    """
    from oqp.utils.mpi_utils import MPIManager
    if MPIManager().world_rank != 0 and getattr(mol, "usempi", False):
        return ""
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
    flipped = _apply_spin_flip(coeff, dets, norb, det_index)
    return float(s2 + float(coeff @ flipped))


def _apply_spin_flip(coeff, dets, norb, det_index):
    """Return the spin-flip part of ``S^2`` applied to one CI vector.

    ``<a|S^2|b> = ms(ms+1) <a|b> + a . _apply_spin_flip(b)``, so factoring the
    application out of :func:`compute_s2` is what makes the OFF-diagonal
    elements reachable. Only the diagonal was ever needed while a root's
    multiplicity was read from its own expectation value; resolving a
    spin-mixed degenerate subspace needs the whole matrix.
    """
    out = np.zeros(len(dets), dtype=float)
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
                    out[row] += (
                        phase_q * phase_alpha_q * phase_p * phase_beta_p * c_col
                    )
    return out


def s2_matrix(vectors, determinants, norb, nelec):
    """``<a|S^2|b>`` over the columns of ``vectors``, in the determinant basis."""
    dets = list(determinants)
    vecs = _real_array(vectors, "CI vectors")
    if vecs.ndim == 1:
        vecs = vecs[:, None]
    if vecs.shape[0] != len(dets):
        raise ValueError("CI vectors must have one row per determinant")
    nalpha, nbeta = _as_nelec_pair(nelec)
    ms = 0.5 * (nalpha - nbeta)
    det_index = {det: idx for idx, det in enumerate(dets)}
    applied = np.column_stack([
        _apply_spin_flip(vecs[:, k], dets, norb, det_index)
        for k in range(vecs.shape[1])
    ])
    mat = vecs.T @ applied + ms * (ms + 1.0) * (vecs.T @ vecs)
    # S^2 is Hermitian; symmetrise away the round-off so the eigensolver
    # cannot return complex pairs from an asymmetric input.
    return 0.5 * (mat + mat.T)


#: Energies closer than this are treated as one degenerate cluster. A CI solve
#: converges eigenpairs to eig_tol (1e-10 by default), so 1e-8 is loose enough
#: to catch a genuine degeneracy resolved only to solver precision and tight
#: enough not to merge two physically distinct states.
DEGENERACY_TOLERANCE = 1.0e-8


def degenerate_clusters(energies, tol=DEGENERACY_TOLERANCE):
    """Indices of the energy-ordered roots, grouped into degenerate clusters."""
    e = np.atleast_1d(np.asarray(energies, dtype=float))
    clusters = []
    current = [0] if e.size else []
    for k in range(1, e.size):
        if abs(e[k] - e[current[0]]) <= tol:
            current.append(k)
        else:
            clusters.append(current)
            current = [k]
    if current:
        clusters.append(current)
    return clusters


def spin_purify_degenerate_clusters(energies, coeffs, determinants, norb, nelec,
                                    *, tol=DEGENERACY_TOLERANCE):
    """Re-express each degenerate cluster in the S^2 eigenbasis.

    Within a degenerate -- or numerically near-degenerate -- energy subspace the
    eigensolver is free to return any orthogonal mixture of the degenerate
    states, and a mixture of different-spin states is not an S^2 eigenstate at
    all. Its <S^2> is then a weighted average, and the multiplicity read from it
    can be one no electron count admits.

    Rotating WITHIN a degenerate subspace leaves every vector an eigenvector of
    H with the same eigenvalue, so this changes no energy. What it changes is
    which combination each root is -- and therefore its <S^2>, its multiplicity
    label, and whether a target_spin filter picks it.

    Returns ``(coeffs, s2, multiplicity, changed)``. ``changed`` is False when
    every cluster was already spin-pure, in which case the vectors are returned
    untouched rather than re-expressed in a numerically equivalent basis.
    """
    e = np.atleast_1d(np.asarray(energies, dtype=float))
    vecs = _real_array(coeffs, "CI vectors")
    if vecs.ndim == 1:
        vecs = vecs[:, None]
    dets = list(determinants)
    out = np.array(vecs, dtype=float, copy=True)
    changed = False

    for cluster in degenerate_clusters(e, tol=tol):
        if len(cluster) < 2:
            continue
        block = vecs[:, cluster]
        mat = s2_matrix(block, dets, norb, nelec)
        # Already diagonal to solver precision: leave the vectors alone rather
        # than rotate them by an arbitrary eigenvector convention.
        off = mat - np.diag(np.diag(mat))
        if np.max(np.abs(off)) <= 1.0e-8:
            continue
        vals, vecs_sub = _symmetric_eigh(mat)
        # Deterministic order inside the cluster: ascending S^2. The previous
        # order was whatever the eigensolver happened to produce, so there is
        # nothing to preserve, and a stable rule keeps references reproducible.
        order = np.argsort(vals, kind="stable")
        rotated = block @ vecs_sub[:, order]
        # A rotated eigenvector's sign is arbitrary, so pin it with the SAME
        # convention every other CI vector carries -- canonicalize_ci_phase,
        # which canonical_phase() in fci_driver.F90 mirrors. Rolling a local
        # "largest coefficient positive" here would disagree with it exactly
        # where it matters: ties are the rule in a symmetric molecule, and a
        # degenerate cluster is where symmetry-equivalent determinants carry
        # equal weight, so the tie-break is what decides the sign.
        canonicalize_ci_phase(rotated)
        out[:, cluster] = rotated
        changed = True

    # Relabel through the shared diagnostics so the engine's fast <S^2> path is
    # used where it is available; the Python reference is O(ndet * norb^2) per
    # root and this runs on every root, not just the rotated ones.
    s2, multiplicity = fci_spin_diagnostics(out, dets, norb, nelec)
    return out, s2, multiplicity, changed


def _lib_spin_square(coeffs, dets, norb, nelec):
    """Per-root ``<S^2>`` through the liboqp engine, or ``None`` when the
    symbol is unavailable and the caller should use :func:`compute_s2`.

    The Python reference is O(nroot * ndet * norb^2) with a popcount and a dict
    lookup per term -- 48 ms per root at CAS(8,8), 843 ms at CAS(10,10) -- and
    rebuilds the determinant index for every call.
    """
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    if not hasattr(lib, "fci_spin_square"):
        return None
    nalpha, nbeta = _as_nelec_pair(nelec)
    ndet = len(dets)
    nvec = int(coeffs.shape[1])
    if ndet <= 0 or nvec <= 0 or 2 * int(norb) > 62:
        return None
    # compute_s2 rejects a zero-norm root; keep that contract rather than let
    # the engine quietly report the ms-only value.
    if not np.all(np.linalg.norm(coeffs, axis=0) > 0.0):
        return None
    det_arr = _as_i64c(dets)
    civec = _as_f64c(coeffs)
    s2 = np.zeros(nvec, dtype=np.float64)
    status = int(lib.fci_spin_square(
        int(norb), ndet, nvec,
        ffi.cast("int64_t *", det_arr.ctypes.data),
        ffi.cast("double *", civec.ctypes.data),
        int(nalpha) - int(nbeta),
        ffi.cast("double *", s2.ctypes.data),
        _fci_lib_threads()))
    if status != 0:
        return None
    return s2


def fci_spin_diagnostics(
    ci_vectors: np.ndarray,
    determinants: list[int],
    norb: int,
    nelec: int | tuple[int, int] | list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-root ``(<S^2>, multiplicity)`` diagnostics.

    Routed to the liboqp engine (``fci_setup.F90``), which shares one sort of
    the determinant list across the roots; :func:`compute_s2` stays as the
    per-root Python reference and the fallback.
    """
    dets = list(determinants)
    coeffs = _real_array(ci_vectors, "CI vectors")
    if coeffs.ndim == 1:
        coeffs = coeffs[:, None]
    if coeffs.ndim != 2 or coeffs.shape[0] != len(dets):
        raise ValueError("CI vectors must have shape (ndet, nroot)")

    s2 = _lib_spin_square(coeffs, dets, norb, nelec)
    if s2 is None:
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
    return _as_f64c(s2), np.ascontiguousarray(
        multiplicity,
        dtype=np.int64,
    )


#: How far <S^2> may sit from the nearest S(S+1) before the multiplicity label
#: derived from it stops being trustworthy.  A converged spin eigenstate lands
#: on S(S+1) to solver precision; a mixture of two spin states in a degenerate
#: subspace is off by O(1), so anything in between is generous.
SPIN_LABEL_TOLERANCE = 1.0e-2


class SpinLabelAmbiguityError(RuntimeError):
    """A CI root's <S^2> is not that of any spin eigenstate.

    Deliberately NOT a ValueError: the target_spin call sites retry with more
    roots when the filter raises ValueError, and for a mixed root whose full
    degenerate manifold is already inside the window purification has run and
    failed, so more roots cannot help.  The ONE case where widening does help
    -- the window boundary cut a degenerate manifold, so purification saw a
    truncated cluster -- is recognized at the call site by the mixed root
    lying at the window's edge, and retried there explicitly.
    """


def _total_electrons(nelec) -> int:
    if isinstance(nelec, (tuple, list, np.ndarray)):
        return int(sum(int(x) for x in nelec))
    return int(nelec)


def spin_label_diagnosis(s2, multiplicity, nelec, tol: float = SPIN_LABEL_TOLERANCE):
    """Roots whose multiplicity label cannot be trusted, as (root, s2, mult, why).

    The label comes from rounding ``sqrt(1 + 4*<S^2>)``, which identifies the
    spin only when the eigenvector IS an S^2 eigenstate.  Inside a degenerate --
    or numerically near-degenerate -- energy subspace the eigensolver may return
    an arbitrary orthogonal mixture of the degenerate states, and a mixture of
    different-spin states is not an S^2 eigenstate at all.  Two things then give
    it away:

    * the implied <S^2> does not match S(S+1) for the rounded label, and
    * the label can be impossible for the electron count -- an open-shell
      two-electron Ms=0 determinant has <S^2> = 1, so the formula reports
      sqrt(5) = 2.236 -> a "doublet", which two electrons cannot form.
    """
    s2 = np.atleast_1d(np.asarray(s2, dtype=float))
    mult = np.atleast_1d(np.asarray(multiplicity, dtype=np.int64))
    total = _total_electrons(nelec)
    problems = []
    for root in range(s2.size):
        m = int(mult[root])
        value = float(s2[root])
        spin = 0.5 * (m - 1)
        pure = spin * (spin + 1.0)
        if m > total + 1:
            problems.append((root, value, m,
                             f"multiplicity {m} exceeds the maximum {total + 1} "
                             f"for {total} electrons"))
        elif (m % 2) != ((total + 1) % 2):
            problems.append((root, value, m,
                             f"multiplicity {m} has the wrong parity for "
                             f"{total} electrons (allowed: "
                             f"{'odd' if (total + 1) % 2 else 'even'})"))
        elif abs(value - pure) > tol:
            problems.append((root, value, m,
                             f"<S^2> = {value:.6f} is {abs(value - pure):.2e} "
                             f"away from S(S+1) = {pure:.6f} for multiplicity {m}"))
    return problems


def _impossible_multiplicity(multiplicity: int, nelec) -> str:
    """Why ``multiplicity`` cannot occur for this electron count, or ``''``.

    Parity and magnitude are fixed by the electron count alone, so this is
    decidable before any root is looked at.
    """
    total = _total_electrons(nelec)
    if multiplicity > total + 1:
        return (f"multiplicity {multiplicity} exceeds the maximum {total + 1} "
                f"for {total} electrons")
    if (multiplicity % 2) != ((total + 1) % 2):
        return (f"multiplicity {multiplicity} has the wrong parity for "
                f"{total} electrons (allowed: "
                f"{'odd' if (total + 1) % 2 else 'even'})")
    return ""


def _format_spin_problems(problems, ci_label: str) -> str:
    detail = "; ".join(f"root {root}: {why}" for root, _s2, _m, why in problems)
    return (
        f"{ci_label} spin labels are unreliable for {len(problems)} root(s): "
        f"{detail}. The multiplicity is read from <S^2>, which identifies the "
        f"spin only for an S^2 eigenstate; within a degenerate energy subspace "
        f"the solver may return a spin-mixed combination. Reported energies are "
        f"then the energies of mixtures, not of spin-pure states."
    )


def warn_unreliable_spin_labels(s2, multiplicity, nelec, *, ci_label: str,
                                tol: float = SPIN_LABEL_TOLERANCE) -> list:
    """Emit a warning for any root whose spin label cannot be trusted."""
    problems = spin_label_diagnosis(s2, multiplicity, nelec, tol=tol)
    if problems:
        warnings.warn(_format_spin_problems(problems, ci_label), RuntimeWarning,
                      stacklevel=2)
    return problems


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
        "weights": _as_f64c(weights),
        "roots": np.ascontiguousarray(roots, dtype=np.int64),
        "root_indices": np.ascontiguousarray(selected_root_indices, dtype=np.int64),
    }


def _state_average_real_vector(values, label: str) -> np.ndarray:
    array = _real_array(values, label)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional array")
    return _as_f64c(array)


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


def _spin_problems_can_decide_selection(
    problems,
    keep: np.ndarray,
    energies: np.ndarray,
    target_multiplicity: int,
    requested_nroot: int,
    nelec,
    *,
    degeneracy_tol: float = 1.0e-8,
) -> bool:
    """Whether an unreliable spin label could change the target-spin answer.

    A mislabelled root is harmless when it cannot enter or displace the
    selection: the selection takes the ``requested_nroot`` lowest roots whose
    label matches the target, so a mixed root strictly above that energy
    window, carrying a non-target label, decides nothing.  It does matter when

    * the target multiplicity is achievable for the electron count AND
    * its own label claims the target multiplicity (it may be selected), or
      there are not enough matching labels anyway (the missing roots may be
      hiding inside the mixtures), or it lies at or below the energy of the
      last selected root (the true ordering may differ).
    """
    total = _total_electrons(nelec)
    achievable = (target_multiplicity <= total + 1
                  and (target_multiplicity % 2) == ((total + 1) % 2))
    if not achievable:
        return False
    if keep.size < requested_nroot:
        return True
    cutoff = float(energies[keep[requested_nroot - 1]])
    scale = max(1.0, abs(cutoff))
    for root, _s2_value, label, _why in problems:
        if int(label) == int(target_multiplicity):
            return True
        if float(energies[int(root)]) <= cutoff + degeneracy_tol * scale:
            return True
    return False


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
    nelec=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return root arrays filtered to the requested spin multiplicity.

    When ``nelec`` is supplied the labels are checked before they are used to
    select anything: filtering on a label that is not a spin eigenvalue picks
    the wrong root, or reports none, without saying so.  Warn either way;
    refuse only when a label is actually about to decide the answer.
    """
    target_multiplicity = _target_spin_multiplicity(target_spin)
    root_indices = np.arange(np.asarray(energies).shape[0], dtype=np.int64)
    problems = []
    if nelec is not None:
        problems = warn_unreliable_spin_labels(
            s2, multiplicity, nelec, ci_label=ci_label)
    if target_multiplicity is None:
        return energies, coeffs, s2, multiplicity, root_indices

    # A target the electron count cannot form is answered as unsatisfiable
    # rather than as bad labels -- two electrons have no quintet -- by
    # _spin_problems_can_decide_selection declining to refuse for an
    # unachievable target, which lets the "no matching roots" error below
    # speak.
    keep = np.flatnonzero(np.asarray(multiplicity, dtype=np.int64) == target_multiplicity)
    if problems and _spin_problems_can_decide_selection(
            problems, keep, np.asarray(energies, dtype=float),
            target_multiplicity, int(requested_nroot), nelec):
        raise SpinLabelAmbiguityError(
            _format_spin_problems(problems, ci_label)
            + f" Refusing to apply {ci_section} target_spin={target_spin} "
              "to labels that are not spin eigenvalues; run with "
              "target_spin=any and select the root yourself, or lift the "
              "degeneracy."
        )
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


def _spin_orbital_integrals_reference(
    h1e: np.ndarray,
    eri: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Element-by-element spin-orbital expansion -- the numerical pin.

    Kept out of the execution path (see :func:`_spin_orbital_integrals`) and
    exercised only by the tests, which assert the two agree bit for bit.
    """
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


def _spin_orbital_integrals_numpy(
    h1e: np.ndarray,
    eri: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Strided-block expansion -- the numpy fallback for :func:`_spin_orbital_integrals`.

    ``gspin[p, q, r, s] = (p r | q s)`` whenever ``spin(p) == spin(r)`` and
    ``spin(q) == spin(s)``, and zero otherwise -- i.e. exactly four of the
    sixteen spin blocks are populated, each one the same spatial tensor with
    its second and third indices swapped.  Writing those four blocks as strided
    slice assignments replaces the ``(2*norb)**4`` Python loop of
    :func:`_spin_orbital_integrals_reference`; the result is bit-identical
    because no arithmetic happens, only a permuted copy.
    """
    norb = h1e.shape[0]
    nspin = 2 * norb
    hspin = np.zeros((nspin, nspin), dtype=float)
    hspin[:norb, :norb] = h1e
    hspin[norb:, norb:] = h1e

    # base[p, q, r, s] = eri[p, r, q, s] over spatial indices
    base = np.asarray(eri, dtype=float).transpose(0, 2, 1, 3)
    gspin = np.zeros((nspin, nspin, nspin, nspin), dtype=float)
    for pr_spin in (0, 1):          # spin shared by p and r
        pr = slice(pr_spin * norb, (pr_spin + 1) * norb)
        for qs_spin in (0, 1):      # spin shared by q and s
            qs = slice(qs_spin * norb, (qs_spin + 1) * norb)
            gspin[pr, qs, pr, qs] = base

    return hspin, gspin


def _spin_orbital_integrals(
    h1e: np.ndarray,
    eri: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Expand spatial ``h1e``/``(pq|rs)`` onto the spin-orbital basis.

    Routed to the liboqp engine (``fci_setup.F90``); the numpy build is the
    fallback when the symbol is unavailable and
    :func:`_spin_orbital_integrals_reference` is the element-by-element
    numerical pin.  All three are bit-identical -- the expansion is a permuted
    copy with no arithmetic in it.
    """
    h1e = _as_f64c(h1e)
    eri = _as_f64c(eri)
    norb = int(h1e.shape[0])
    backend = _lib_backend()
    if backend is None or norb <= 0:
        return _spin_orbital_integrals_numpy(h1e, eri)
    lib, ffi = backend
    if not hasattr(lib, "fci_spin_orbital_integrals"):
        return _spin_orbital_integrals_numpy(h1e, eri)
    nspin = 2 * norb
    # np.empty, not np.zeros: fci_spin_orbital_integrals writes every element of
    # both buffers, the zeros of the twelve unpopulated spin blocks included.
    # Zeroing here as well made the native path memset (2*norb)**4 doubles twice
    # -- once in C, once in Fortran -- for an expansion that only writes the
    # array once, which was the whole of its deficit against the numpy fallback.
    hspin = np.empty((nspin, nspin), dtype=np.float64)
    gspin = np.empty((nspin,) * 4, dtype=np.float64)
    lib.fci_spin_orbital_integrals(
        norb,
        _dblp(ffi, h1e), _dblp(ffi, eri),
        _dblp(ffi, hspin), _dblp(ffi, gspin),
        _fci_lib_threads())
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


def _apply_work_bytes_cap(max_memory: int) -> None:
    """Hand the declared ``max_memory`` budget to the native blocked kernels.

    ``[cas]``/``[fci]``/``[pt2] max_memory`` is in MiB and already governs the
    Python-side dense-vs-Davidson choice.  The string-driven sigma and RDM
    (``fci_sigma_strings.F90``) block their alpha-string range against a
    working-set ceiling that used to be a compile-time 512 MiB constant, so the
    keyword reached them not at all: a large budget could not buy bigger blocks,
    and a small one was not respected.

    Half the budget, not all of it -- the ceiling covers only the blocked
    scratch, while the same budget also has to hold the CI vectors, the RDMs and
    the folded active integrals.  A missing symbol is not an error; an older
    liboqp simply keeps its built-in default.
    """
    backend = _lib_backend()
    if backend is None:
        return
    lib, _ = backend
    setter = getattr(lib, "fci_set_work_bytes_cap", None)
    if setter is None:
        return
    try:
        budget = int(max_memory) * 1024 * 1024
    except (TypeError, ValueError):
        return
    setter(max(budget // 2, 0))


def _fci_default_threads() -> int:
    """Default team width for the native FCI/CASSCF/RDM kernels.

    These kernels are dominated by the gather/scatter passes over the
    ``(ndet, nact^2)`` working set, which are memory-bandwidth bound rather than
    flop bound.  They therefore stop scaling well before the logical CPU count,
    and *past* the saturation point extra threads make them slower, not faster.
    Measured here on a 2-socket 22-core Xeon E5-2699A v4 (88 logical CPUs),
    RDM build at CAS(12,12), and a full CAS(8,8)/cc-pVDZ CASSCF run:

        threads    1      8     22     44     88
        RDM     1.68s  1.22s  1.08s  1.57s  1.56s
        CASSCF     -   2.85s  2.49s  2.42s  3.37s

    So the useful default is neither the old hard-coded 8 nor the full 88: it is
    the physical cores of one NUMA node -- hyperthread siblings share the load/
    store units these loops saturate, and crossing the socket boundary puts the
    working set on remote memory.  Where the topology cannot be read, fall back
    to the previous behaviour rather than guess.
    """
    cpus = os.cpu_count() or 1
    try:
        with open("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list") as fh:
            siblings = fh.read().strip()
        per_core = len([s for s in siblings.replace("-", ",").split(",") if s]) or 1
        nodes = len(glob.glob("/sys/devices/system/node/node[0-9]*")) or 1
        return max(1, min(cpus, cpus // per_core // nodes))
    except OSError:
        return max(1, min(8, cpus))


def _fci_lib_threads() -> int:
    """Thread count handed to every native FCI/CASSCF/RDM kernel.

    This is the `nthreads` argument of the liboqp entry points, and it becomes
    the `num_threads(...)` clause on their OpenMP regions.  It used to be
    `min(8, cpu_count)`, which silently pinned the CI sigma, the RDM build, the
    dense Hamiltonian and the spin-orbital transform to eight threads no matter
    how wide the machine was -- on an 88-core node that left the CI hot path
    running on under a tenth of it.

    The number is resolved the same way the rest of PyOQP resolves threading, so
    one control governs the whole job:

      * ``OQP_FCI_THREADS`` -- explicit escape hatch, for pinning these kernels
        alone (benchmarking, or a machine where the gather/scatter stops scaling
        before the core count).
      * ``OMP_NUM_THREADS``, but only when it is the caller's own request --
        either already in the environment before ``import oqp``, or applied from
        ``[input] omp_threads``.  ``oqp/__init__.py`` *defaults* the variable to
        ``'1'`` when it is unset, and that default must not be mistaken for a
        request to run the CI hot path serially; ``oqp.OMP_THREADS_FROM_ENV``
        is what tells the two apart.
      * otherwise ``_fci_default_threads()``, derived from the CPU topology
        rather than hard-coded.

    A positive value is always returned: the sigma and RDM kernels read
    ``max(1, nthreads)``, so a zero would mean *serial*, not *let OpenMP decide*.
    """
    try:
        pinned = int(os.environ.get("OQP_FCI_THREADS", "0"))
    except ValueError:
        pinned = 0
    if pinned > 0:
        return pinned
    if getattr(oqp, "OMP_THREADS_FROM_ENV", False):
        try:
            requested = int(os.environ.get("OMP_NUM_THREADS", "0"))
        except ValueError:
            requested = 0
        if requested > 0:
            return requested
    return _fci_default_threads()


def _lib_dense_hamiltonian(dets, hspin, gspin, nspin, cutoff):
    """Dense determinant Hamiltonian through the liboqp OpenMP engine, or None
    when the symbol is missing (Python enumeration remains the fallback)."""
    backend = _lib_backend()
    if backend is None or nspin > 62:
        return None
    lib, ffi = backend
    if not hasattr(lib, "fci_dense_hamiltonian"):
        return None
    det_arr = _as_i64c(dets)
    ndet = int(det_arr.size)
    h = np.zeros((ndet, ndet), dtype=np.float64)
    hs = _as_f64c(hspin)
    gs = _as_f64c(gspin)
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
        self._dets = _as_i64c(dets)
        self._hs = _as_f64c(hspin)
        self._gs = _as_f64c(gspin)
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


def _determinant_index(dets, det_index=None):
    """Determinant -> row map, built only when a Python enumeration needs it.

    The native builders never look at it, so building it eagerly next to the
    determinant list was a dict of ``ndet`` entries thrown away on every CI
    solve -- 0.2 s of an 84 s H2O/cc-pVDZ CAS(6,6) CASSCF over 6185 solves,
    all of it pure waste.
    """
    if det_index is not None:
        return det_index
    return {det: idx for idx, det in enumerate(dets)}


def _build_dense_hamiltonian(dets, det_index, hspin, gspin, nspin, cutoff):
    """Assemble the explicit dense FCI Hamiltonian (8*ndet**2 bytes)."""
    lib_h = _lib_dense_hamiltonian(dets, hspin, gspin, nspin, cutoff)
    if lib_h is not None:
        return lib_h
    ndet = len(dets)
    hamiltonian = np.zeros((ndet, ndet), dtype=float)
    for row, col, value in _iter_hamiltonian_elements(
        dets, _determinant_index(dets, det_index), hspin, gspin, nspin, cutoff
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
        dets, _determinant_index(dets, det_index), hspin, gspin, nspin, cutoff
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
        _as_f64c(eigvals[order]),
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


_F64 = np.dtype(np.float64)
_I64 = np.dtype(np.int64)
_BACKEND_CACHE = None


def _as_f64c(a):
    """A C-contiguous float64 array, without paying for a no-op conversion.

    ``np.ascontiguousarray`` costs ~4.5 us even when it has nothing to do, and
    the native wrappers call it on every crossing -- the eigensolver path alone
    ran it ~2300 times per CASSCF run on arrays that were already contiguous
    float64.  Checking first is several times cheaper than calling it."""
    if type(a) is np.ndarray and a.dtype is _F64 and a.flags.c_contiguous:
        return a
    return np.ascontiguousarray(a, dtype=np.float64)


_CT_DOUBLE_P = None


def _dblp(ffi, a):
    """``a``'s buffer as a ``double *``, the cheap way.

    ``ffi.cast("double *", a.ctypes.data)`` costs ~1.08 us per array: most of it
    is ``a.ctypes``, which builds a fresh ctypes proxy object every time, and the
    rest is re-parsing the type string.  Going through ``ffi.from_buffer`` with a
    cached ctype is ~0.42 us.  That is not noise at this granularity -- the
    spin-orbital expansion casts four arrays around a 7 us kernel, so the casts
    alone were a third of what the native path cost.

    The returned pointer does not own a reference to ``a``, exactly like the
    integer ``a.ctypes.data`` it replaces, so ``a`` must stay alive across the
    native call.  Every caller here holds it in a local for that reason.
    """
    global _CT_DOUBLE_P
    if _CT_DOUBLE_P is None:
        _CT_DOUBLE_P = ffi.typeof("double *")
    return ffi.cast(_CT_DOUBLE_P, ffi.from_buffer(a))


def _as_i64c(a):
    """A C-contiguous int64 array, without paying for a no-op conversion.

    The int64 twin of :func:`_as_f64c`; determinant lists cross the boundary on
    every native call and were being run through ``asarray`` plus
    ``ascontiguousarray`` each time."""
    if type(a) is np.ndarray and a.dtype is _I64 and a.flags.c_contiguous:
        return a
    return np.ascontiguousarray(np.asarray(a, dtype=np.int64))


def _lib_backend():
    """(lib, ffi) of the loaded native core, or None when unavailable.

    Cached after the first success: this is called on every native crossing
    (877 times in a single CASSCF run) and the import machinery is not free.
    A failure is not cached, so a backend that appears later is still picked
    up."""
    global _BACKEND_CACHE
    if _BACKEND_CACHE is not None:
        return _BACKEND_CACHE
    try:
        from oqp import lib, ffi
    except Exception:
        return None
    _BACKEND_CACHE = (lib, ffi)
    return _BACKEND_CACHE


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
    # LAPACK requires LDA >= max(1, N); a 0x0 block (empty active/virtual
    # space) would trip "parameter 5 had an illegal value" on stderr before
    # returning nonzero. 1x1 is not worth a LAPACK call either.
    if n < 2:
        if n == 0:
            return (np.zeros(0, dtype=np.float64),
                    np.zeros((0, 0), dtype=np.float64))
        return (np.asarray(sym, dtype=np.float64).reshape(1).copy(),
                np.ones((1, 1), dtype=np.float64))
    a = np.array(sym, dtype=np.float64, order="C", copy=True)
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
    sym = _as_f64c(0.5 * (matrix + matrix.T))
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
        return _as_f64c(eigvals), _as_f64c(eigvecs)

    # liboqp's own ILP64 LAPACK is the primary path: it is the same threaded
    # library the rest of the native code computes with, and it is ABI-proof by
    # construction (numpy links its own, possibly LP64, LAPACK).
    lib_result = _lib_dsyevd(sym)
    if lib_result is not None:
        result = _verified(*lib_result)
        if result is not None:
            return result
    # Fall back to numpy's LAPACK when the native symbol is unavailable (e.g.
    # a pure-Python test run with no built backend), then to Jacobi.
    try:
        result = _verified(*np.linalg.eigh(sym))
        if result is not None:
            return result
    except Exception:
        pass
    return _symmetric_eigh_jacobi(sym)


class _DenseMatrixOperator:
    """``_davidson``'s operator duck type backed by an explicit dense matrix.

    Unlike :class:`_LibHamiltonianOperator` this does not re-enumerate the
    determinant Hamiltonian per application -- the matrix is already built, so
    an application is one BLAS GEMM.
    """

    __slots__ = ("_h", "shape", "_diag")

    def __init__(self, matrix):
        self._h = _as_f64c(matrix)
        self.shape = self._h.shape
        self._diag = np.ascontiguousarray(np.diag(self._h))

    def diagonal(self):
        return self._diag

    def dot(self, x):
        return self._h @ x

    def toarray(self):
        return self._h


def _lowest_dense_roots(hamiltonian, nroot, tol, max_iter):
    """Lowest ``nroot`` eigenpairs of an explicit dense symmetric matrix, or
    ``None`` to let the caller fall back to the full LAPACK solve.

    A full ``dsyevd`` produces all ``ndet`` eigenpairs, which is pure waste when
    ``nroot << ndet``: it is O(ndet^3) regardless, so a CAS(8,8) CI solve spends
    ~12 s producing 4900 eigenvectors to return one.  Davidson on the same
    already-built dense matrix returns that root in 89 ms, agreeing to 6e-14 in
    the eigenvalue and 1e-15 in the Ritz-vector overlap -- i.e. within LAPACK's
    own rounding.

    Davidson is converged an order of magnitude tighter than ``eig_tol`` so the
    caller's residual check has margin, and any failure to reach that (a hard
    tolerance, a pathological spectrum) returns ``None`` rather than raising.
    """
    ndet = int(hamiltonian.shape[0])
    # Measured crossover for the whole solve_fci dense path (nroot=1): at
    # ndet=100 the full solve already wins (0.9x), at ndet=400 it is 3.1x, at
    # 1225 16x, at 4900 62x.  Below the floor, or when the requested window is
    # a large fraction of the spectrum, take the full solve.
    if ndet < 256 or 8 * int(nroot) >= ndet:
        return None
    try:
        eigvals, eigvecs = _davidson(
            _DenseMatrixOperator(hamiltonian),
            np.ascontiguousarray(np.diag(hamiltonian)),
            int(nroot),
            tol=0.1 * float(tol),
            max_iter=int(max_iter),
        )
    except Exception:
        return None
    # _davidson can also bail out early (an exhausted correction space returns
    # the current Ritz pairs), so re-check here rather than let the caller's
    # eig_tol assertion turn a solvable case into a hard failure.
    if not (np.all(np.isfinite(eigvals)) and np.all(np.isfinite(eigvecs))):
        return None
    worst = float(np.max(np.linalg.norm(
        hamiltonian @ eigvecs - eigvecs * eigvals[None, :], axis=0)))
    if worst > float(tol):
        return None
    return eigvals, eigvecs


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
            # Normalize BEFORE orthogonalizing so the acceptance test below is
            # relative.  A preconditioned correction has norm ~||r||/gap, so an
            # absolute cutoff rejects every correction once the residual gets
            # small and the iteration silently returns unconverged: with the
            # default eig_tol=1e-10 it stalled at a residual of 3.7e-8 on a
            # 400-determinant CAS(6,6) because the correction norm had reached
            # 3.8e-9.  What matters is whether the direction survives
            # orthogonalization, which is scale-free.
            scale = np.linalg.norm(vec)
            if scale <= 0.0:
                continue
            vec = vec / scale
            for _ in range(2):  # orthogonalize twice for numerical stability
                vec = vec - current @ (current.T @ vec)
            norm = np.linalg.norm(vec)
            if norm > 1.0e-6:
                vec = vec / norm
                current = np.column_stack([current, vec])
                new_cols.append(vec)
        if not new_cols:
            # Every preconditioned correction collapsed onto the existing
            # subspace.  Returning the Ritz pairs here published them as
            # CONVERGED even with a residual above eig_tol -- and unlike the
            # dense-root path, the iterative caller does no residual check of
            # its own, so an explicit solver=davidson job (or a large
            # solver=auto one) could report an energy that violates the
            # requested tolerance without any diagnostic.  The subspace cannot
            # grow, so this is a stall, not convergence.
            if np.max(residual_norms) <= tol:
                return theta, ritz
            raise ValueError(
                f"FCI Davidson stalled: the search subspace can no longer be "
                f"expanded (every preconditioned correction is linearly "
                f"dependent on it) at a max residual of "
                f"{np.max(residual_norms):.3e} > eig_tol={tol:.3e}.  Loosen "
                f"[ci] eig_tol, raise max_subspace, or improve the guess."
            )
        sigma = np.column_stack([sigma, matvec(np.column_stack(new_cols))])
        basis = current

    raise ValueError(
        f"FCI Davidson did not converge in {max_iter} iterations "
        f"(max residual {np.max(residual_norms):.3e} > eig_tol={tol:.3e})"
    )


@dataclass(frozen=True)
class CISolveSpec:
    """Validated, resolved CI solver settings.

    One place decides what the solver is actually asked to do, so the Python
    driver and the native :func:`fci_solve` entry point cannot disagree about a
    determinant count, a solver choice or a memory budget -- and every
    user-facing message keeps being raised from here, in Python, whichever
    executes the solve.
    """

    nalpha: int
    nbeta: int
    ndet: int
    nroot: int
    solver: str
    ecore: float
    eig_tol: float
    integral_cutoff: float
    max_det: int
    max_memory: int
    davidson_maxiter: int
    davidson_subspace: int
    target_spin: str
    target_multiplicity: int | None
    budget_bytes: int


def resolve_ci_solve(
    norb: int,
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
) -> CISolveSpec:
    """Validate CI options and resolve ``solver="auto"``."""
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
    _apply_work_bytes_cap(max_memory)
    if solver == "auto":
        # The Python solver holds the (2*norb)^4 spin tensor whichever branch
        # runs, so choosing dense on the Hamiltonian alone could pick a path the
        # later guard then refuses -- auto failing a job that explicit
        # solver=davidson runs comfortably.  Weigh the combined dense working
        # set here.  This only ever steers auto TOWARD davidson; it is not a
        # rejection, so the native decline path is unaffected.
        _auto_spin = 8 * (2 * int(norb)) ** 4
        # The Python eigensolver keeps the Hamiltonian while it forms a
        # symmetric copy and then gives LAPACK a writable copy.  All three
        # ndet x ndet arrays overlap the resident spin tensor at the peak.
        solver = ("dense" if 3 * dense_bytes + _auto_spin <= budget_bytes
                  else "davidson")
    if solver == "dense" and 2 * dense_bytes > budget_bytes:
        raise ValueError(
            f"FCI dense Hamiltonian for {ndet} determinants needs "
            f"~{2 * dense_bytes / 1024 ** 3:.2f} GiB (the matrix plus the "
            f"eigenvector buffer held alongside it), exceeding the {active_section} max_memory budget "
            f"of {max_memory} MiB. Reduce the active space, raise {active_section} max_memory, "
            f"or use the iterative solver ({ci_section} solver=davidson)."
        )

    return CISolveSpec(
        nalpha=nalpha,
        nbeta=nbeta,
        ndet=ndet,
        nroot=nroot,
        solver=solver,
        ecore=ecore,
        eig_tol=eig_tol,
        integral_cutoff=integral_cutoff,
        max_det=max_det,
        max_memory=max_memory,
        davidson_maxiter=davidson_maxiter,
        davidson_subspace=davidson_subspace,
        target_spin=target_spin,
        target_multiplicity=target_multiplicity,
        budget_bytes=budget_bytes,
    )


#: Relative window that decides which amplitudes count as "the largest" when
#: :func:`canonicalize_ci_phase` picks the amplitude whose sign it fixes.
#: Symmetry-equivalent determinants carry mathematically equal weights that a
#: floating-point solve reproduces only to ~1e-14 relative, so an exact
#: ``argmax`` over ``|c|`` would let round-off elect a different representative
#: -- and therefore a different sign -- from one run to the next, which is the
#: whole freedom this convention exists to remove.  Mirrored by
#: ``FCI_PHASE_TIE_RTOL`` in ``source/modules/fci_driver.F90``.
_CI_PHASE_TIE_RTOL = 1.0e-8


def canonicalize_ci_phase(civecs: np.ndarray) -> np.ndarray:
    """Fix the arbitrary overall sign of each CI vector, in place.

    An eigenvector is determined only up to ``|Psi> -> -|Psi>``, and which of
    the two a diagonalization hands back is not a property of the calculation:
    it follows the Davidson start vector, the LAPACK implementation, and the
    order the OpenMP reductions happen to accumulate in.  H4_MCQDPT2 returned
    ``+6.085296e-03`` for the same Heff off-diagonal on one and two threads of
    one build and ``-6.085296e-03`` on four, with the energies identical to
    1e-10.  Everything LINEAR in the CI vector inherits that freedom; the
    multistate CASPT2 effective Hamiltonian, whose off-diagonal carries the
    product of two root phases, is where it surfaced -- as a 2 x 6.085e-03
    Hartree "regression" on a calculation that had not moved.

    Convention: the amplitude of largest magnitude is positive.  Ties -- the
    rule rather than the exception in a symmetric molecule, where
    symmetry-equivalent determinants carry equal weight -- go to the lowest
    determinant index, selected within :data:`_CI_PHASE_TIE_RTOL` so that
    round-off between two mathematically equal amplitudes cannot elect a
    different representative from one run to the next.  ``canonical_phase()``
    in ``source/modules/fci_driver.F90`` applies the identical rule to the
    vectors the native engine returns, and
    ``tests/test_fci_solve.py::test_native_solve_matches_python_driver`` pins
    the two together by requiring a SIGNED overlap of +1.

    ``civecs`` is ``(ndet, nroot)``.  The convention is idempotent, so applying
    it to already-canonical vectors is a no-op.
    """
    vecs = np.asarray(civecs)
    if vecs.ndim != 2 or vecs.size == 0:
        return civecs
    if not vecs.flags.writeable:
        vecs = vecs.copy()
    magnitudes = np.abs(vecs)
    # A wholly zero column falls out of the same expression: its window admits
    # every row, the leading amplitude is 0, and 0 is not negative.
    window = magnitudes.max(axis=0) * (1.0 - _CI_PHASE_TIE_RTOL)
    lead = np.argmax(magnitudes >= window, axis=0)
    flip = vecs[lead, np.arange(vecs.shape[1])] < 0.0
    if flip.any():
        vecs[:, flip] *= -1.0
    return vecs


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

    spec = resolve_ci_solve(
        norb,
        nelec,
        ecore=ecore,
        nroot=nroot,
        max_det=max_det,
        max_memory=max_memory,
        eig_tol=eig_tol,
        integral_cutoff=integral_cutoff,
        solver=solver,
        davidson_maxiter=davidson_maxiter,
        davidson_subspace=davidson_subspace,
        target_spin=target_spin,
        active_section=active_section,
        ci_section=ci_section,
    )
    nalpha, nbeta = spec.nalpha, spec.nbeta
    ndet = spec.ndet
    nroot = spec.nroot
    solver = spec.solver
    ecore = spec.ecore
    eig_tol = spec.eig_tol
    integral_cutoff = spec.integral_cutoff
    davidson_maxiter = spec.davidson_maxiter
    davidson_subspace = spec.davidson_subspace
    target_multiplicity = spec.target_multiplicity
    budget_bytes = spec.budget_bytes

    dets = _determinants(norb, (nalpha, nbeta))
    # The determinant -> row map is only consulted by the Python enumeration
    # fallbacks, so it is left to _determinant_index to build on demand.
    det_index = None
    nspin = 2 * norb
    # This tensor is (2*norb)**4 doubles and does NOT depend on the determinant
    # count, so the ndet x ndet budget above can pass while this allocation is
    # what actually exhausts memory: 2 electrons in 50 orbitals is only 2500
    # determinants (~48 MiB Hamiltonian, inside max_memory=128) while gspin
    # alone is ~763 MiB.  It is built here for BOTH the dense and the Davidson
    # branch below, so it is checked once, here -- at the allocation, not in
    # resolve_ci_solve, because the native driver never reaches this code and
    # must stay free to accept a wide active space under a small budget.
    _spin_bytes = 8 * nspin ** 4
    # gspin stays LIVE while the dense Hamiltonian is built below, so the two
    # have to be weighed together, not one at a time: 2 electrons in 50
    # orbitals clears a ~763 MiB spin tensor and a ~48 MiB Hamiltonian
    # separately under max_memory=768, while their sum does not.  (50 orbitals
    # also forces this Python fallback, since the native packed-determinant
    # driver caps at 31.)
    # _symmetric_eigh retains the Hamiltonian while forming a symmetric copy
    # and a writable LAPACK buffer.  Counting only the Hamiltonian allowed a
    # nominally in-budget dense solve to exceed its ceiling by two full
    # determinant matrices.
    _dense_peak_bytes = 3 * 8 * ndet * ndet if solver == "dense" else 0
    _live_bytes = _spin_bytes + _dense_peak_bytes
    if _live_bytes > budget_bytes:
        _extra = ("" if solver != "dense" else
                  f" plus the {_dense_peak_bytes / 1024 ** 3:.2f} GiB dense "
                  f"matrix working set held alongside it")
        raise ValueError(
            f"The Python CI solver needs ~{_spin_bytes / 1024 ** 3:.2f} GiB "
            f"for the spin-orbital integral tensor at norb={norb}{_extra} "
            f"-- ~{_live_bytes / 1024 ** 3:.2f} GiB live at once, exceeding "
            f"the configured max_memory budget of "
            f"{budget_bytes / 1024 ** 2:.0f} MiB.  The spin tensor does not "
            f"depend on the determinant count, so fewer electrons will not "
            f"help: reduce active_orbitals, raise max_memory, or use "
            f"solver=davidson to drop the dense Hamiltonian."
        )
    hspin, gspin = _spin_orbital_integrals(h1e, eri)

    if solver == "dense":
        hamiltonian = _build_dense_hamiltonian(
            dets, det_index, hspin, gspin, nspin, integral_cutoff
        )
        eigvals = eigvecs = None
        if target_multiplicity is None:
            # Only the lowest `nroot` are ever consumed here, so a full
            # ndet-root dsyevd is waste; Davidson on the already-built dense
            # matrix returns them to LAPACK accuracy for a fraction of the cost
            # (None falls back to the full solve below).
            lowest = _lowest_dense_roots(
                hamiltonian, nroot, eig_tol, davidson_maxiter
            )
            if lowest is not None:
                eigvals, eigvecs = lowest
        if eigvals is None:
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
                        nelec=(nalpha, nbeta),
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
            # gspin was built above and stays resident through the solve, so
            # the workspace cannot be weighed against the whole budget -- the
            # native Davidson guard already adds it.
            work_bytes = (8 * ndet * (2 * max_subspace + 4 * solve_nroot)
                          + 8 * nspin ** 4)
            if work_bytes > budget_bytes:
                raise ValueError(
                    f"FCI Davidson working set for {ndet} determinants needs "
                    f"~{work_bytes / 1024 ** 3:.2f} GiB including the resident "
                    f"spin-orbital tensor, exceeding the {active_section} max_memory "
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
                        nelec=(nalpha, nbeta),
                    )
                )
                break
            except ValueError:
                if solve_nroot >= ndet:
                    raise
                solve_nroot = min(ndet, max(solve_nroot + 1, 2 * solve_nroot))

    return selected_eigvals + ecore, canonicalize_ci_phase(selected_eigvecs)


def _unpack_lower_triangle(packed: np.ndarray, n: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            matrix[i, j] = packed[idx]
            matrix[j, i] = packed[idx]
            idx += 1
    return matrix


def _transform_integrals_reference(
    hcore_ao: np.ndarray,
    eri_ao: np.ndarray,
    coeff: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Single-einsum AO -> MO transformation -- the numerical pin.

    Out of the execution path (see :func:`_transform_integrals`) but kept in
    the tree as the reference the Fortran engine is validated against, and as
    the fallback when the native symbol is unavailable.
    """
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


def _transform_integrals(
    hcore_ao: np.ndarray,
    eri_ao: np.ndarray,
    coeff: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """AO -> MO transformation of the core Hamiltonian and the ERI tensor.

    Routed to the liboqp engine (``mo_transform.F90``), which does the four
    quarter transformations as four DGEMMs with transposed output so the index
    roll between them is free.  :func:`_transform_integrals_reference` is the
    numpy pin and the fallback.

    Note what this port does and does not buy.  The einsum was never the naive
    O(nbf^8) contraction -- ``optimize=True`` already reduces it to a GEMM
    chain running at 46-72 GFlop/s -- so this is an architectural move (compute
    belongs in liboqp) rather than a rescue.  A *naive* explicit chain in numpy
    is in fact slower than the einsum (0.48-0.78x at nbf=24..58) because the
    ``ascontiguousarray`` rolls between steps cost more than the einsum's own
    bookkeeping; avoiding those copies is what the Fortran formulation is for.
    """
    hcore_ao = _as_f64c(hcore_ao)
    eri_ao = _as_f64c(eri_ao)
    coeff = _as_f64c(coeff)
    nbf = int(coeff.shape[0])
    backend = _lib_backend()
    if (backend is None or nbf <= 0
            or not hasattr(backend[0], "mo_transform_eri")):
        return _transform_integrals_reference(hcore_ao, eri_ao, coeff)
    lib, ffi = backend

    h1e = np.zeros((nbf, nbf), dtype=np.float64)
    lib.mo_transform_h1e(
        nbf,
        ffi.cast("double *", hcore_ao.ctypes.data),
        ffi.cast("double *", coeff.ctypes.data),
        ffi.cast("double *", h1e.ctypes.data))

    eri = np.zeros((nbf,) * 4, dtype=np.float64)
    status = int(lib.mo_transform_eri(
        nbf,
        ffi.cast("double *", eri_ao.ctypes.data),
        ffi.cast("double *", coeff.ctypes.data),
        ffi.cast("double *", eri.ctypes.data)))
    if status != 0:
        # the engine could not allocate its nbf^4 work buffer
        return _transform_integrals_reference(hcore_ao, eri_ao, coeff)
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
        target_irrep=str(raw.get("irrep", "any")).strip() or "any",
        irrep_min_purity=float(raw.get("irrep_min_purity", 0.5)),
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
        target_irrep=str(ci.get("irrep", "any")).strip() or "any",
        irrep_min_purity=float(ci.get("irrep_min_purity", 0.5)),
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


def _fold_core_reference(h1e, eri, active_list, core_list, ecore=0.0):
    """Frozen-core fold -- the numerical pin.

    Out of the execution path (see :func:`_fold_core`).  The accumulation order
    here is the contract: the inactive energy is summed into the caller's
    ``ecore`` in the order (2 h_ii), then (2 (ii|jj) - (ij|ji)), and each
    exchange-corrected term is formed before it is added.  The engine
    reproduces that order so the two agree bit for bit.
    """
    active = np.asarray(active_list, dtype=int)
    h_active = h1e[np.ix_(active, active)].copy()
    ecore_active = float(ecore)
    if not len(core_list):
        return h_active, ecore_active
    for i in core_list:
        ecore_active += 2.0 * h1e[i, i]
    for i in core_list:
        for j in core_list:
            ecore_active += 2.0 * eri[i, i, j, j] - eri[i, j, j, i]
    for p, pp in enumerate(active):
        for q, qq in enumerate(active):
            for i in core_list:
                h_active[p, q] += 2.0 * eri[pp, qq, i, i] - eri[pp, i, i, qq]
    return h_active, ecore_active


def _fold_core(h1e, eri, active_list, core_list, ecore=0.0):
    """Folded active one-electron integrals and the inactive energy.

    Routed to the liboqp engine (``fci_setup.F90``); the Python loop of
    :func:`_fold_core_reference` is the pin and the fallback.  Both index lists
    are explicit so the sequential and the explicitly selected CAS paths share
    one implementation.
    """
    h1e = _as_f64c(h1e)
    eri = _as_f64c(eri)
    norb = int(h1e.shape[0])
    nact = len(active_list)
    backend = _lib_backend()
    if (backend is None or nact <= 0
            or not hasattr(backend[0], "fci_fold_core")):
        return _fold_core_reference(h1e, eri, active_list, core_list, ecore)
    lib, ffi = backend
    active = np.ascontiguousarray(np.asarray(active_list, dtype=np.int32))
    core = np.ascontiguousarray(np.asarray(core_list, dtype=np.int32))
    if core.size == 0:                       # cffi will not cast an empty buffer
        core = np.zeros(1, dtype=np.int32)
    h_active = np.zeros((nact, nact), dtype=np.float64)
    energy = np.array([float(ecore)], dtype=np.float64)
    lib.fci_fold_core(
        norb, nact, len(core_list),
        ffi.cast("int32_t *", active.ctypes.data),
        ffi.cast("int32_t *", core.ctypes.data),
        ffi.cast("double *", h1e.ctypes.data),
        ffi.cast("double *", eri.ctypes.data),
        ffi.cast("double *", h_active.ctypes.data),
        ffi.cast("double *", energy.ctypes.data))
    return h_active, float(energy[0])


@dataclass(frozen=True)
class ActiveSpacePlan:
    """Which MOs are active, which are frozen core, and how many electrons.

    Pure integer bookkeeping derived from the ``[cas]``/``[fci]`` settings, with
    no dependence on the integrals, so it is computed once and reused across
    every CI solve of a CASSCF run instead of being re-derived per
    macroiteration.  It is also exactly what the native :func:`fci_solve` needs
    to do the frozen-core fold and the active gather itself.
    """

    norb: int
    active: tuple[int, ...]
    core: tuple[int, ...]
    nelec: tuple[int, int]
    metadata: dict

    @property
    def nact(self) -> int:
        return len(self.active)

    @property
    def ncore(self) -> int:
        return len(self.core)


def active_space_plan(
    norb: int,
    nelec: tuple[int, int],
    settings: FCISettings,
) -> ActiveSpacePlan:
    """Resolve the active/core orbital selection from the settings."""
    nalpha, nbeta = _as_nelec_pair(nelec)
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
    return ActiveSpacePlan(
        norb=int(norb),
        active=tuple(int(i) for i in active_list),
        core=tuple(int(i) for i in core_list),
        nelec=(int(active_nelec[0]), int(active_nelec[1])),
        metadata=metadata,
    )


def contiguous_active_space(
    norb: int,
    nelec: tuple[int, int],
    settings,
    label: str,
) -> tuple[int, int, tuple[int, int], ActiveSpacePlan]:
    """Resolve ``(ncore, nact, active_nelec, plan)`` for the drivers that
    assume a contiguous core/active partition.

    CASSCF and the PT2 reference derived the active electron count as
    ``nelec - frozen_core`` and so ignored ``[cas] active_electrons``
    entirely: ``active_electrons=4, active_orbitals=4`` on H2O left
    ``frozen_core`` at its default and asked for CAS(10,4), which cannot hold
    the electrons; a merely *inconsistent* pair silently ran a different
    active space than the one requested.  ``active_space_plan`` is the one
    place that reconciles the two, so route through it.

    Both drivers build their orbital blocks as ``range(ncore)`` /
    ``range(ncore, ncore + nact)`` and cannot honour a scattered
    ``active_orbital_indices`` selection -- previously they did not even look
    at it.  Reject that explicitly rather than computing a different space.
    """
    if not settings.active_orbitals and not tuple(
            getattr(settings, "active_orbital_indices", ()) or ()):
        raise ValueError(f"{label} needs a valid [cas] active_orbitals / frozen_core")
    plan = active_space_plan(norb, nelec, settings)
    ncore, nact = len(plan.core), len(plan.active)
    if (tuple(plan.core) != tuple(range(ncore))
            or tuple(plan.active) != tuple(range(ncore, ncore + nact))):
        raise ValueError(
            f"{label} requires a contiguous core/active orbital partition, but "
            "[cas] active_orbital_indices / core_orbital_indices select a "
            "scattered set.  Reorder the orbitals ([cas] sort_orbitals or an "
            "orbital file) so the active space is contiguous, or run "
            "method=casci, which supports arbitrary selections."
        )
    return ncore, nact, tuple(plan.nelec), plan


def check_ao_eri_budget(nbf: int, max_memory, section: str,
                        live_tensors: int = 3) -> None:
    """Guard the dense ERI allocations (nbf**4 doubles each) before building them.

    ``live_tensors`` is how many nbf**4 tensors are simultaneously resident at
    the peak, not how many are created in total.  Budgeting a single tensor let
    a job pass and then be OOM-killed anyway: every caller here follows the AO
    build with ``_transform_integrals``, which holds the AO tensor while
    allocating the MO tensor, and the native engine allocates its own nbf**4
    work buffer on top -- three at once.  (Native CASSCF likewise keeps the AO
    record live while casscf_driver allocates ctx%eri.)
    """
    one = 8 * int(nbf) ** 4
    peak_bytes = one * max(1, int(live_tensors))
    budget_bytes = max(1, int(max_memory)) * 1024 * 1024
    if peak_bytes > budget_bytes:
        raise ValueError(
            f"Dense ERI handling for nbf={nbf} needs ~{peak_bytes / 1024 ** 3:.2f} GiB "
            f"({max(1, int(live_tensors))} x {one / 1024 ** 3:.2f} GiB tensors live at "
            f"once: AO, MO and the transform work buffer), exceeding the {section} "
            f"max_memory budget of {int(max_memory)} MiB. "
            "Dense FCI is only intended for small basis sets."
        )


def apply_active_space(
    h1e: np.ndarray,
    eri: np.ndarray,
    plan: ActiveSpacePlan,
    ecore: float,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int], float]:
    """Gather the active ERI block and fold the frozen core into ``h``."""
    active = np.array(plan.active, dtype=int)
    eri_active = eri[np.ix_(active, active, active, active)].copy()
    h_active, ecore_active = _fold_core(
        h1e, eri, list(plan.active), list(plan.core), ecore=float(ecore))
    return h_active, eri_active, plan.nelec, ecore_active


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
    ecore = _float_setting(ecore, "ecore")
    plan = active_space_plan(norb, nelec, settings)
    h_active, eri_active, active_nelec, ecore_active = apply_active_space(
        h1e, eri, plan, ecore)
    return h_active, eri_active, active_nelec, ecore_active, dict(plan.metadata)


# ------------------------------------------------------- native CI entry point
# Mirror of the option schema in source/modules/fci_driver.F90 (FCI_I_* /
# FCI_D_* parameter block).  The Fortran file is authoritative; the names below
# are matched against it by tests/test_fci_solve.py, so the two cannot drift.
_FCI_IOPT = (
    "norb", "nact", "ncore", "nalpha", "nbeta", "nroot", "solver", "maxiter",
    "subspace", "mult", "maxmemory", "nthreads", "want_s2", "guess",
    # Correlated-state irrep selection.  0 = any, which is what a caller that
    # never heard of this leaves in the slot, so the native root selection is
    # byte-identical without it.  When it is non-zero the per-irrep and
    # per-active-orbital XOR codes ride in the tail of the same array.
    "irrep", "nirrep",
)
_FCI_DOPT = ("ecore", "eig_tol", "cutoff", "min_purity")
#: The lengths the RELEASED ``fci_solve`` reads, mirroring FCI_NIOPT_V1 /
#: FCI_NDOPT_V1 in fci_driver.F90 and include/oqp.h.  Everything added after
#: v1.3.x is reached through ``fci_solve_ex``, which is told how much the
#: caller allocated; a binary built against the v1.3.1 header allocates only
#: this much, so the old symbol must never read past it.
_FCI_NIOPT_V1 = 14
_FCI_NDOPT_V1 = 3
_FCI_IOPT_INDEX = {name: i for i, name in enumerate(_FCI_IOPT)}
_FCI_DOPT_INDEX = {name: i for i, name in enumerate(_FCI_DOPT)}
_FCI_SOLVER_CODE = {"auto": 0, "dense": 1, "davidson": 2}
# Determinants are packed into one 64-bit key, so the active space caps out at
# 31 spatial orbitals (62 spin orbitals).
_FCI_MAX_NSPIN = 62


class IrrepSelectionUnavailable(RuntimeError):
    """An irrep was requested but the symmetry data to honour it is missing."""


def resolve_irrep_selection(mol, plan, target_irrep, min_purity=0.5):
    """Resolve ``[ci] irrep`` into the tables the native selector needs.

    Returns ``None`` when no irrep was requested. Otherwise returns
    ``(index, xor_codes, active_orbital_codes, min_purity)`` where the codes are
    the bitmask encoding ``Molecule.stage_mo_irreps`` writes, so the direct
    product of orbital irreps is a bitwise XOR.

    Raises rather than falling back: silently ignoring a symmetry request would
    return roots of the wrong symmetry with no diagnostic, which is the failure
    this is meant to prevent.
    """
    name = str(target_irrep or "any").strip()
    if not name or name.lower() == "any":
        return None

    meta = (getattr(mol, "symmetry_metadata", None) or {}).get("mo_irreps") or {}
    if meta.get("status") != "active":
        raise IrrepSelectionUnavailable(
            f"[ci] irrep={name} needs MO irrep labels, which are not available "
            f"(mo_irreps status: {meta.get('status', 'not staged')}). Enable "
            "[symmetry] and use a geometry whose point group is detected."
        )

    irreps = list(meta.get("irreps") or ())
    lowered = {str(n).strip().lower(): i + 1 for i, n in enumerate(irreps)}
    index = lowered.get(name.lower())
    if index is None:
        raise ValueError(
            f"[ci] irrep={name} is not an irrep of the detected point group; "
            f"available: {', '.join(irreps)}"
        )

    xor_codes = [int(c) for c in (meta.get("xor_codes") or ())]
    if len(xor_codes) != len(irreps):
        raise IrrepSelectionUnavailable(
            "staged irrep table is inconsistent with its XOR codes")

    try:
        mo_index = np.asarray(mol.data["OQP::sym_mo_irrep_a"], dtype=int).ravel()
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise IrrepSelectionUnavailable(
            "OQP::sym_mo_irrep_a is not staged; cannot classify determinants"
        ) from exc

    codes = [0] + xor_codes                     # slot 0 = unclassified
    orb_codes = []
    for p in plan.active:
        p = int(p)
        if p < 0 or p >= mo_index.size:
            raise IrrepSelectionUnavailable(
                f"active orbital {p} is outside the staged MO irrep table")
        idx = int(mo_index[p])
        if idx <= 0:
            # An unclassified orbital would be given code 0, which is also the
            # totally symmetric code -- every determinant containing it would
            # be mislabelled rather than declined.
            raise IrrepSelectionUnavailable(
                f"active orbital {p} has no irrep label (the SCF returned a "
                "mixed or near-degenerate orbital), so determinant irreps "
                "cannot be assigned; run without [ci] irrep")
        orb_codes.append(codes[idx])

    return (index, tuple(xor_codes), tuple(orb_codes), float(min_purity))


def _fci_solve_backend():
    """liboqp ``(lib, ffi)`` for the one-call CI driver, or ``None``."""
    backend = _lib_backend()
    if backend is None:
        return None
    lib, ffi = backend
    # fci_solve_ex is what the driver calls: it is the length-negotiated entry
    # point, and the only one that reports root provenance.  A library old
    # enough to lack it is treated as no backend at all, which is the existing
    # "engine declined" path into the Python driver.
    if not hasattr(lib, "fci_solve_ex"):
        return None
    return lib, ffi


def _lib_fci_solve(h1e, eri, plan, spec, *, nthreads, want_s2, use_target_spin):
    """A complete CI solve inside liboqp, or ``None`` to use the Python driver.

    One boundary crossing replaces the whole orchestration loop: the frozen-core
    fold, the active ERI gather, the spin-orbital expansion, determinant
    generation, the dense/Davidson dispatch, every Davidson application and
    subspace diagonalization, the spin diagnostics and the root selection all
    happen without returning here.

    Returns ``(energies, civecs, s2 | None)``.  ``None`` means the engine
    declined -- a missing symbol, an active space wider than a 64-bit
    determinant key, an allocation failure, a non-converged Davidson, or a spin
    filter that ran out of roots -- and the caller re-runs :func:`solve_fci`,
    which raises the established message where there is one.
    """
    backend = _fci_solve_backend()
    if backend is None:
        return None
    lib, ffi = backend
    nact = plan.nact
    if nact <= 0 or 2 * nact > _FCI_MAX_NSPIN:
        return None

    # An irrep request rides in plan.metadata: it is resolved once at the
    # driver, where the Molecule (and therefore the staged symmetry tables) is
    # in scope, rather than threaded through every solver signature.
    irrep = (plan.metadata or {}).get("irrep")
    tail = 0
    if irrep is not None:
        _idx, _xor, _orb, _pure = irrep
        tail = len(_xor) + len(_orb)

    iopt = np.zeros(len(_FCI_IOPT) + tail, dtype=np.int32)
    iopt[_FCI_IOPT_INDEX["norb"]] = plan.norb
    iopt[_FCI_IOPT_INDEX["nact"]] = nact
    iopt[_FCI_IOPT_INDEX["ncore"]] = plan.ncore
    iopt[_FCI_IOPT_INDEX["nalpha"]] = spec.nalpha
    iopt[_FCI_IOPT_INDEX["nbeta"]] = spec.nbeta
    iopt[_FCI_IOPT_INDEX["nroot"]] = spec.nroot
    iopt[_FCI_IOPT_INDEX["solver"]] = _FCI_SOLVER_CODE[spec.solver]
    iopt[_FCI_IOPT_INDEX["maxiter"]] = spec.davidson_maxiter
    iopt[_FCI_IOPT_INDEX["subspace"]] = spec.davidson_subspace
    iopt[_FCI_IOPT_INDEX["mult"]] = (
        spec.target_multiplicity or 0) if use_target_spin else 0
    iopt[_FCI_IOPT_INDEX["maxmemory"]] = spec.max_memory
    iopt[_FCI_IOPT_INDEX["nthreads"]] = nthreads
    iopt[_FCI_IOPT_INDEX["want_s2"]] = 1 if want_s2 else 0

    if irrep is not None:
        iopt[_FCI_IOPT_INDEX["irrep"]] = _idx
        iopt[_FCI_IOPT_INDEX["nirrep"]] = len(_xor)
        base = len(_FCI_IOPT)
        iopt[base:base + len(_xor)] = np.asarray(_xor, dtype=np.int32)
        iopt[base + len(_xor):] = np.asarray(_orb, dtype=np.int32)

    dopt = np.zeros(len(_FCI_DOPT), dtype=np.float64)
    dopt[_FCI_DOPT_INDEX["ecore"]] = spec.ecore
    dopt[_FCI_DOPT_INDEX["eig_tol"]] = spec.eig_tol
    dopt[_FCI_DOPT_INDEX["cutoff"]] = spec.integral_cutoff
    if irrep is not None:
        dopt[_FCI_DOPT_INDEX["min_purity"]] = _pure

    active = np.ascontiguousarray(plan.active, dtype=np.int32)
    # cffi will not cast an empty buffer; the engine reads none of it anyway
    core = (np.ascontiguousarray(plan.core, dtype=np.int32) if plan.core
            else np.zeros(1, dtype=np.int32))
    hh = _as_f64c(h1e)
    vv = _as_f64c(eri)
    ndet = comb(nact, spec.nalpha) * comb(nact, spec.nbeta)
    energies = np.zeros(spec.nroot, dtype=np.float64)
    civecs = np.zeros((ndet, spec.nroot), dtype=np.float64)
    s2 = np.zeros(spec.nroot, dtype=np.float64)

    roots = np.zeros(spec.nroot, dtype=np.int32)

    status = int(lib.fci_solve_ex(
        np.int32(iopt.size),
        np.int32(dopt.size),
        ffi.cast("int32_t *", iopt.ctypes.data),
        ffi.cast("double *", dopt.ctypes.data),
        ffi.cast("int32_t *", active.ctypes.data),
        ffi.cast("int32_t *", core.ctypes.data),
        ffi.cast("double *", hh.ctypes.data),
        ffi.cast("double *", vv.ctypes.data),
        ffi.cast("double *", energies.ctypes.data),
        ffi.cast("double *", civecs.ctypes.data),
        ffi.cast("double *", s2.ctypes.data),
        ffi.cast("int32_t *", roots.ctypes.data)))
    if status < 0:
        return None
    # No canonicalize_ci_phase() here: fci_driver.F90 applies the convention on
    # its own copy-out, and repeating it in Python would hide a regression on
    # the Fortran side from the very test that checks for one.
    return (energies, civecs, (s2 if want_s2 else None),
            np.ascontiguousarray(roots, dtype=np.int64))


def fix_ci_phase(ci_vectors):
    """Give each CI vector a deterministic sign.

    Nothing in a CI solve fixes the sign of an eigenvector, so every signed
    quantity built from one -- the off-diagonal elements of a multistate
    effective Hamiltonian, non-adiabatic couplings, spin-orbit matrix elements
    -- carries a sign that is not a property of the calculation.  It changes
    with the BLAS, the thread count and the machine, which is why those
    regression keys are marked ``phase_invariant``: comparing magnitudes hides
    the symptom rather than removing it, and it also discards sign information
    that IS physical (an individual ``H_IJ`` carries the arbitrary product
    ``s_I s_J``, but a loop product ``H_IJ H_JK H_KI`` carries ``+1`` and is
    meaningful).

    The rule is the largest-magnitude coefficient positive, with near-ties
    resolved by the LOWEST index: an exact tie is measure-zero, but two
    coefficients agreeing to round-off are not, and picking by index keeps the
    choice stable on platforms where ``argmax`` alone would not be.
    """
    coeffs = _real_array(ci_vectors, "CI vectors")
    single = coeffs.ndim == 1
    if single:
        coeffs = coeffs[:, None]
    out = np.array(coeffs, dtype=np.float64, copy=True)
    for k in range(out.shape[1]):
        col = out[:, k]
        amax = float(np.abs(col).max()) if col.size else 0.0
        if amax <= 0.0:
            continue
        pin = int(np.argmax(np.abs(col) >= amax * (1.0 - 1.0e-10)))
        if col[pin] < 0.0:
            out[:, k] = -col
    return out[:, 0] if single else _as_f64c(out)


def solve_active_ci(
    h1e: np.ndarray,
    eri: np.ndarray,
    plan: ActiveSpacePlan,
    ecore: float,
    settings: FCISettings,
    *,
    nroot: int,
    want_s2: bool = False,
    use_target_spin: bool = False,
    integral_cutoff: float = 0.0,
    active_section: str = "[fci]",
    ci_section: str = "[fci]",
    want_roots: bool = False,
):
    """One complete CI solve from the FULL MO integrals.

    This is the seam the native driver replaces.  Options are validated here --
    in Python, once, with the established messages -- and the compute path then
    runs entirely inside liboqp.  :func:`solve_fci` stays as the numerical pin
    and the fallback, exactly like ``_gfock_backend()`` in ``casscf.py``.

    ``use_target_spin`` hands the spin filter to the engine (which has to
    re-solve to widen its root window, so bouncing back here for it would defeat
    the point).  ``FCI.energy`` leaves it off because it reports the root
    indices *within the unfiltered window*, which is a reporting concern rather
    than a compute one.

    Returns ``(energies, civecs, s2 | None)``, or, with ``want_roots``, a
    fourth element: each returned root's 0-based index among the roots the
    solve computed, or ``None`` when this path cannot know it.  An irrep filter
    runs inside the engine, so without that the caller would have to invent
    ``0, 1, 2, ...`` -- and publish a B1 root that is really state 1 as root 0.
    """
    spec = resolve_ci_solve(
        plan.nact,
        plan.nelec,
        ecore=ecore,
        nroot=nroot,
        max_det=settings.max_det,
        max_memory=settings.max_memory,
        eig_tol=settings.eig_tol,
        integral_cutoff=integral_cutoff,
        solver=settings.solver,
        davidson_maxiter=settings.davidson_maxiter,
        davidson_subspace=(
            max(settings.davidson_subspace, int(nroot))
            if settings.davidson_subspace else 0
        ),
        target_spin=settings.target_spin if use_target_spin else "any",
        active_section=active_section,
        ci_section=ci_section,
    )

    native = _lib_fci_solve(
        h1e, eri, plan, spec,
        nthreads=_fci_lib_threads(), want_s2=want_s2,
        use_target_spin=use_target_spin,
    )
    if native is not None:
        energies, coeffs, s2, roots = native
        active_nelec = (spec.nalpha, spec.nbeta)
    else:
        # The Python driver has no irrep classification, so honouring an irrep
        # request here is impossible.  Refuse rather than silently return the
        # lowest roots of any symmetry: the two paths disagreeing about which
        # root a run returns is worse than not offering the feature on the
        # fallback.
        if (plan.metadata or {}).get("irrep") is not None:
            raise IrrepSelectionUnavailable(
                "[ci] irrep selection requires the native CI engine, which "
                "declined this solve (missing symbol, active space wider than a "
                "64-bit determinant key, allocation failure, or no root of that "
                "irrep). Re-run without [ci] irrep to use the Python driver."
            )

        h_act, eri_act, active_nelec, ecore_act = apply_active_space(
            h1e, eri, plan, ecore)
        energies, coeffs = solve_fci(
            h_act, eri_act, active_nelec,
            ecore=ecore_act,
            nroot=spec.nroot,
            max_det=spec.max_det,
            max_memory=spec.max_memory,
            eig_tol=spec.eig_tol,
            integral_cutoff=spec.integral_cutoff,
            solver=settings.solver,
            davidson_maxiter=spec.davidson_maxiter,
            davidson_subspace=spec.davidson_subspace,
            target_spin=spec.target_spin,
            active_section=active_section,
            ci_section=ci_section,
        )
        s2 = None
        # The Python driver returns the lowest roots of the window it solved,
        # so the identity map is the truth here -- EXCEPT when it did its own
        # spin filtering, where it drops the indices it selected on.  Say so
        # rather than hand back a plausible-looking lie.
        roots = (None if use_target_spin
                 else np.arange(len(energies), dtype=np.int64))

    # Both engines converge here, so this is the one place a spin-mixed
    # degenerate subspace can be resolved once and have the native and Python
    # paths agree by construction rather than by keeping two implementations in
    # step.  The native civecs are in the same determinant order _determinants
    # produces (pinned by tests/test_fci_solve.py), so one purification serves
    # both.  Skipped entirely when no two roots are degenerate, which is the
    # common case and costs one pass over the energies.
    dets_active = None
    if any(len(c) > 1 for c in degenerate_clusters(energies)):
        dets_active = _determinants(plan.nact, active_nelec)
        coeffs, purified_s2, _mult, changed = spin_purify_degenerate_clusters(
            energies, coeffs, dets_active, plan.nact, active_nelec)
        if changed:
            s2 = purified_s2 if want_s2 else None
    if want_s2 and s2 is None:
        if dets_active is None:
            dets_active = _determinants(plan.nact, active_nelec)
        s2, _multiplicity = fci_spin_diagnostics(
            coeffs, dets_active, plan.nact, active_nelec)
    # Both engines converge here, so this is also the one place a deterministic
    # sign can be imposed once and have the native and Python paths agree.
    coeffs = fix_ci_phase(coeffs)
    if not want_roots:
        return energies, coeffs, (s2 if want_s2 else None)
    return energies, coeffs, (s2 if want_s2 else None), roots


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

        h1e_mo, eri_mo, plan, ecore, metadata = self._native_mo_integrals()
        nact = plan.nact
        nelec = plan.nelec
        target_multiplicity = _target_spin_multiplicity(self.settings.target_spin)
        determinant_count = int(metadata["determinants"])
        # Enforce the cap BEFORE enumerating.  _determinants materializes the
        # full alpha x beta product as Python tuples, so a half-filled
        # 20-orbital active space tried to build tens of billions of objects
        # and was OOM-killed on the way to the ValueError that resolve_ci_solve
        # would have raised.  metadata["determinants"] is the same binomial
        # count that cap is compared against, and it is already available here.
        if determinant_count > int(self.settings.max_det):
            raise ValueError(
                f"{self.data_prefix} determinant space has {determinant_count} "
                f"determinants, exceeding max_det={self.settings.max_det}. "
                f"Reduce the active space ({self.active_section} "
                f"frozen_core/active_orbitals) or raise "
                f"{self.active_section} max_det."
            )
        determinants = _determinants(nact, nelec)

        def solve_and_diagnose(window_nroot: int):
            # The spin filter stays here rather than in the engine: the root
            # indices reported below are positions in the UNFILTERED window, so
            # this is a reporting concern, not part of the compute path.
            window_energies, window_coeffs, window_s2, window_roots = solve_active_ci(
                h1e_mo,
                eri_mo,
                plan,
                ecore,
                self.settings,
                nroot=window_nroot,
                want_s2=True,
                use_target_spin=False,
                integral_cutoff=self.settings.integral_cutoff,
                active_section=self.active_section,
                ci_section=self.ci_section,
                want_roots=True,
            )
            window_multiplicity = np.maximum(
                np.rint(np.sqrt(np.maximum(0.0, 1.0 + 4.0 * window_s2))).astype(np.int64),
                1,
            )
            if window_roots is None:
                window_roots = np.arange(len(window_energies), dtype=np.int64)
            return (window_energies, window_coeffs, _as_f64c(window_s2),
                    np.ascontiguousarray(window_multiplicity, dtype=np.int64),
                    np.ascontiguousarray(window_roots, dtype=np.int64))

        if target_multiplicity is None:
            energies, coeffs, s2, multiplicity, root_indices = solve_and_diagnose(
                self.settings.nroot)
            # target_spin=any skips the filter, so nothing else would look at
            # these labels -- but they are still reported, and a label that is
            # not a spin eigenvalue must not be reported silently.
            warn_unreliable_spin_labels(s2, multiplicity, nelec,
                                        ci_label=self.data_prefix)
        else:
            # Checked before the retry loop, not inside it: widening the window
            # cannot conjure a multiplicity the electron count forbids, and the
            # loop would grow to the full determinant space to find that out.
            impossible = _impossible_multiplicity(target_multiplicity, nelec)
            if impossible:
                raise ValueError(
                    f"{self.data_prefix} {self.ci_section} "
                    f"target_spin={self.settings.target_spin} is impossible "
                    f"here: {impossible}."
                )
            solve_nroot = min(determinant_count, max(1, int(self.settings.nroot)))
            while True:
                energies, coeffs, s2, multiplicity, window_roots = solve_and_diagnose(
                    solve_nroot)
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
                        nelec=nelec,
                    )
                    # _filter_roots_by_target_spin reports positions inside the
                    # window; map them back onto the solve's own root numbering,
                    # which an irrep filter may already have made non-contiguous.
                    root_indices = np.ascontiguousarray(
                        np.asarray(window_roots, dtype=np.int64)[root_indices],
                        dtype=np.int64)
                    break
                except SpinLabelAmbiguityError:
                    # A mixed root strictly inside the window means its whole
                    # degenerate manifold was solved and purification still
                    # failed -- more roots cannot help.  A mixed root AT the
                    # window's edge means the window may have cut the manifold,
                    # so purification saw a truncated cluster; widening closes
                    # it.
                    problems = spin_label_diagnosis(s2, multiplicity, nelec)
                    at_edge = any(
                        abs(float(energies[root]) - float(energies[-1]))
                        <= DEGENERACY_TOLERANCE
                        for root, _s2v, _m, _why in problems
                    )
                    if not at_edge or solve_nroot >= determinant_count:
                        raise
                    solve_nroot = min(determinant_count, max(solve_nroot + 1, 2 * solve_nroot))
                except ValueError:
                    if solve_nroot >= determinant_count:
                        raise
                    solve_nroot = min(determinant_count, max(solve_nroot + 1, 2 * solve_nroot))

        self.mol.energies = energies.tolist()
        self.mol.data[f"OQP::{self.data_prefix}_ENERGIES"] = _as_f64c(energies)
        self.mol.data[f"OQP::{self.data_prefix}_CI_VECTORS"] = _as_f64c(coeffs)
        self.mol.data[f"OQP::{self.data_prefix}_S2"] = _as_f64c(s2)
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
                "energies": _as_f64c(energies),
                "ci_vectors": _as_f64c(coeffs),
                "s2": _as_f64c(s2),
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
                nact,
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
                nact,
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
        # Same budget the CASCI path now applies, on the separate FCI code path:
        # nroot copies of an norb^4 RDM2, with the comprehension's list still
        # live while np.stack allocates the copy, so the peak is twice the
        # stack.  FCI CAS(2,30) at nroot=100 clears the ~105 MiB CI allocation
        # checks and then asks ~618 MiB, over 1.2 GiB at the peak.
        _nr = int(coeffs.shape[1])
        _rdm_bytes = 8 * _nr * (int(norb) ** 4 + int(norb) ** 2)
        _budget = max(1, int(self.settings.max_memory)) * 1024 ** 2
        if 2 * _rdm_bytes > _budget:
            raise ValueError(
                "%s per-root RDM records for %d root(s) at norb=%d need "
                "~%.2f GiB (~%.2f GiB at the stacking peak), exceeding the "
                "%s max_memory budget of %d MiB.  Reduce the root count or the "
                "active space, or raise %s max_memory."
                % (self.data_prefix, _nr, int(norb), _rdm_bytes / 1024 ** 3,
                   2 * _rdm_bytes / 1024 ** 3, self.active_section,
                   int(self.settings.max_memory), self.active_section))
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
            rdm1=_as_f64c(rdm1_roots[0]),
            rdm2=_as_f64c(rdm2_roots[0]),
            rdm1_roots=_as_f64c(rdm1_roots),
            rdm2_roots=_as_f64c(rdm2_roots),
            natural_occupations=np.ascontiguousarray(
                natural_occupation_roots[0],
                dtype=np.float64,
            ),
            natural_occupation_roots=np.ascontiguousarray(
                natural_occupation_roots,
                dtype=np.float64,
            ),
            energies=_as_f64c(energies),
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
        check_ao_eri_budget(nbf, self.settings.max_memory, self.active_section)

    def _stage_irrep_selection(self, plan):
        """Resolve an ``[ci] irrep`` request onto ``plan.metadata``.

        Done here, in the driver, because this is where the Molecule -- and so
        the tables ``stage_mo_irreps`` wrote -- is in scope.  Riding in
        ``plan.metadata`` keeps it out of every solver signature between here
        and ``_lib_fci_solve``.

        Inherited by CASCI, which builds its own plan from its own orbitals.
        """
        selection = resolve_irrep_selection(
            self.mol, plan,
            getattr(self.settings, "target_irrep", "any"),
            getattr(self.settings, "irrep_min_purity", 0.5),
        )
        if selection is not None:
            plan.metadata["irrep"] = selection
        return selection

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
        # The FULL MO integrals are returned: the active gather and the
        # frozen-core fold now happen inside the native CI driver, so doing
        # them here would only be work thrown away.
        plan = active_space_plan(h1e.shape[0], nelec, self.settings)
        self._stage_irrep_selection(plan)
        metadata = dict(plan.metadata)
        metadata["orbital_source"] = "rhf"
        return h1e, eri, plan, ecore, metadata
