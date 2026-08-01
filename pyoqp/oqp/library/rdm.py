"""Active-space reduced density matrices from determinant CI vectors."""

from __future__ import annotations

from itertools import combinations

import numpy as np


try:  # Python >= 3.10
    _POPCOUNT = int.bit_count
except AttributeError:  # pragma: no cover - Python 3.9 fallback
    def _POPCOUNT(value: int) -> int:
        return bin(value).count("1")


def _bit_count(value: int) -> int:
    # Hot path: this is called hundreds of millions of times when building the
    # 3- and 4-particle RDMs, so it goes through int.bit_count where available
    # rather than materialising a bin() string and counting characters.
    return _POPCOUNT(int(value))


def annihilate(det: int, orb: int) -> tuple[int, int] | None:
    """Apply ``a_orb`` to a determinant bitstring.

    The bit ordering matches :mod:`oqp.library.fci`: lower bits are alpha
    spin orbitals, followed by beta spin orbitals for spatial active spaces.
    """
    bit = 1 << int(orb)
    det = int(det)
    if not det & bit:
        return None
    phase = -1 if _bit_count(det & (bit - 1)) % 2 else 1
    return det ^ bit, phase


def create(det: int, orb: int) -> tuple[int, int] | None:
    """Apply ``a_orb^+`` to a determinant bitstring."""
    bit = 1 << int(orb)
    det = int(det)
    if det & bit:
        return None
    phase = -1 if _bit_count(det & (bit - 1)) % 2 else 1
    return det | bit, phase


def determinant_basis(norb: int, nelec: int | tuple[int, int] | list[int]) -> list[int]:
    """Return determinant bitstrings in the native FCI/CASCI ordering."""
    nalpha, nbeta = _as_nelec_pair(nelec)
    alpha = _string_basis(norb, nalpha)
    beta = _string_basis(norb, nbeta)
    return [a | (b << norb) for a in alpha for b in beta]


def make_rdm1_spinorb(
    ci_vec: np.ndarray,
    determinants: list[int] | tuple[int, ...],
    n_spinorb: int,
) -> np.ndarray:
    """Return the spin-orbital 1-RDM ``D[p,q] = <a_p^+ a_q>``."""
    coeff, dets, n_spinorb = _validate_ci_inputs(ci_vec, determinants, n_spinorb)
    det_index = _determinant_index(dets)
    dtype = np.result_type(coeff.dtype, np.float64)
    rdm1 = np.zeros((n_spinorb, n_spinorb), dtype=dtype)

    for col, det in enumerate(dets):
        c_ket = coeff[col]
        if c_ket == 0:
            continue
        for q in _occupied(det, n_spinorb):
            ann_q = annihilate(det, q)
            if ann_q is None:
                continue
            det_q, phase_q = ann_q
            for p in range(n_spinorb):
                cre_p = create(det_q, p)
                if cre_p is None:
                    continue
                det_pq, phase_p = cre_p
                row = det_index.get(det_pq)
                if row is not None:
                    rdm1[p, q] += np.conjugate(coeff[row]) * c_ket * phase_q * phase_p

    return _finalize(rdm1)


def make_rdm2_spinorb(
    ci_vec: np.ndarray,
    determinants: list[int] | tuple[int, ...],
    n_spinorb: int,
) -> np.ndarray:
    """Return the spin-orbital 2-RDM ``D[p,q,r,s] = <a_p^+ a_q^+ a_s a_r>``.

    Built by factorising through the doubly-annihilated intermediate space
    instead of walking creations back up onto every determinant.  Writing
    ``X[(a,b), t] = sum_d c_d <t| a_b a_a |d>``, the bra pair (p,q) and the ket
    pair (r,s) of ``<a_p^+ a_q^+ a_s a_r>`` meet on that shared intermediate,
    so the whole 2-RDM is one Gram matrix::

        D[p,q,r,s] = sum_t conj(X[(p,q), t]) * X[(r,s), t]

    That is a single GEMM.  The direct nested form instead called ``create``
    once per (p, q) per intermediate and threw away every hit that landed on an
    occupied orbital -- 370M calls for a CAS(8,8) reference, versus roughly
    ndet * nelec^2 cheap insertions here.
    """
    coeff, dets, n_spinorb = _validate_ci_inputs(ci_vec, determinants, n_spinorb)
    dtype = np.result_type(coeff.dtype, np.float64)
    npair = n_spinorb * n_spinorb

    # X, assembled as {intermediate determinant -> column}; only intermediates
    # that are actually reachable get a column.
    inter_index: dict[int, int] = {}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []

    for col, det in enumerate(dets):
        c_ket = coeff[col]
        if c_ket == 0:
            continue
        occ = _occupied(det, n_spinorb)
        for a in occ:
            ann_a = annihilate(det, a)
            if ann_a is None:
                continue
            det_a, phase_a = ann_a
            for b in _occupied(det_a, n_spinorb):
                ann_b = annihilate(det_a, b)
                if ann_b is None:
                    continue
                det_ab, phase_b = ann_b
                t = inter_index.get(det_ab)
                if t is None:
                    t = len(inter_index)
                    inter_index[det_ab] = t
                rows.append(a * n_spinorb + b)
                cols.append(t)
                vals.append(phase_a * phase_b * c_ket)

    if not vals:
        return _finalize(
            np.zeros((n_spinorb,) * 4, dtype=dtype)
        )

    x = np.zeros((npair, len(inter_index)), dtype=dtype)
    np.add.at(x, (np.asarray(rows), np.asarray(cols)), np.asarray(vals, dtype=dtype))

    rdm2 = (x.conj() @ x.T).reshape(n_spinorb, n_spinorb, n_spinorb, n_spinorb)
    return _finalize(rdm2)


def make_rdm1_spatial(
    ci_vec: np.ndarray,
    determinants: list[int] | tuple[int, ...],
    n_spatial_orb: int,
) -> np.ndarray:
    """Return the spin-summed spatial 1-RDM for alpha-then-beta bitstrings."""
    n_spatial_orb = int(n_spatial_orb)
    spin = make_rdm1_spinorb(ci_vec, determinants, 2 * n_spatial_orb)
    spatial = spin[:n_spatial_orb, :n_spatial_orb] + spin[n_spatial_orb:, n_spatial_orb:]
    return _finalize(spatial)


def make_rdm2_spatial(
    ci_vec: np.ndarray,
    determinants: list[int] | tuple[int, ...],
    n_spatial_orb: int,
) -> np.ndarray:
    """Return spin-summed spatial 2-RDM in chemists' integral order.

    The returned tensor is
    ``D[p,q,r,s] = sum_{sigma,tau} <a_p_sigma^+ a_r_tau^+ a_s_tau a_q_sigma>``.
    It contracts with spatial integrals as
    ``0.5 * einsum("pqrs,pqrs", eri, D)`` for ``eri[p,q,r,s] = (pq|rs)``.
    """
    n_spatial_orb = int(n_spatial_orb)
    spin = make_rdm2_spinorb(ci_vec, determinants, 2 * n_spatial_orb)
    spatial = np.zeros(
        (n_spatial_orb, n_spatial_orb, n_spatial_orb, n_spatial_orb),
        dtype=spin.dtype,
    )
    for sigma_offset in (0, n_spatial_orb):
        for tau_offset in (0, n_spatial_orb):
            block = spin[
                sigma_offset : sigma_offset + n_spatial_orb,
                tau_offset : tau_offset + n_spatial_orb,
                sigma_offset : sigma_offset + n_spatial_orb,
                tau_offset : tau_offset + n_spatial_orb,
            ]
            spatial += np.transpose(block, (0, 2, 1, 3))
    return _finalize(spatial)


def natural_orbital_occupations(rdm1_spatial: np.ndarray) -> np.ndarray:
    """Return spin-summed natural orbital occupations sorted descending."""
    matrix = _real_square_matrix(rdm1_spatial, "spatial 1-RDM")
    if not np.allclose(matrix, matrix.T, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("spatial 1-RDM must be symmetric")
    from oqp.library.fci import _symmetric_eigh

    occupations = _symmetric_eigh(0.5 * (matrix + matrix.T))[0]
    return np.ascontiguousarray(occupations[::-1], dtype=np.float64)


def _as_nelec_pair(nelec: int | tuple[int, int] | list[int]) -> tuple[int, int]:
    if isinstance(nelec, (tuple, list)):
        if len(nelec) != 2:
            raise ValueError("nelec must be an integer or an (nalpha, nbeta) pair")
        return int(nelec[0]), int(nelec[1])
    if int(nelec) % 2:
        raise ValueError("integer nelec implies a closed-shell singlet and must be even")
    return int(nelec) // 2, int(nelec) // 2


def _string_basis(norb: int, nelec: int) -> list[int]:
    return [sum(1 << orb for orb in occ) for occ in combinations(range(int(norb)), int(nelec))]


def _occupied(det: int, norb: int) -> list[int]:
    return [orb for orb in range(norb) if int(det) & (1 << orb)]


def _validate_ci_inputs(
    ci_vec: np.ndarray,
    determinants: list[int] | tuple[int, ...],
    n_spinorb: int,
) -> tuple[np.ndarray, list[int], int]:
    n_spinorb = int(n_spinorb)
    if n_spinorb <= 0:
        raise ValueError("n_spinorb must be positive")

    dets = [int(det) for det in determinants]
    if not dets:
        raise ValueError("determinants must be non-empty")
    if any(det < 0 for det in dets):
        raise ValueError("determinants must be non-negative bitstrings")
    if any(det >> n_spinorb for det in dets):
        raise ValueError("determinants contain occupied orbitals outside n_spinorb")

    coeff = np.asarray(ci_vec)
    if coeff.ndim != 1:
        raise ValueError("ci_vec must be a one-dimensional CI vector")
    if coeff.shape[0] != len(dets):
        raise ValueError("CI vector length must match the determinant count")
    return coeff, dets, n_spinorb


def _determinant_index(determinants: list[int]) -> dict[int, int]:
    det_index = {det: idx for idx, det in enumerate(determinants)}
    if len(det_index) != len(determinants):
        raise ValueError("determinants must not contain duplicates")
    return det_index


def _real_square_matrix(values: np.ndarray, label: str) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must contain finite real values") from exc
    if raw.dtype.kind not in {"i", "u", "f"}:
        raise ValueError(f"{label} must contain finite real values")
    matrix = np.ascontiguousarray(raw, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError(f"{label} must be a non-empty square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{label} must be finite")
    return matrix


def _finalize(array: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(np.real_if_close(array, tol=1000))
