"""Active-space reduced density matrices from determinant CI vectors."""

from __future__ import annotations

from itertools import combinations

import numpy as np


def _bit_count(value: int) -> int:
    return bin(int(value)).count("1")


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
    """Return the spin-orbital 2-RDM ``D[p,q,r,s] = <a_p^+ a_q^+ a_s a_r>``."""
    coeff, dets, n_spinorb = _validate_ci_inputs(ci_vec, determinants, n_spinorb)
    det_index = _determinant_index(dets)
    dtype = np.result_type(coeff.dtype, np.float64)
    rdm2 = np.zeros((n_spinorb, n_spinorb, n_spinorb, n_spinorb), dtype=dtype)

    for col, det in enumerate(dets):
        c_ket = coeff[col]
        if c_ket == 0:
            continue
        for r in _occupied(det, n_spinorb):
            ann_r = annihilate(det, r)
            if ann_r is None:
                continue
            det_r, phase_r = ann_r
            for s in _occupied(det_r, n_spinorb):
                ann_s = annihilate(det_r, s)
                if ann_s is None:
                    continue
                det_rs, phase_s = ann_s
                for q in range(n_spinorb):
                    cre_q = create(det_rs, q)
                    if cre_q is None:
                        continue
                    det_qrs, phase_q = cre_q
                    for p in range(n_spinorb):
                        cre_p = create(det_qrs, p)
                        if cre_p is None:
                            continue
                        det_pqrs, phase_p = cre_p
                        row = det_index.get(det_pqrs)
                        if row is not None:
                            rdm2[p, q, r, s] += (
                                np.conjugate(coeff[row])
                                * c_ket
                                * phase_r
                                * phase_s
                                * phase_q
                                * phase_p
                            )

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
