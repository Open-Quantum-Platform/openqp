from pathlib import Path
import importlib.util

import numpy as np


_NAC_UTILS_PATH = Path(__file__).resolve().parents[1] / "pyoqp" / "oqp" / "library" / "nac_utils.py"
_spec = importlib.util.spec_from_file_location("nac_utils_under_test", _NAC_UTILS_PATH)
assert _spec is not None and _spec.loader is not None
nac_utils = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nac_utils)

align_state_overlap = nac_utils.align_state_overlap
closest_orthogonal = nac_utils.closest_orthogonal
stable_nac_from_overlap = nac_utils.stable_nac_from_overlap
symmetrize_derivative_coupling = nac_utils.symmetrize_derivative_coupling
verify_nac_conventions = nac_utils.verify_nac_conventions


def rotation_from_generator(d, dt):
    """Small matrix exponential for a real antisymmetric generator."""
    vals, vecs = np.linalg.eig(d * dt)
    return np.real_if_close(vecs @ np.diag(np.exp(vals)) @ np.linalg.inv(vecs)).real


def test_phase_and_permutation_alignment_recovers_canonical_overlap():
    dt = 0.05
    d = np.array(
        [
            [0.0, 0.12, -0.03],
            [-0.12, 0.0, 0.08],
            [0.03, -0.08, 0.0],
        ]
    )
    canonical = rotation_from_generator(d, dt)
    permutation = np.array([2, 0, 1])
    signs = np.array([-1.0, 1.0, -1.0])

    raw = canonical[:, permutation] * signs.reshape((1, -1))
    aligned, diagnostics = align_state_overlap(raw)

    assert diagnostics["permutation"].tolist() == [1, 2, 0]
    assert np.allclose(aligned, canonical, atol=1.0e-12)
    assert not diagnostics["ambiguous"]


def test_stable_nac_is_invariant_to_random_state_signs():
    dt = 0.01
    d = np.array(
        [
            [0.0, 0.20, -0.04],
            [-0.20, 0.0, 0.11],
            [0.04, -0.11, 0.0],
        ]
    )
    overlap = rotation_from_generator(d, dt)
    flipped = overlap * np.array([-1.0, 1.0, -1.0]).reshape((1, -1))

    dc_ref, _, _, ref_diag = stable_nac_from_overlap(overlap, dt, method="polar")
    dc_flip, _, _, flip_diag = stable_nac_from_overlap(flipped, dt, method="polar")

    assert np.allclose(dc_flip, dc_ref, atol=1.0e-12)
    assert ref_diag["antisymmetry_error"] < 1.0e-14
    assert flip_diag["antisymmetry_error"] < 1.0e-14


def test_polar_projection_suppresses_symmetric_overlap_noise():
    dt = 0.02
    d = np.array(
        [
            [0.0, 0.30, 0.02],
            [-0.30, 0.0, -0.07],
            [-0.02, 0.07, 0.0],
        ]
    )
    clean = rotation_from_generator(d, dt)
    symmetric_noise = np.array(
        [
            [0.0, 0.04, -0.02],
            [0.04, 0.0, 0.03],
            [-0.02, 0.03, 0.0],
        ]
    )
    noisy = clean + symmetric_noise

    dc_noisy, _, _, diag = stable_nac_from_overlap(noisy, dt, method="polar")
    dc_clean, _, _, _ = stable_nac_from_overlap(clean, dt, method="polar")

    raw_error = np.linalg.norm(noisy.T @ noisy - np.eye(3))
    used_error = diag["orthogonality_error_used"]
    assert used_error < raw_error
    assert np.allclose(closest_orthogonal(noisy).T @ closest_orthogonal(noisy), np.eye(3), atol=1.0e-12)
    assert np.allclose(dc_noisy + dc_noisy.T, 0.0, atol=1.0e-14)


def test_numerical_nac_vector_cleanup_projects_state_pair_noise():
    dcv = np.zeros((3, 3, 1, 2))
    dcv[0, 1, 0, 0] = 1.0
    dcv[1, 0, 0, 0] = -0.8  # symmetric noise contaminates the pair
    dcv[0, 2, 0, 1] = -0.4
    dcv[2, 0, 0, 1] = 0.6

    cleaned = symmetrize_derivative_coupling(dcv)

    assert np.allclose(cleaned + np.swapaxes(cleaned, 0, 1), 0.0)
    assert np.allclose(cleaned[0, 0], 0.0)
    assert cleaned[0, 1, 0, 0] == -cleaned[1, 0, 0, 0]


def test_verify_nac_conventions_can_warn_without_asserting():
    dc = np.array([[0.0, 1.0], [-0.9, 0.0]])
    nac = np.array([[0.0, 2.0], [2.2, 0.0]])

    d_res, h_res = verify_nac_conventions(dc, nac, atol=1.0e-12, strict=False)

    assert d_res > 0.0
    assert h_res > 0.0
