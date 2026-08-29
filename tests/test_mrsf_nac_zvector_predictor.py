"""Numerical and structural gates for the resident NAC Z-vector predictor."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INTERCHANGE = (ROOT / "source/modules/mrsf_nac_interchange.F90").read_text()
CPHF = (ROOT / "source/modules/cphf.F90").read_text()
NAMD = (ROOT / "pyoqp/oqp/library/namd.py").read_text()


def _orthogonal(rng, n):
    q, r = np.linalg.qr(rng.normal(size=(n, n)))
    q *= np.sign(np.diag(r))[None, :]
    return q


def test_block_procrustes_transport_reproduces_a_known_orbital_rotation():
    rng = np.random.default_rng(20260830)
    old = rng.normal(size=(3, 2))
    old_to_current_left = _orthogonal(rng, 3).T
    old_to_current_right = _orthogonal(rng, 2).T

    # get_basis_overlap stores <old|current>.  The nearest orthogonal factor
    # of its transpose is therefore the old-to-current coordinate map.
    overlap_left = old_to_current_left.T
    overlap_right = old_to_current_right.T
    ul, _, vtl = np.linalg.svd(overlap_left.T)
    ur, _, vtr = np.linalg.svd(overlap_right.T)
    ql = ul @ vtl
    qr = ur @ vtr
    transported = ql @ old @ qr.T
    oracle = old_to_current_left @ old @ old_to_current_right.T
    np.testing.assert_allclose(transported, oracle, atol=2.0e-15)


def test_linear_predictor_is_transport_of_same_frame_extrapolation():
    rng = np.random.default_rng(4896372)
    recent = rng.normal(size=(4, 3))
    earlier = rng.normal(size=(4, 3))
    ql = _orthogonal(rng, 4)
    qr = _orthogonal(rng, 3)
    eta = 0.7
    candidate = ql @ (recent + eta * (recent - earlier)) @ qr.T
    transported = (
        ql @ recent @ qr.T
        + eta * (ql @ recent @ qr.T - ql @ earlier @ qr.T)
    )
    np.testing.assert_allclose(candidate, transported, atol=2.0e-15)


def test_predictor_never_weakens_the_certified_minres_solution():
    assert "residual_sq < sum(bvec(:,irhs)**2)" in CPHF
    assert "initial_guess_accepted(irhs) = .true." in CPHF
    assert "tol=1.0e-20_dp" in INTERCHANGE
    assert "scaled_residual_sq > fallback_rel_sq" in INTERCHANGE
    assert "OQP_MRSF_NAC_ZV_MAX_DISP" in INTERCHANGE
    assert "OQP_MRSF_NAC_ZV_OVERLAP_MIN" in INTERCHANGE
    assert "recent_current +" in INTERCHANGE
    assert "eta*(recent_current-earlier_current)" in INTERCHANGE


def test_approximate_mode_replaces_only_z_and_has_exact_refreshes():
    assert "trim(mode) == 'linear_approx'" in INTERCHANGE
    assert "solution = guess" in INTERCHANGE
    assert "residual = -1.0_dp" in INTERCHANGE
    assert "iterations = 0" in INTERCHANGE
    assert "OQP_MRSF_NAC_ZV_EXACT_EVERY" in INTERCHANGE
    assert "call cphf_solve_rohf" in INTERCHANGE
    assert INTERCHANGE.index("if (use_approximation) then") < INTERCHANGE.index(
        "call cphf_solve_rohf"
    )


def test_namd_records_full_vector_and_velocity_contraction_predictor_errors():
    for token in (
        ".namd.zpredict.tsv",
        "production_is_predictor",
        "d_relative_l2",
        "d_cosine_min",
        "h_relative_l2",
        "vd_relative_l2",
        "tracking_overlap_min",
    ):
        assert token in NAMD


def test_gradient_and_nac_can_share_one_native_multi_rhs_solve():
    driver = (ROOT / "source/modules/mrsf_nac_driver.F90").read_text()
    zvector = (ROOT / "source/modules/tdhf_mrsf_z_vector.F90").read_text()
    single_point = (ROOT / "pyoqp/oqp/library/single_point.py").read_text()
    assert "fused_rhs(:,1) = 2.0_dp*gradient_rhs" in driver
    assert "solution_batch = fused_solution(:,2:npair+1)" in driver
    assert "rhs_tolerances=fused_tolerance" in driver
    assert "NAC_GRADIENT_Z_FUSION" in driver
    assert "OQP_MRSF_NAC_ZV_FUSE_GRADIENT" in zvector
    assert "mrsf_nac_fusion_set_rhs(rhs,cnvtol)" in zvector
    assert "mrsf_nac_fusion_get_tolerance(gradient_tolerance)" in driver
    assert "fused_tolerance(1) = max(1.0e-20_dp, gradient_tolerance)" in driver
    assert "mrsf_nac_lagrangian_fused_external(infos)" in zvector
    assert "mrsf_nac_fusion_take_solution(xk)" in zvector
    assert "_nac_fused_gradient_ready" in single_point
    assert "self._state_overlap(istep, update_analytic=False)" in NAMD
