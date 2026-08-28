"""Verification of the finite-grid XC nuclear-gradient contributions.

The MRSF energy derivative contains the linear XC probe P=T+Z.  On an
atom-centred quadrature its derivative requires both the normalized
partition-weight response and motion of the points with the owner atom.
"""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = (ROOT / "tests" / "fixtures" / "mrsf_xc_moving_grid"
           / "low_sym_h2o_gradient_fd.json")


def _source(relative):
    return (ROOT / relative).read_text()


def test_mrsf_gradient_requests_probe_only_moving_grid_response():
    source = _source("source/modules/tdhf_mrsf_gradient.F90")
    body = source.split("subroutine tdhf_mrsf_gradient(infos)", 1)[1]
    body = body.split("end subroutine tdhf_mrsf_gradient", 1)[0]

    assert "grid_correction" in body
    assert "include_ground_state=.false." in body
    assert "include_weight_derivative=.true." in body
    assert "weight_derivative_only=.true." in body
    assert "include_ground_state=.true." in body
    assert "grid_p = 0.0_dp" in body


def test_linear_probe_has_partition_and_owner_motion_terms():
    source = _source("source/dftlib/dft_gridint_tdxc_grad.F90")
    body = source.split("subroutine add_partition_weight_gradient", 1)[1]
    body = body.split("end subroutine add_partition_weight_gradient", 1)[0]

    assert "partfunc%deriv(mu)" in body
    assert "dlog(:,b,owner)" in body
    assert "merge(1,0,b == owner)" in body
    assert "sum(tmpGrad, dim=1)" in source
    assert "xce%currAtom" in source


def test_moving_grid_mode_is_restricted_to_one_linear_probe():
    source = _source("source/dftlib/dft_gridint_tdxc_grad.F90")

    assert "requested_weight_derivative .and. doFxc" in source
    assert "requested_weight_derivative .and. nMtx /= 1" in source
    assert "if (.not. requested_weight_only) then" in source
    assert "xc%exc(i)*xce%xyzw(i,4)" in source


def test_ground_state_xc_gradient_has_partition_and_owner_motion_terms():
    source = _source("source/dftlib/dft_gridint_grad.F90")
    body = source.split("subroutine add_partition_weight_gradient", 1)[1]
    body = body.split("end subroutine add_partition_weight_gradient", 1)[0]

    assert "xc%exc(1:numPts)*xce%xyzw(1:numPts,4)" in source
    assert "partfunc%deriv(mu)" in body
    assert "dlog(:,b,owner)" in body
    assert "sum(tmpGrad, dim=1)" in source
    assert "dedft = dedft + dat%nucGrad(:,:,1)" in source


def test_grid_retains_partition_metadata_for_all_surface_shift_modes():
    grid = _source("source/dftlib/dft_molgrid.F90")
    builder = _source("source/dftlib/dft.F90")

    assert "partFunType" in grid
    assert "hasSurfaceShift" in grid
    assert "surfaceShift" in grid
    assert "molGrid%partFunType = dft_partfun" in builder
    assert builder.count("molGrid%surfaceShift = aij") == 2


def test_partition_function_derivatives_include_coordinate_chain_rule():
    source = _source("source/dftlib/dft_partfunc.F90")

    assert "(1.0_fp+x*x)*frac*frac" in source
    assert source.count("df = -SCALEF") == 4


def test_low_symmetry_h2o_matches_five_point_energy_derivatives():
    data = json.loads(FIXTURE.read_text())
    assert data["finite_difference"] == "five-point central"
    assert data["functional"] == "bhhlyp"
    assert data["basis"] == "6-31g*"

    roks = [case for case in data["cases"] if "ROKS" in case["state"]]
    mrsf = [case for case in data["cases"] if "MRSF" in case["state"]]
    assert len(roks) == 1
    assert len(mrsf) == 3
    assert roks[0]["max_abs_error_before_hartree_per_bohr"] > 1.0e-4
    assert roks[0]["max_abs_error_after_hartree_per_bohr"] < 1.0e-8
    assert all(case["max_abs_error_before_hartree_per_bohr"] > 1.0e-4
               for case in mrsf)
    assert all(case["max_abs_error_after_hartree_per_bohr"] < 1.0e-6
               for case in mrsf)
    spread = max(case["max_abs_error_after_hartree_per_bohr"] for case in mrsf) \
        - min(case["max_abs_error_after_hartree_per_bohr"] for case in mrsf)
    assert spread < 5.0e-10
