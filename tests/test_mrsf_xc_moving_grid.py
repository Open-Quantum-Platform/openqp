"""Contracts for the relaxed-density XC moving-grid contribution.

The MRSF energy derivative contains the linear XC probe P=T+Z.  On an
atom-centred quadrature its derivative requires both the normalized
partition-weight response and motion of the points with the owner atom.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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
    assert "dat%do_weight_derivative .and. dat%do_ground_state" in source
    assert "if (.not. requested_weight_only) then" in source


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
