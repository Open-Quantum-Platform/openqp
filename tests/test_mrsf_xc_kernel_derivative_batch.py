"""Structural gates for exact batched MRSF XC-kernel derivatives."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/dftlib/dft_gridint_mrsf_xc_kernel_derivative.F90"


def test_point_workspace_is_reused_across_grid_points():
    text = SOURCE.read_text()
    start = text.index("  subroutine accumulate_kernel_point")
    stop = text.index("  end subroutine accumulate_kernel_point", start)
    point = text[start:stop]
    assert "workspace" in point
    assert "allocate(" not in point
    assert "deallocate(" not in point


def test_all_density_and_response_fields_share_ao_pair_passes():
    text = SOURCE.read_text()
    start = text.index("  subroutine accumulate_kernel_point")
    stop = text.index("  end subroutine accumulate_kernel_point", start)
    point = text[start:stop]
    assert "call field_value_gradient_batch" in point
    assert "call gga_density_nuclear_point_first_batch" in point
    assert "call response_field_value_gradient_batch" in point
    assert "call total_density_derivative" not in point
