import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_sg1_forces_original_radial_grid_definition():
    dft = (ROOT / "source" / "dftlib" / "dft.F90").read_text()
    sg1_block = dft[dft.index('case ("SG1")'):dft.index('case ("SG2", "SG3")')]

    assert "integer, parameter :: SG1_NRAD = 50" in dft
    assert "use dft_radial_grid_types, only: dft_radial_grid_mhl" in dft
    assert "pruned%nrad = SG1_NRAD" in sg1_block
    assert "infos%dft%rad_grid_type = dft_radial_grid_mhl" in sg1_block


def test_sg1_uses_fixed_shell_partitions_for_fifty_point_mhl_grid():
    dft = (ROOT / "source" / "dftlib" / "dft.F90").read_text()
    sg1_block = dft[dft.index('case ("SG1")'):dft.index('case ("SG2", "SG3")')]

    assert "integer, parameter :: sg1shells(5,3)" in dft
    for shell_counts in (
        "17, 4, 4, 9, 16",
        "14, 7, 3, 9, 17",
        "12, 7, 5, 7, 19",
    ):
        assert shell_counts in dft

    assert re.search(
        r"allocate\(.*pruned%nradPerRegion\(pruned%ngrids,\s*ntyps\)",
        sg1_block,
        re.DOTALL,
    )
    assert "pruned%nradPerRegion(:, 1:3) = sg1shells" in sg1_block
