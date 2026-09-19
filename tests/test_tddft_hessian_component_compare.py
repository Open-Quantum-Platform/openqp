import importlib.util
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).parent / "tools" / "compare_tddft_hessian_components.py"
SPEC = importlib.util.spec_from_file_location("compare_tddft_hessian_components", MODULE_PATH)
COMPARATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(COMPARATOR)


def _write_fixture(tmp_path: Path) -> tuple[Path, Path, dict[str, np.ndarray]]:
    n, nat = 6, 2
    matrices = {
        "TDDXR": np.arange(36.0).reshape(n, n) / 100.0,
        "TDD1G": np.arange(36.0, 72.0).reshape(n, n) / 10.0,
        "TDD2G": -np.arange(72.0, 108.0).reshape(n, n) / 20.0,
    }
    blocks = [
        np.diag([1.0, 2.0, 3.0]),
        np.arange(9.0).reshape(3, 3) / 10.0,
        np.diag([4.0, 5.0, 6.0]),
    ]
    hxc = np.zeros((n, n))
    hxc[:3, :3], hxc[3:, :3], hxc[3:, 3:] = blocks
    hxc[:3, 3:] = blocks[1].T
    matrices["TDH2XC"] = hxc

    gamess_lines = []
    for tag in ("TDDXR", "TDD2G", "TDD1G"):
        label = {"TDDXR": "DXR", "TDD1G": "D1G", "TDD2G": "D2G"}[tag]
        gamess_lines.append(f" ${tag} NXYZ= {n}")
        for k in range(n):
            gamess_lines.append(f" {label} COORD {k + 1}")
            column = matrices[tag][:, k]
            gamess_lines.extend(" ".join(f"{value: .15E}" for value in column[i : i + 3]) for i in range(0, n, 3))
        gamess_lines.append(" $END")
    gamess_lines.append(" $TDH2XC NAT,NAT2= 2 3")
    for block in blocks:
        gamess_lines.extend(" ".join(f"{value: .15E}" for value in row) for row in block)
    gamess_lines.append(" $END")
    gamess = tmp_path / "components.dat"
    gamess.write_text("\n".join(gamess_lines) + "\n", encoding="utf-8")

    openqp_lines = []
    for tag, label in COMPARATOR._OPENQP_LABEL.items():
        openqp_lines.append(label)
        openqp_lines.extend(" ".join(f"{value: .15E}" for value in row) for row in matrices[tag])
    openqp = tmp_path / "run.log.components"
    openqp.write_text("\n".join(openqp_lines) + "\n", encoding="utf-8")
    return gamess, openqp, matrices


def test_component_parsers_preserve_gamess_orientation_and_blocks(tmp_path):
    gamess, openqp, matrices = _write_fixture(tmp_path)
    for tag in ("TDDXR", "TDD1G", "TDD2G"):
        np.testing.assert_allclose(COMPARATOR.load_gamess_coordinate_component(gamess, tag), matrices[tag])
    np.testing.assert_allclose(COMPARATOR.load_gamess_tdh2xc(gamess), matrices["TDH2XC"])
    np.testing.assert_allclose(COMPARATOR.load_openqp_component(openqp, "rows_one", 6), matrices["TDD1G"])


def test_component_comparison_reports_direct_and_transposed_controls(tmp_path):
    gamess, openqp, _ = _write_fixture(tmp_path)
    result = COMPARATOR.compare_components(gamess, openqp)
    assert result["dimension"] == 6
    assert all(metrics["max_abs"] == 0.0 for metrics in result["components"].values())
    assert result["combined_rows"]["max_abs"] == 0.0
    assert result["combined_rows_and_hxc"]["max_abs"] == 0.0
    assert result["transposed_control"]["TDDXR"]["max_abs"] > 0.0


def test_metrics_expose_scale_and_largest_signed_residual():
    reference = np.eye(2)
    candidate = 2.0 * reference
    candidate[0, 1] = -0.25
    metrics = COMPARATOR.matrix_metrics(reference, candidate)
    assert metrics["least_squares_scale"] == 2.0
    assert metrics["max_abs_index_1based"] == [1, 1]
    assert metrics["candidate_minus_reference_at_max"] == 1.0
    assert metrics["scaled_max_abs"] == 0.25
