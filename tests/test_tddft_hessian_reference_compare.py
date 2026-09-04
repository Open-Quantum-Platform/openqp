import importlib.util
import json
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).parent / "tools" / "compare_tddft_hessian_references.py"
SPEC = importlib.util.spec_from_file_location("compare_tddft_hessian_references", MODULE_PATH)
COMPARATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(COMPARATOR)


def _write_openqp(path, hessian, frequencies):
    path.write_text(
        json.dumps({"atoms": ["H", "H"], "hessian": hessian.tolist(), "freqs": frequencies}),
        encoding="utf-8",
    )


def test_gamess_parser_and_case_metrics(tmp_path):
    hessian = np.arange(36.0).reshape(6, 6)
    hessian = hessian + hessian.T
    dat = tmp_path / "reference.dat"
    lines = [" $HESS", "ENERGY IS 0.0"]
    for row, values in enumerate(hessian, 1):
        lines.append(f"{row:2d}  1" + "".join(f"{value:16.8E}" for value in values[:5]))
        lines.append(f"{row:2d}  2" + "".join(f"{value:16.8E}" for value in values[5:]))
    lines += [" $END", " MODE 1 FREQUENCY= 123.00000 (CM**-1)"]
    dat.write_text("\n".join(lines), encoding="utf-8")

    analytic = tmp_path / "analytic.json"
    numerical = tmp_path / "numerical.json"
    _write_openqp(analytic, hessian, [-123.0])
    perturbed = hessian.copy()
    perturbed[2, 4] += 0.25
    _write_openqp(numerical, perturbed, [-120.0])

    parsed = COMPARATOR.load_gamess_dat(dat)
    np.testing.assert_allclose(parsed["hessian"], hessian)
    assert parsed["frequencies"] == [123.0]

    result = COMPARATOR.compare_case("synthetic", dat, analytic, numerical)
    metrics = result["comparisons"]["analytic_vs_numerical"]
    assert metrics["max_abs"] == 0.25
    assert metrics["max_abs_index_1based"] == [3, 5]
    assert np.isclose(metrics["rms"], 0.25 / 6.0)
    assert result["symmetry"]["openqp_analytic"]["max_abs"] == 0.0
    assert result["symmetry"]["openqp_numerical"]["max_abs"] == 0.25
    assert result["frequencies"]["analytic_vs_gamess"]["max_abs"] == 0.0


def test_translational_residual_is_zero_for_block_laplacian():
    atom_laplacian = np.array([[1.0, -1.0], [-1.0, 1.0]])
    hessian = np.kron(atom_laplacian, np.eye(3))
    residual = COMPARATOR.translational_metrics(hessian, natom=2)
    assert residual["right_max_abs"] == 0.0
    assert residual["left_max_abs"] == 0.0
