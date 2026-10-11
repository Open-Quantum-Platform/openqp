"""Spherical and core-removing ECP MP2 gradients against energy derivatives."""
import configparser
import json
from pathlib import Path

import numpy as np
import pytest
from oqp.pyoqp import Runner

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    "MP2/h2o_rmp2_cc-pvdz_spherical_grad",
    "MP2/h2o_rmp2_cc-pvtz_spherical_grad",
    "ECP/HBr_RMP2_LANL2DZ_GRAD",
)


def _run(tmp_path, case, name, coords=None, gradient=True):
    parser = configparser.ConfigParser()
    parser.read(ROOT / "examples" / (case + ".inp"))
    config = {section: dict(parser[section]) for section in parser.sections()}
    config["input"]["runtype"] = "grad" if gradient else "energy"
    runner = Runner(project=name, input_file=None, log=str(tmp_path / (name + ".log")),
                    input_dict=config, silent=1, usempi=False)
    if coords is not None:
        runner.mol.update_system(coords)
    runner.run(test_mod=True)
    return runner


@pytest.mark.parametrize("case", CASES)
def test_spherical_and_ecp_gradient_matches_central_difference(tmp_path, case):
    runner = _run(tmp_path, case, "analytic")
    x = np.asarray(runner.mol.get_system(), dtype=float).reshape(-1).copy()
    grad = np.asarray(runner.mol.grads[0], dtype=float).reshape(-1)
    reference = json.loads((ROOT / "examples" / (case + ".json")).read_text())
    assert abs(float(runner.mol.energies[0]) - reference["energy"]) < 1e-8
    np.testing.assert_allclose(grad, reference["grad"], atol=1e-8, rtol=0)
    if case.startswith("MP2/"):
        expected_nbf = 58 if "cc-pvtz" in case else 24
        assert int(runner.mol.data.get_basis()["nbf"]) == expected_nbf
    # All components, including symmetry zeros, use the energy derivative.
    step = 1e-3  # Bohr
    for k in range(x.size):
        energies = []
        for sign in (-1, 1):
            displaced = x.copy()
            displaced[k] += sign * step
            r = _run(tmp_path, case, f"fd_{k}_{sign}", displaced, gradient=False)
            energies.append(float(r.mol.energies[0]))
        fd = (energies[1] - energies[0]) / (2 * step)
        assert abs(grad[k] - fd) < 1e-6, (case, k, grad[k], fd)
    np.testing.assert_allclose(grad.reshape(-1, 3).sum(axis=0), 0, atol=1e-10)
    if case.startswith("ECP/"):
        # Br/LANL2DZ removes 28 core electrons: HBr has eight valence electrons.
        assert runner.mol.mol_prop.nocc == 4


def test_gradient_repeated_calls_with_increasing_and_decreasing_basis(tmp_path):
    first = _run(tmp_path, CASES[0], "small_first")
    expected = np.asarray(first.mol.grads[0], dtype=float).copy()
    for k, case in enumerate((CASES[1], CASES[2], CASES[0])):
        runner = _run(tmp_path, case, f"repeat_{k}")
    np.testing.assert_allclose(runner.mol.grads[0], expected, atol=1e-10, rtol=0)
