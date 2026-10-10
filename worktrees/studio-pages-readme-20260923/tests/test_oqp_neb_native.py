"""Focused pure-Python coverage for native OpenQP NEB behavior."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


NEB = _load(
    "oqp_neb_native_features_under_test",
    "pyoqp/oqp/library/oqp_neb.py",
).NEB
NEB_UTILS = _load(
    "oqp_neb_utils_native_features_under_test",
    "pyoqp/oqp/library/neb_utils.py",
)


def test_kabsch_align_removes_rigid_rotation_and_translation():
    reference = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.2, 1.1, 0.3], [-0.1, 0.2, 1.4]]
    )
    angle = np.deg2rad(67.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]]
    )
    mobile = reference @ rotation.T + np.array([4.2, -2.7, 1.3])

    aligned = np.asarray(NEB_UTILS.kabsch_align(reference, mobile))

    assert np.allclose(aligned, reference, atol=1.0e-12)


def test_aligned_endpoint_interpolation_removes_pure_rigid_motion(tmp_path):
    reactant = tmp_path / "reactant.xyz"
    product = tmp_path / "product.xyz"
    reactant.write_text(
        "3\nreactant\nC 0 0 0\nH 1 0 0\nH 0 1 0\n"
    )
    # Same geometry, rotated by 90 degrees and translated.
    product.write_text(
        "3\nproduct\nC 2 3 1\nH 2 4 1\nH 1 3 1\n"
    )

    images = NEB_UTILS.interpolate_xyz_endpoints(
        reactant, product, nimage=3, align=True
    )

    expected = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert np.allclose(images[-1].coordinates_angstrom, expected, atol=1.0e-12)
    assert np.allclose(images[1].coordinates_angstrom, expected, atol=1.0e-12)


def _one_dimensional_quadratic(coordinates):
    coordinate = float(coordinates[0])
    return coordinate * coordinate, np.array([2.0 * coordinate])


def test_uneven_straight_band_does_not_converge_with_nonzero_spring_force():
    def flat_surface(coordinates):
        return 0.0, np.zeros_like(coordinates)

    images = [np.array([0.0]), np.array([0.2]), np.array([2.0])]
    neb = NEB(images, k_spring=1.0, climbing=False)

    result = neb.run(
        flat_surface,
        fmax_tol=1.0e-6,
        frms_tol=1.0e-6,
        maxiter=1,
        dt=0.1,
        maxmove=1.0,
    )

    assert not result["converged"]
    assert result["fmax"] > 1.0e-6
    assert result["frms"] > 1.0e-6
    assert result["images"][1][0] != pytest.approx(0.2)


def test_rms_force_threshold_is_independently_required():
    def transverse_force(coordinates):
        return float(coordinates[1]), np.array([0.0, 0.5, 0.0])

    images = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
    ]
    strict_rms = NEB(images, k_spring=0.0, climbing=False).run(
        transverse_force,
        fmax_tol=0.6,
        frms_tol=0.2,
        maxiter=1,
        dt=1.0e-8,
    )
    loose_rms = NEB(images, k_spring=0.0, climbing=False).run(
        transverse_force,
        fmax_tol=0.6,
        frms_tol=0.3,
        maxiter=1,
        dt=1.0e-8,
    )

    assert strict_rms["fmax"] == pytest.approx(0.5)
    assert strict_rms["frms"] == pytest.approx(0.5 / np.sqrt(3.0))
    assert not strict_rms["converged"]
    assert loose_rms["converged"]


def test_result_data_match_coordinates_after_the_last_fire_move():
    images = [np.array([0.0]), np.array([0.2]), np.array([2.0])]
    result = NEB(images, k_spring=1.0, climbing=False).run(
        _one_dimensional_quadratic,
        fmax_tol=0.0,
        frms_tol=0.0,
        maxiter=1,
        dt=0.1,
        maxmove=1.0,
    )

    assert set(result) >= {
        "images", "energies", "gradients", "fmax", "frms", "converged", "iters"
    }
    for image, energy, gradient in zip(
        result["images"], result["energies"], result["gradients"]
    ):
        expected_energy, expected_gradient = _one_dimensional_quadratic(image)
        assert energy == pytest.approx(expected_energy)
        assert np.allclose(gradient, expected_gradient)


def test_climbing_activation_threshold_must_precede_final_convergence():
    images = [np.array([0.0]), np.array([1.0]), np.array([2.0])]
    neb = NEB(images, k_spring=0.1, climbing=True, climb_fmax=0.005)

    with pytest.raises(ValueError, match="climb_fmax.*greater than or equal"):
        neb.run(
            _one_dimensional_quadratic,
            fmax_tol=0.01,
            frms_tol=0.01,
            maxiter=1,
        )


def test_write_neb_xyz_emits_angstrom_multiframe_path(tmp_path):
    bohr = 1.0 / NEB_UTILS.BOHR_TO_ANGSTROM
    output = tmp_path / "nested" / "final_neb.xyz"

    returned = NEB_UTILS.write_neb_xyz(
        output,
        symbols=["H", "O"],
        images_bohr=[
            np.array([0.0, 0.0, 0.0, bohr, 0.0, 0.0]),
            np.array([[0.0, 0.0, 0.0], [0.0, bohr, 0.0]]),
        ],
        energies=[-75.0, -74.9],
    )

    lines = output.read_text().splitlines()
    assert returned == output
    assert lines[0] == "2"
    assert lines[1] == "NEB image 1/2 energy_hartree=-75.0000000000000000"
    assert lines[4] == "2"
    assert lines[5] == "NEB image 2/2 energy_hartree=-74.9000000000000057"
    assert np.allclose([float(value) for value in lines[3].split()[1:]], [1.0, 0.0, 0.0])
    assert np.allclose([float(value) for value in lines[7].split()[1:]], [0.0, 1.0, 0.0])


def test_write_neb_xyz_rejects_mismatched_energy_count(tmp_path):
    with pytest.raises(ValueError, match="energies must match"):
        NEB_UTILS.write_neb_xyz(
            tmp_path / "bad.xyz",
            symbols=["H"],
            images_bohr=[np.zeros(3), np.ones(3)],
            energies=[0.0],
        )
