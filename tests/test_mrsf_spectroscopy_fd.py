import importlib.util
from pathlib import Path
import sys
import types

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "pyoqp/oqp/library"


def load_modules():
    oqp_package = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    oqp_package.__path__ = [str(ROOT / "pyoqp/oqp")]
    library_package = sys.modules.setdefault(
        "oqp.library", types.ModuleType("oqp.library")
    )
    library_package.__path__ = [str(LIBRARY)]

    vibronic_name = "oqp.library.vibronic"
    if vibronic_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            vibronic_name, LIBRARY / "vibronic.py"
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[vibronic_name] = module
        spec.loader.exec_module(module)

    name = "oqp.library.mrsf_spectroscopy_fd"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            name, LIBRARY / "mrsf_spectroscopy_fd.py"
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


@pytest.fixture(scope="module")
def fd():
    return load_modules()


class FakeStates:
    na = 4
    nb = 2

    def __init__(self):
        self.energies = np.array([0.0, 0.5, 1.0, 2.0, 3.0])

    def amplitude_matrix(self, root):
        assert root == 0
        # rows C,C,O,O; columns O,O,V,V.  Every physical class is present.
        return np.array(
            [
                [1.0, 0.0, 2.0, 0.0],
                [0.0, 1.0, 0.0, 2.0],
                [3.0, 0.0, 4.0, 0.0],
                [0.0, 3.0, 0.0, 4.0],
            ]
        )

    def transition_dipole(self, first, second):
        pair = tuple(sorted((first, second)))
        return {
            (0, 1): np.array([1.0, 0.0, 0.0]),
            (0, 2): np.array([0.0, 1.0, 0.0]),
            (0, 3): np.zeros(3),
            (0, 4): np.zeros(3),
        }[pair]


def test_request_schema_requires_spin_adapted_two_somo_topology(fd):
    record = {
        "schema_version": 1,
        "electronic_method": "MRSF-TDDFT",
        "electronic_state": "MRSF S0",
        "electronic_state_role": "target_excited_state",
        "state_index": 1,
        "response_representation": fd.SPIN_ADAPTED_LAYOUT,
        "coordinate_basis": "excited_normal",
        "coordinate_unit": "sqrt(amu)*bohr",
        "coordinate_phase_convention": "fixed synthetic phases",
        "normal_modes": [[1.0, 0.0, 0.0]],
        "displacement": 0.001,
    }
    request = fd.MRSFPropertyFDRequest.from_dict(record)
    assert request.state_index == 1
    assert request.displacement.tolist() == [0.001]

    bad = dict(record, response_representation="determinant_projection")
    with pytest.raises(ValueError, match="spin-adapted"):
        fd.MRSFPropertyFDRequest.from_dict(bad)

    bad_steps = dict(record, finite_difference_step_scales=[1.0, 0.5])
    with pytest.raises(ValueError, match="exactly"):
        fd.MRSFPropertyFDRequest.from_dict(bad_steps)

    weak_tracking = dict(record, minimum_state_overlap=0.98)
    with pytest.raises(ValueError, match="projector treatment"):
        fd.MRSFPropertyFDRequest.from_dict(weak_tracking)

    with pytest.raises(ValueError, match="must be real"):
        fd._finite_real("ordinary MRSF quantity", [1.0 + 1.0e-4j])


def test_response_weights_retain_all_co_ov_cv_oo_classes(fd):
    weights = fd.response_block_weights(FakeStates(), 0)
    assert set(weights) == {"CO", "OV", "CV", "OO"}
    assert sum(weights.values()) == pytest.approx(1.0)
    assert all(value > 0.0 for value in weights.values())


def test_normal_coordinate_unit_requires_mass_normalized_modes(fd):
    masses = np.array([1.0, 4.0])
    modes = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.0, 0.0],
        ]
    )
    assert fd._validate_mass_normalized_modes(modes, masses) < 1.0e-15
    with pytest.raises(ValueError, match="not mass-normalized"):
        fd._validate_mass_normalized_modes(2.0 * modes, masses)


def test_full_state_dipole_contracts_complete_mrsf_density_and_nuclei(fd):
    states = types.SimpleNamespace(
        n_elec=3,
        S=np.ones((1, 1)),
        R=[np.array([[0.2]]), np.zeros((1, 1)), np.zeros((1, 1))],
        state_density_ao=lambda root: np.array([[3.0]]) if root == 0 else None,
    )
    mol = types.SimpleNamespace(
        config={"input": {"charge": 0}},
        data={"ecp_zn": np.array([0, 0])},
        get_atoms=lambda: np.array([1, 2]),
        get_mass=lambda: np.array([1.0, 1.0]),
        get_system=lambda: np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    )
    # Nuclear moment about the COM is +1; the full three-electron state
    # density contributes -3*0.2 = -0.6.
    np.testing.assert_allclose(
        fd.full_mrsf_state_dipole(states, 0, mol), [0.4, 0.0, 0.0]
    )


def test_truncated_sos_tensor_and_tail_convergence_gate(fd):
    tensor, diagnostics = fd.truncated_sos_polarizability(
        FakeStates(),
        0,
        tail_states=2,
        tail_relative_tolerance=1.0e-8,
        minimum_gap_hartree=1.0e-5,
    )
    np.testing.assert_allclose(tensor, np.diag([4.0, 2.0, 0.0]))
    assert diagnostics["backend"] == "truncated_mrsf_state_sos"
    assert diagnostics["analytic"] is False
    assert diagnostics["finite_field"] is False

    states = FakeStates()
    original = states.transition_dipole

    def tail_dipole(first, second):
        if tuple(sorted((first, second))) == (0, 4):
            return np.array([10.0, 0.0, 0.0])
        return original(first, second)

    states.transition_dipole = tail_dipole
    with pytest.raises(ValueError, match="not converged"):
        fd.truncated_sos_polarizability(
            states,
            0,
            tail_states=2,
            tail_relative_tolerance=0.01,
            minimum_gap_hartree=1.0e-5,
        )

    unsorted = FakeStates()
    unsorted.energies = np.array([0.0, 0.5, 2.0, 1.0, 3.0])
    with pytest.raises(ValueError, match="nondecreasing order"):
        fd.truncated_sos_polarizability(
            unsorted,
            0,
            tail_states=2,
            tail_relative_tolerance=0.05,
            minimum_gap_hartree=1.0e-5,
        )


def test_normal_coordinate_driver_produces_ir_and_raman_tensors(fd):
    step = 1.0e-3
    request = fd.MRSFPropertyFDRequest.create(
        electronic_method="MRSF-TDDFT",
        electronic_state="MRSF S0",
        state_index=1,
        normal_modes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        displacement=step,
        coordinate_phase_convention="synthetic fixed modes",
        minimum_state_overlap=0.99,
        minimum_tracking_margin=0.05,
    )

    def evaluator(mode, sign, geometry):
        displacement_value = float(np.asarray(geometry).reshape(-1)[mode])
        dipole = np.zeros(3)
        dipole[mode] = displacement_value * (mode + 2.0)
        polar = np.eye(3) * displacement_value * (mode + 1.0)
        return fd.MRSFTrackedPropertySnapshot.create(
            state_dipole_au=dipole,
            polarizability_au=polar,
            matched_overlap=0.995,
            tracking_margin=0.2,
            phase_to_central=-1.0 if sign < 0 else 1.0,
            raw_root_index=mode,
            state_energy_hartree=-75.0 + displacement_value,
            target_multiplicity=1,
            expected_s2=0.0,
            state_irrep="a1",
            response_block_weights={
                "CO": 0.1,
                "OV": 0.2,
                "CV": 0.3,
                "OO": 0.4,
            },
            polarizability_diagnostics={"backend": "synthetic_truncated_sos"},
            provenance={"fixture": True},
        )

    result = fd.assemble_mrsf_spectroscopy_derivatives(
        request,
        evaluator,
        [0.0, 0.0, 0.0],
        electronic_method="MRSF-TDDFT",
    )
    np.testing.assert_allclose(
        result.dipole_derivative.values,
        [[2.0, 0.0], [0.0, 3.0], [0.0, 0.0]],
    )
    np.testing.assert_allclose(
        result.polarizability_derivative.values[:, :, 0], np.eye(3)
    )
    np.testing.assert_allclose(
        result.polarizability_derivative.values[:, :, 1], 2.0 * np.eye(3)
    )
    assert result.infrared.intensities_km_mol.shape == (2,)
    np.testing.assert_allclose(
        result.infrared.intensities_km_mol,
        974.88011 * np.array([4.0, 9.0]),
    )
    assert result.raman.activities_au.shape == (2,)
    assert len(result.displacement_records) == 12
    assert result.provenance["dipole_fd_convergence"]["accepted_step_scale"] == 0.25
    assert result.displacement_records[0]["state_dipole_unit"] == "e*bohr"
    assert result.displacement_records[0]["polarizability_unit"] == "bohr^3"
    assert result.provenance["response_representation"] == fd.SPIN_ADAPTED_LAYOUT


def test_incomplete_sos_keeps_ir_but_publishes_no_raman(fd):
    request = fd.MRSFPropertyFDRequest.create(
        electronic_method="MRSF-TDHF",
        electronic_state="MRSF S1",
        state_index=2,
        normal_modes=[[1.0, 0.0, 0.0]],
        displacement=0.001,
        coordinate_phase_convention="synthetic",
    )

    def evaluator(mode, sign, geometry):
        del mode, sign
        displacement_value = float(np.asarray(geometry).reshape(-1)[0])
        return fd.MRSFTrackedPropertySnapshot.create(
            state_dipole_au=[displacement_value, 0.0, 0.0],
            polarizability_au=None,
            matched_overlap=0.99,
            tracking_margin=0.1,
            phase_to_central=1.0,
            raw_root_index=1,
            state_energy_hartree=-10.0,
            target_multiplicity=1,
            expected_s2=0.0,
            state_irrep=None,
            response_block_weights={
                "CO": 0.25,
                "OV": 0.25,
                "CV": 0.25,
                "OO": 0.25,
            },
            polarizability_diagnostics={
                "status": "unavailable",
                "replacement_used": False,
            },
        )

    result = fd.assemble_mrsf_spectroscopy_derivatives(
        request,
        evaluator,
        [0.0, 0.0, 0.0],
        electronic_method="MRSF-TDHF",
    )
    assert result.infrared.intensities_km_mol.size == 1
    assert result.polarizability_derivative is None
    assert result.raman is None


def test_unconverged_polarizability_fd_keeps_independent_ir(fd):
    request = fd.MRSFPropertyFDRequest.create(
        electronic_method="MRSF-TDHF",
        electronic_state="MRSF S0",
        state_index=1,
        normal_modes=[[1.0, 0.0, 0.0]],
        displacement=0.001,
        coordinate_phase_convention="synthetic",
    )

    def evaluator(mode, sign, geometry):
        del mode, sign
        displacement_value = float(np.asarray(geometry).reshape(-1)[0])
        return fd.MRSFTrackedPropertySnapshot.create(
            state_dipole_au=[displacement_value, 0.0, 0.0],
            polarizability_au=np.eye(3) * (
                displacement_value + 1.0e9 * displacement_value**3
            ),
            matched_overlap=0.99,
            tracking_margin=0.1,
            phase_to_central=1.0,
            raw_root_index=0,
            state_energy_hartree=-10.0,
            target_multiplicity=1,
            expected_s2=0.0,
            state_irrep=None,
            response_block_weights={
                "CO": 0.25,
                "OV": 0.25,
                "CV": 0.25,
                "OO": 0.25,
            },
            polarizability_diagnostics={"status": "synthetic"},
        )

    result = fd.assemble_mrsf_spectroscopy_derivatives(
        request,
        evaluator,
        [0.0, 0.0, 0.0],
        electronic_method="MRSF-TDHF",
    )
    assert result.infrared.intensities_km_mol.size == 1
    assert result.raman is None
    assert "not converged" in result.provenance["polarizability_fd_failure"]


def test_fd_convergence_is_required_for_each_normal_mode(fd):
    request = fd.MRSFPropertyFDRequest.create(
        electronic_method="MRSF-TDHF",
        electronic_state="MRSF S0",
        state_index=1,
        normal_modes=[
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        ],
        displacement=0.001,
        coordinate_phase_convention="synthetic",
    )

    def evaluator(mode, sign, geometry):
        del sign
        displacement_value = float(np.asarray(geometry).reshape(-1)[mode])
        dipole = np.zeros(3)
        if mode == 0:
            # This large, exactly linear derivative would dominate the former
            # aggregate Frobenius-norm convergence check.
            dipole[0] = 1.0e6 * displacement_value
        else:
            # Central derivatives at h, h/2, and h/4 are 100001, 25001, and
            # 6251.  The individual mode is plainly unconverged even though its
            # medium-to-fine change is below 5% of the aggregate derivative.
            dipole[1] = displacement_value + 1.0e8 * displacement_value**3
        return fd.MRSFTrackedPropertySnapshot.create(
            state_dipole_au=dipole,
            polarizability_au=None,
            matched_overlap=0.995,
            tracking_margin=0.2,
            phase_to_central=1.0,
            raw_root_index=0,
            state_energy_hartree=-10.0,
            target_multiplicity=1,
            expected_s2=0.0,
            state_irrep=None,
            response_block_weights={
                "CO": 0.25,
                "OV": 0.25,
                "CV": 0.25,
                "OO": 0.25,
            },
        )

    with pytest.raises(ValueError, match=r"normal modes 2 "):
        fd.assemble_mrsf_spectroscopy_derivatives(
            request,
            evaluator,
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            electronic_method="MRSF-TDHF",
        )


def test_finite_field_request_fails_without_substituting_rohf_cphf(fd):
    request = fd.MRSFPropertyFDRequest.create(
        electronic_method="MRSF-TDDFT",
        electronic_state="MRSF S0",
        state_index=1,
        normal_modes=[[1.0, 0.0, 0.0]],
        displacement=0.001,
        coordinate_phase_convention="synthetic",
        polarizability_backend="finite_field",
    )
    with pytest.raises(NotImplementedError, match="no ROHF/CPHF"):
        fd.assemble_mrsf_spectroscopy_derivatives(
            request,
            lambda *_: None,
            [0.0, 0.0, 0.0],
            electronic_method="MRSF-TDDFT",
        )


def test_displaced_openqp_child_config_is_configparser_safe(fd):
    evaluator = object.__new__(fd.OpenQPMRSFNormalModeEvaluator)
    evaluator.config = {
        "input": {"runtype": "hess", "method": "tdhf"},
        "guess": {"type": "json", "save_mol": True},
        "properties": {"back_door": False, "scf_prop": ["dipole", "quadrupole"]},
        "nac": {"align": "no"},
        "symmetry": {
            "use_integral_symmetry": True,
            "move_to_standard_frame": True,
        },
        "hess": {"temperature": np.array([100.0, 298.15])},
    }
    child = evaluator._child_config()
    assert all(
        isinstance(value, str)
        for section in child.values()
        for value in section.values()
    )
    assert child["guess"]["save_mol"] == "false"
    assert child["properties"]["scf_prop"] == ""
    assert child["hess"]["temperature"] == "100.0,298.15"


def test_mrsf_hessian_never_reaches_generic_rohf_cphf_property_path():
    source = (LIBRARY / "single_point.py").read_text(encoding="utf-8")
    function = source[source.index("    def _compute_vibrational_intensities") :]
    function = function[: function.index("    def analytical_hess")]
    guard = function.index("if td_type == 'mrsf':")
    generic = function.index("displacement = 1.0e-3", guard)
    dedicated = function.index("run_openqp_mrsf_spectroscopy_fd", guard)
    assert guard < dedicated < generic
    branch = function[guard:generic]
    assert "generic_rohf_cphf_fallback': False" in branch
    assert "return" in branch
    assert "cphf_static_polarizability" not in branch


def test_implementation_contains_no_determinant_projection():
    source = (LIBRARY / "mrsf_spectroscopy_fd.py").read_text(encoding="utf-8").lower()
    assert "slater" not in source
    assert "determinant_projection" not in source
    assert "electric_dipole_au" not in source
    assert "cphf_static_polarizability" not in source
    assert "track_isolated_mrsf_hessian_root" in source
    assert "two_somo_spin_adapted_co_ov_cv_oo" in source
