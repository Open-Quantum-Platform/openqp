"""Native ODP umbrella equations, invariants, restart, and provenance gates."""

import json
import math
import os
import shutil
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from oqp_module_isolation import child_env

ROOT = Path(__file__).resolve().parents[1]


def test_odp_kernel_is_fortran_resident_and_c_interoperable():
    source = (ROOT / "source" / "modules" / "odp_umbrella.F90").read_text()
    header = (ROOT / "include" / "oqp.h").read_text()
    driver = (ROOT / "pyoqp" / "oqp" / "library" / "odp.py").read_text()

    assert 'bind(C, name="oqp_odp_umbrella_evaluate")' in source
    assert "int oqp_odp_umbrella_evaluate(" in header
    assert "oqp.oqp_odp_umbrella_evaluate(" in driver
    assert "gradient = gradient - odp['force']" in (
        ROOT / "pyoqp" / "oqp" / "library" / "namd.py"
    ).read_text()
    assert "abs(" not in source.split("xi =", 1)[1].splitlines()[0]


def test_odp_configuration_requires_explicit_metric_and_continuous_cvs():
    from oqp.library.odp import ODPUmbrella, parse_odp_cv_specification

    types, atoms, labels = parse_odp_cv_specification(
        "distance(1,2); asymmetric_distance(1,2,1,3); angle(4,5,6)"
    )
    assert types.tolist() == [1, 2, 3]
    assert atoms.tolist() == [[0, 1, -1, -1], [0, 1, 0, 2], [3, 4, 5, -1]]
    assert labels == (
        "distance(1,2)", "asymmetric_distance(1,2,1,3)", "angle(4,5,6)")
    with pytest.raises(ValueError, match="scale requires exactly 3"):
        ODPUmbrella({
            "enabled": True, "cv": "; ".join(labels), "scale": [1.0],
            "reference_r": [1.0, 0.0, 1.0],
            "reference_p": [2.0, 0.5, 2.0], "k_parallel": 0.1,
        })
    with pytest.raises(ValueError, match="unsupported CV"):
        parse_odp_cv_specification("nearest_atom(1,2,3)")
    for specification in (
        "asymmetric_distance(1,2,1,2)",
        "asymmetric_distance(1,2,2,1)",
    ):
        with pytest.raises(ValueError, match="atom pairs must differ"):
            parse_odp_cv_specification(specification)


def test_runtime_checker_rejects_enabled_odp_outside_namd():
    from oqp.utils.input_checker import check_input_values

    config = {
        "input": {
            "system": "\n1 0.0 0.0 0.0\n1 0.0 0.0 1.4",
            "charge": 0, "basis": "sto-3g", "functional": "",
            "method": "hf", "runtype": "energy",
        },
        "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
        "odp": {"enabled": True},
    }
    report = check_input_values(config, raise_error=False, emit=False)
    assert any(item.path == "odp.enabled" for item in report.errors)

    config["odp"]["enabled"] = False
    report = check_input_values(config, raise_error=False, emit=False)
    assert not any(item.path == "odp.enabled" for item in report.errors)


def test_odp_accepts_langevin_nvt_in_legacy_and_semantic_inputs():
    from oqp.utils.input_checker import CheckReport, _check_runtype
    from oqp.utils.oqp_input import parse_canonical_oqp

    config = {
        "input": {"method": "tdhf", "runtype": "namd"},
        "md": {"ensemble": "nvt", "thermostat": "langevin"},
        "odp": {"enabled": True},
    }
    report = CheckReport()
    _check_runtype(config, report)
    assert not any(item.path in ("md.ensemble", "odp.enabled")
                   for item in report.errors)

    parsed = parse_canonical_oqp(
        'mrsf(nstate=2)/bhhlyp/sto-3g '
        'namd(S1,scheme=Overlap,nstep=2,ensemble=nvt,thermostat=langevin) '
        'odp(enabled=true,cv="distance(1,2)",scale="1",'
        'reference_r="1",reference_p="2",k_parallel=0.1) '
        'geom="h2.xyz"'
    )
    assert any(call.name == "odp" for call in parsed.modifiers)


def test_namd_has_no_odp_ensemble_or_boundary_rejection():
    source = (ROOT / "pyoqp" / "oqp" / "library" / "namd.py").read_text()
    assert "_validate_odp_boundary_conditions" not in source
    assert "same-spin NVE NAMD only" not in source
    # every ODP evaluation in the NAMD drivers goes through the box hook
    assert "self.odp.evaluate(coordinates, box=self._odp_box())" in source
    assert source.count("self.odp.evaluate(") == 1


def test_wham_rejects_mismatched_nvt_sampling_temperature(tmp_path, monkeypatch):
    from oqp.library import namd as namd_module
    from oqp.library.odp import odp_wham

    trajectory = tmp_path / "window.namd.trj"
    trajectory.write_bytes(b"temperature-provenance")
    monkeypatch.setattr(
        namd_module, "read_odp_wham_series",
        lambda _path: {
            "ensemble": "NVT",
            "thermostat_temperature_kelvin": 310.0,
        },
    )
    with pytest.raises(ValueError, match="sampled at 310 K"):
        odp_wham([trajectory], 300.0)


def test_built_odp_randomized_fd_invariance_signed_progress_and_singularities():
    script = r'''
import json
import math
import numpy as np
from oqp.library.odp import ODPUmbrella

section = {
    "enabled": True,
    "cv": "distance(1,2); asymmetric_distance(1,2,1,3); angle(4,5,6)",
    "scale": [0.7, 0.6, 0.9],
    "reference_r": [1.0, 0.1, 1.3],
    "reference_p": [1.8, -0.5, 2.0],
    "center": 0.35,
    "k_parallel": 0.08,
    "k_perpendicular": 0.03,
    "window": 7,
}
odp = ODPUmbrella(section)

def raw_cv(x):
    d12 = np.linalg.norm(x[0] - x[1])
    asym = d12 - np.linalg.norm(x[0] - x[2])
    a = x[3] - x[4]
    b = x[5] - x[4]
    cosine = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))
    angle = math.acos(np.clip(cosine, -1.0, 1.0))
    return np.array([d12, asym, angle])

def oracle_energy(x):
    raw = raw_cv(x)
    scaled = np.array(section["scale"])*raw
    reactant = np.array(section["scale"])*np.array(section["reference_r"])
    direction = np.array(section["scale"])*(np.array(section["reference_p"])
                                              - np.array(section["reference_r"]))
    displacement = scaled - reactant
    xi = np.dot(displacement, direction)/np.dot(direction, direction)
    perpendicular = displacement - xi*direction
    energy_parallel = 0.5*section["k_parallel"]*(xi - section["center"])**2
    energy_perpendicular = 0.5*section["k_perpendicular"]*np.dot(
        perpendicular, perpendicular)
    return energy_parallel + energy_perpendicular, xi, raw, perpendicular

rng = np.random.default_rng(20260803)
base = np.array([
    [0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.2, 1.1, 0.3],
    [3.0, 0.0, 0.0], [4.0, 0.2, 0.0], [3.2, 1.1, 0.4],
])
max_energy_error = 0.0
max_force_error = 0.0
max_translation_force = 0.0
max_torque = 0.0
max_rotation_energy_error = 0.0
max_rotation_force_error = 0.0
h = 1.0e-6
for _ in range(16):
    xyz = base + rng.normal(scale=0.12, size=base.shape)
    native = odp.evaluate(xyz)
    expected, expected_xi, expected_raw, expected_perp = oracle_energy(xyz)
    max_energy_error = max(max_energy_error, abs(native["energy"] - expected),
                           abs(native["xi"] - expected_xi),
                           np.max(np.abs(native["cv_raw"] - expected_raw)),
                           np.max(np.abs(native["cv_perpendicular"] - expected_perp)))
    fd_force = np.zeros_like(xyz)
    for atom in range(len(xyz)):
        for axis in range(3):
            plus = xyz.copy(); plus[atom, axis] += h
            minus = xyz.copy(); minus[atom, axis] -= h
            fd_force[atom, axis] = -(odp.evaluate(plus)["energy"]
                                     - odp.evaluate(minus)["energy"])/(2*h)
    native = odp.evaluate(xyz)
    max_force_error = max(max_force_error,
                          np.max(np.abs(native["force"] - fd_force)))
    max_translation_force = max(max_translation_force,
                                np.linalg.norm(native["force"].sum(axis=0)))
    max_torque = max(max_torque,
                     np.linalg.norm(np.cross(xyz, native["force"]).sum(axis=0)))

    matrix = rng.normal(size=(3, 3))
    q, _ = np.linalg.qr(matrix)
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    shift = rng.normal(size=3)
    transformed = xyz @ q.T + shift
    rotated = odp.evaluate(transformed)
    max_rotation_energy_error = max(
        max_rotation_energy_error, abs(rotated["energy"] - native["energy"]))
    max_rotation_force_error = max(
        max_rotation_force_error,
        np.max(np.abs(rotated["force"] - native["force"] @ q.T)))

distance = ODPUmbrella({
    "enabled": True, "cv": "distance(1,2)", "scale": [0.5],
    "reference_r": [2.0], "reference_p": [3.0], "center": 0.5,
    "k_parallel": 0.1, "k_perpendicular": 0.0,
})
negative = distance.evaluate(np.array([[0., 0., 0.], [1., 0., 0.]]))["xi"]
beyond = distance.evaluate(np.array([[0., 0., 0.], [4., 0., 0.]]))["xi"]

singular_distance = None
try:
    distance.evaluate(np.zeros((2, 3)))
except RuntimeError as exc:
    singular_distance = str(exc)
singular_angle = None
angle = ODPUmbrella({
    "enabled": True, "cv": "angle(1,2,3)", "scale": [1.0],
    "reference_r": [1.0], "reference_p": [2.0], "center": 0.5,
    "k_parallel": 0.1, "k_perpendicular": 0.0,
})
try:
    angle.evaluate(np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]))
except RuntimeError as exc:
    singular_angle = str(exc)

print("ODP=" + json.dumps({
    "max_energy_error": max_energy_error,
    "max_force_error": max_force_error,
    "max_translation_force": max_translation_force,
    "max_torque": max_torque,
    "max_rotation_energy_error": max_rotation_energy_error,
    "max_rotation_force_error": max_rotation_force_error,
    "negative": negative, "beyond": beyond,
    "singular_distance": singular_distance,
    "singular_angle": singular_angle,
}))
'''
    env = child_env()
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("ODP=")), None)
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("ODP="))
    assert values["max_energy_error"] < 2.0e-14
    assert values["max_force_error"] < 2.0e-9
    assert values["max_translation_force"] < 2.0e-14
    assert values["max_torque"] < 2.0e-13
    assert values["max_rotation_energy_error"] < 2.0e-14
    assert values["max_rotation_force_error"] < 2.0e-13
    assert values["negative"] == pytest.approx(-1.0)
    assert values["beyond"] == pytest.approx(2.0)
    assert "singular distance" in values["singular_distance"]
    assert "collinear angle" in values["singular_angle"]


def test_built_odp_periodic_cvs_use_minimum_images():
    script = r"""
import json
import numpy as np
from oqp.library.odp import ODPUmbrella
from oqp.library.namd import NAMD_QMMM

section = {
    "enabled": True,
    "cv": "distance(2,5); asymmetric_distance(5,2,5,7); angle(2,5,7)",
    "scale": [0.8, 0.6, 1.1],
    "reference_r": [3.4, -0.6, 1.7],
    "reference_p": [2.1, 0.9, 2.2],
    "center": 0.37,
    "k_parallel": 0.07,
    "k_perpendicular": 0.02,
    "window": 3,
}
odp = ODPUmbrella(section)
box = np.array([23.0, 27.0, 31.0])
rng = np.random.default_rng(20261020)
# A compact, non-symmetric cluster of CV atoms (2,5,7 in 1-based labels)
# among spectator atoms; CV atoms are then pushed through different faces.
whole = rng.uniform(-14.0, 14.0, size=(8, 3))
whole[1] = [10.9, -12.6, 14.8]
whole[4] = [12.7, -11.2, 16.1]
whole[6] = [11.3, -13.9, 17.4]
max_energy = max_force = max_xi = max_fd = 0.0
for trial in range(12):
    jitter = rng.normal(scale=0.15, size=whole.shape)
    reference_xyz = whole + jitter
    wrapped = reference_xyz.copy()
    for atom in (1, 4, 6):
        wrapped[atom] += box*rng.integers(-2, 3, size=3)
    for atom in (0, 2, 3, 5, 7):
        wrapped[atom] += box*rng.integers(-1, 2, size=3)
    expected = odp.evaluate(reference_xyz)
    periodic = odp.evaluate(wrapped, box=box)
    max_energy = max(max_energy, abs(periodic["energy"] - expected["energy"]))
    max_xi = max(max_xi, abs(periodic["xi"] - expected["xi"]),
                 float(np.max(np.abs(periodic["cv_raw"] - expected["cv_raw"]))))
    max_force = max(max_force,
                    float(np.max(np.abs(periodic["force"] - expected["force"]))))
    h = 1.0e-6
    for atom, axis in ((1, 0), (4, 2), (6, 1)):
        plus = wrapped.copy(); plus[atom, axis] += h
        minus = wrapped.copy(); minus[atom, axis] -= h
        fd = -(odp.evaluate(plus, box=box)["energy"]
               - odp.evaluate(minus, box=box)["energy"])/(2*h)
        max_fd = max(max_fd, abs(fd - periodic["force"][atom, axis]))

diagonal = np.zeros((8, 3))
diagonal[1] = [0.0, 0.0, 0.0]
diagonal[4] = [0.4*23.0, 0.4*23.0, 0.0]
diagonal[6] = [0.4*23.0 + 1.3, 0.4*23.0 - 0.7, 0.9]
cube = np.array([23.0, 23.0, 23.0])
try:
    odp.evaluate(diagonal, box=cube)
    diagonal_ok = True
except RuntimeError:
    diagonal_ok = False

# the reviewer's angle: two 4.9-bohr arms about vertex 2 in a 10-bohr cell;
# imaging atom 3 against atom 1 would make the second arm 5.1 bohr
angle = ODPUmbrella({
    "enabled": True, "cv": "angle(1,2,3)", "scale": [1.0],
    "reference_r": [1.0], "reference_p": [2.0], "center": 0.5,
    "k_parallel": 0.1, "k_perpendicular": 0.0})
small = np.array([10.0, 10.0, 10.0])
wrapped_angle = np.array([[0.0, 0.0, 0.0], [4.9, 1.0, 0.0], [-0.2, 2.5, 0.7]])
unwrapped_angle = wrapped_angle.copy()
unwrapped_angle[2, 0] += 10.0
angle_periodic = angle.evaluate(wrapped_angle, box=small)
angle_reference = angle.evaluate(unwrapped_angle)
angle_error = abs(angle_periodic["energy"] - angle_reference["energy"]) + float(
    np.max(np.abs(angle_periodic["force"] - angle_reference["force"])))

def image(vector, cell):
    return vector - cell*np.round(vector/cell)

# A ring of distance CVs that winds around the cell (review case: x = 0, 4,
# -4 in a 10-bohr box; minimum-image lengths 4, 2, 4) has no common Cartesian
# embedding, yet every CV is unambiguous: each must take its own
# minimum-image value, and the folded force must be the energy derivative.
triangle = ODPUmbrella({
    "enabled": True, "cv": "distance(1,2); distance(2,3); distance(3,1)",
    "scale": [1.0, 0.7, 1.3], "reference_r": [3.0, 1.5, 3.5],
    "reference_p": [4.5, 2.5, 4.0], "center": 0.4, "k_parallel": 0.1,
    "k_perpendicular": 0.03})
ring_box = np.array([10.0, 10.0, 10.0])
loop = np.array([[0.0, 0.0, 0.0], [4.0, 0.3, -0.2], [-4.0, 0.5, 0.4]])
ring = triangle.evaluate(loop, box=ring_box)
ring_expected = [float(np.linalg.norm(image(loop[j] - loop[i], ring_box)))
                 for i, j in ((0, 1), (1, 2), (2, 0))]
ring_cv = float(np.max(np.abs(ring["cv_raw"] - ring_expected)))
ring_fd = 0.0
for atom in range(3):
    for axis in range(3):
        plus = loop.copy(); plus[atom, axis] += 1.0e-6
        minus = loop.copy(); minus[atom, axis] -= 1.0e-6
        fd = -(triangle.evaluate(plus, box=ring_box)["energy"]
               - triangle.evaluate(minus, box=ring_box)["energy"])/2.0e-6
        ring_fd = max(ring_fd, abs(fd - ring["force"][atom, axis]))

# exactly half the box apart along x: two minimum images, refused
refusal = None
try:
    triangle.evaluate(np.array([[0.0, 0.0, 0.0], [5.0, 0.3, 0.0],
                                [2.0, 1.0, 0.0]]), box=ring_box)
except RuntimeError as exc:
    refusal = str(exc)

# Atoms 2 and 7 are 11.6 bohr apart along x (more than half of 23) but no CV
# uses that pair; each CV, measured along its own pairs 5-2 and 5-7, equals
# the CV of the explicitly unwrapped coordinates.
far_ends = whole.copy()
far_ends[6] = whole[1] + np.array([11.6, 0.0, 0.0])
unwrapped = far_ends.copy()
unwrapped[1] = far_ends[4] + image(far_ends[1] - far_ends[4], box)
unwrapped[6] = far_ends[4] + image(far_ends[6] - far_ends[4], box)
far_ends_shifted = far_ends.copy()
far_ends_shifted[1] += box*np.array([1, 0, -1])
far_ends_shifted[6] += box*np.array([-2, 1, 0])
far_ends_cv = float(np.max(np.abs(
    odp.evaluate(far_ends_shifted, box=box)["cv_raw"]
    - odp.evaluate(unwrapped)["cv_raw"])))

# Two unconnected distances in a 10-bohr cube: each pair is imaged on its own
# (imaging atom 4 against atom 1 put the 3-4 pair 9.6 bohr apart).
apart = ODPUmbrella({"enabled": True, "cv": "distance(1,2); distance(3,4)",
                     "scale": [1.0, 1.0], "reference_r": [1.0, 1.0],
                     "reference_p": [2.0, 2.0], "center": 0.5,
                     "k_parallel": 0.1})
apart_cv = apart.evaluate(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                                    [4.9, 0.0, 0.0], [5.3, 0.0, 0.0]]),
                          box=np.array([10.0, 10.0, 10.0]))["cv_raw"].tolist()

class Driver:
    def __init__(self, cell):
        self.cell = cell
    def _box_lengths_bohr(self):
        return self.cell

namd = NAMD_QMMM.__new__(NAMD_QMMM)
namd.odp = odp
namd.driver = Driver(box)
namd.r_all = wrapped
force, energy = namd._apply_odp_to_force_energy(np.zeros_like(wrapped), -2.5)
direct = odp.evaluate(wrapped, box=box)
cluster = NAMD_QMMM.__new__(NAMD_QMMM)
cluster.odp = odp
cluster.driver = Driver(None)
cluster.r_all = wrapped
cluster_force, cluster_energy = cluster._apply_odp_to_force_energy(
    np.zeros_like(wrapped), -2.5)
raw = odp.evaluate(wrapped)
print("PBC=" + json.dumps({
    "max_energy": max_energy, "max_xi": max_xi, "max_force": max_force,
    "max_fd": max_fd, "refusal": refusal, "diagonal_ok": diagonal_ok,
    "ring_cv": ring_cv, "ring_fd": ring_fd,
    "angle_error": angle_error,
    "far_ends_cv": far_ends_cv, "apart_cv": apart_cv,
    "namd_energy": energy - (-2.5 + direct["energy"]),
    "namd_force": float(np.max(np.abs(force - direct["force"]))),
    "namd_unbiased": namd._unbiased_potential_energy,
    "cluster_energy": cluster_energy - (-2.5 + raw["energy"]),
    "periodic_differs_from_raw": abs(direct["energy"] - raw["energy"]),
}))
"""
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=child_env(),
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("PBC=")), None)
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("PBC="))
    assert values["max_energy"] < 1.0e-12
    assert values["max_xi"] < 1.0e-12
    assert values["max_force"] < 1.0e-12
    assert values["max_fd"] < 2.0e-8
    assert "exactly half the box length" in values["refusal"]
    # a ring of CVs winding around the cell: each CV at its own minimum image
    assert values["ring_cv"] < 1.0e-12
    assert values["ring_fd"] < 2.0e-8
    assert values["angle_error"] < 1.0e-12
    assert values["far_ends_cv"] < 1.0e-12
    assert abs(values["apart_cv"][0] - 1.0) < 1.0e-12
    assert abs(values["apart_cv"][1] - 0.4) < 1.0e-12
    # a diagonal pair separation of (0.4 L, 0.4 L, 0) has a unique minimum
    # image although its norm exceeds half the box
    assert values["diagonal_ok"] is True
    assert values["namd_energy"] == 0.0
    assert values["namd_force"] == 0.0
    assert values["namd_unbiased"] == -2.5
    assert values["cluster_energy"] == 0.0
    # the wrapped frame really crosses faces: without imaging the bias differs
    assert values["periodic_differs_from_raw"] > 1.0e-3


def test_odp_packed_trj_and_restart_preserve_wham_provenance(tmp_path):
    script = r'''
import json
import os
import shutil
from types import SimpleNamespace
import numpy as np
import oqp.library.namd as namd_module
from oqp.library.namd import NAMD, read_namd_trajectory, read_odp_wham_series
from oqp.library.odp import ODPUmbrella
namd_module.dump_log = lambda *_args, **_kwargs: None

root = os.environ["OQP_ODP_TEST_ROOT"]
class Mol:
    log = os.path.join(root, "job.log")
    oqp_canonical_input = (
        "mrsf(nstate=2)/bhhlyp/6-31g*\n"
        "namd(S1,scheme=Overlap,nstep=2,dt=0.5)\n"
        "odp(enabled=true,cv=\"distance(1,2)\",scale=\"1.0\","
        "reference_r=\"1.0\",reference_p=\"2.0\",center=0.5,"
        "k_parallel=0.1,window=4)\ngeom=\"h2.xyz\"\n"
    )
    config = {
        "input": {"method": "tdhf", "functional": "bhhlyp", "basis": "6-31g*"},
        "tdhf": {"type": "mrsf", "nstate": 2, "tlf": 2},
        "md": {"substep": 2, "decoherence": "off", "edc_c": 0.1,
               "thrshe": 1.0, "tdc": "fd", "trivial": False,
               "trivial_thresh": 0.5, "first_hop_step": 2},
    }
    data = {"OQP::td_energies": np.array([0.0, 0.1])}
    @staticmethod
    def get_atoms():
        return np.array([1, 1])
    @staticmethod
    def get_mass():
        return np.array([1.0, 1.0])
    @staticmethod
    def get_state_tracking():
        return None
    def put_data(self, data):
        self.loaded = data

d = NAMD.__new__(NAMD)
d.mol = Mol(); d.nstate = 2; d.dt_fs = 0.5; d.dt_adaptive = False
d._t_fs = 0.0; d.active = 1; d.coef = np.array([1.+0j, 0.+0j])
d.vel = np.zeros((2, 3)); d.trajectory_interval = 1
d.mass = np.ones(2); d.velocity_source = "zero"; d.ensemble = "nve"
d.thermostat = "off"; d.thermostat_temperature = 300.0
d.thermostat_friction = 1.0
d.trajectory_file = os.path.join(root, "odp.namd.trj")
d.restart_file = os.path.join(root, "odp.restart.npz")
d.restart_manifest_file = os.path.join(root, "restart.oqp")
d.restart_interval = 1; d.restart_requested = False
d.seed = 1; d.rng_stream = 0; d.init_temp = 300.0; d._rng_step = 1
d._restart_system_identity = {
    "kind": "test", "natom": 2, "sha256": "odp-system"
}
d._wham_system_identity = {
    "kind": "test", "natom": 2, "sha256": "odp-topology"
}
d._last_hop_random = 0.25; d._last_state_overlap = None
d._last_overlap_tdc = None; d._nacme_reference_tdc = None
d._nacme_reference_mask = None; d._nacme_reference_source = 0
d._nacme_gate_last = None; d._nve_gate_last = None
d._nacme_gate_failures = 0; d._nve_gate_failures = 0
d._ba_energy_left = d._ba_energy_center = d._ba_tdc_left = d._ba_dt_left = None
d._nve_reference_energy = d._nve_previous_energy = None
d.prev_xyz = np.zeros(6); d.prev_data = {"x": np.array([1.0])}
d.odp = ODPUmbrella({
    "enabled": True, "cv": "distance(1,2)", "scale": [1.0],
    "reference_r": [1.0], "reference_p": [2.0], "center": 0.5,
    "k_parallel": 0.1, "k_perpendicular": 0.0, "window": 4,
})
coords = np.array([[0., 0., 0.], [1.4, 0., 0.]])
d._odp_last = d.odp.evaluate(coords)
d._unbiased_potential_energy = -1.2
d._write_md_trajectory(1, coords, -1.2 + d._odp_last["energy"], 0.3, False)
d._save_restart(1, coords, np.zeros_like(coords), np.ones_like(coords))
header, records = read_namd_trajectory(d.trajectory_file)
series = read_odp_wham_series(d.trajectory_file)
with np.load(d.restart_file, allow_pickle=False) as saved:
    restart_provenance = json.loads(str(saved["odp_provenance"][0]))
d2 = NAMD.__new__(NAMD)
d2.mol = Mol(); d2.nstate = 2; d2.dt_fs = 0.5; d2.dt_adaptive = False
d2.seed = 1; d2.rng_stream = 0; d2.restart_requested = True
d2.mass = np.ones(2); d2.velocity_source = "zero"; d2.ensemble = "nve"
d2.thermostat = "off"; d2.thermostat_temperature = 300.0
d2.thermostat_friction = 1.0
d2._restart_system_identity = d._restart_system_identity
d2._wham_system_identity = d._wham_system_identity
d2.restart_file = d.restart_file; d2.trajectory_file = d.trajectory_file
d2.odp = ODPUmbrella({
    "enabled": True, "cv": "distance(1,2)", "scale": [1.0],
    "reference_r": [1.0], "reference_p": [2.0], "center": 0.5,
    "k_parallel": 0.1, "k_perpendicular": 0.0, "window": 4,
})
restored = d2._load_restart()
restored_odp = d2.odp.evaluate(restored["coordinates"])
payload = {
    "header": header,
    "xi": float(records["odp_xi"][0]),
    "window": int(records["odp_window"][0]),
    "bias": float(records["odp_bias_hartree"][0]),
    "unbiased": float(records["e_unbiased_pot_hartree"][0]),
    "series_xi": float(series["xi"][0]),
    "restart_provenance": restart_provenance,
    "restart_nuclear_exact": (
        np.array_equal(restored["coordinates"], coords)
        and np.array_equal(restored["velocities"], np.zeros_like(coords))
        and np.array_equal(restored["acceleration"], np.ones_like(coords))
    ),
    "restart_odp_exact": (
        restored_odp["energy"] == d._odp_last["energy"]
        and restored_odp["xi"] == d._odp_last["xi"]
        and np.array_equal(restored_odp["force"], d._odp_last["force"])
    ),
    "manifest": open(d.restart_manifest_file).read(),
}
del records
print("ODP_IO=" + json.dumps(payload))
'''
    env = child_env()
    env["OQP_ODP_TEST_ROOT"] = str(tmp_path)
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("ODP_IO=")), None)
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("ODP_IO="))
    assert values["header"]["ensemble"] == "NVE"
    assert values["header"]["odp"]["window"] == 4
    assert values["header"]["wham"]["reaction_coordinate_field"] == "odp_xi"
    assert values["window"] == 4
    assert values["xi"] == pytest.approx(0.4)
    assert values["series_xi"] == pytest.approx(0.4)
    assert values["unbiased"] == pytest.approx(-1.2)
    assert values["bias"] == pytest.approx(0.0005)
    assert values["restart_provenance"] == values["header"]["odp"]
    assert values["restart_nuclear_exact"] is True
    assert values["restart_odp_exact"] is True
    assert "odp(" in values["manifest"]
    assert "restart=true" in values["manifest"]


def test_python_wham_recovers_synthetic_unbiased_distribution(tmp_path):
    script = r'''
import hashlib
import json
import os
import shutil
import shutil
import struct
import numpy as np
import oqp.library.namd as namd_module
from oqp.library.namd import (
    NAMD_TRAJECTORY_MAGIC, NAMD_TRAJECTORY_SCHEMA_VERSION,
    _namd_trajectory_dtype,
)
from oqp.library.odp import KB_HARTREE_PER_KELVIN, odp_wham

root = os.environ["OQP_ODP_WHAM_ROOT"]
temperature = 300.0
beta = 1.0/(KB_HARTREE_PER_KELVIN*temperature)
k_unbiased = 0.002
k_umbrella = 0.006
centers = [-1.2, -0.6, 0.0, 0.6, 1.2]
rng = np.random.default_rng(24681357)
paths = []
dtype = _namd_trajectory_dtype(1, 1, 1)
system_identity = {
    "method": "synthetic", "charge": 0, "functional": "", "basis": "none",
    "d4": False,
    "scf_type": "rhf", "scf_multiplicity": 1,
    "tdhf_type": "rpa", "tdhf_multiplicity": 1, "nstate": 1, "tlf": 2,
    "pcm": {"enabled": True, "model": "ddpcm", "epsilon": 78.3553},
    "trajectory_representation": "same_spin_adiabatic",
    "system": {"kind": "synthetic", "natom": 1, "sha256": "system-a"},
}
for window, center in enumerate(centers):
    mean = k_umbrella*center/(k_unbiased + k_umbrella)
    sigma = np.sqrt(1.0/(beta*(k_unbiased + k_umbrella)))
    xi = rng.normal(mean, sigma, size=12000)
    umbrella = 0.5*k_umbrella*(xi - center)**2
    unbiased = 0.5*k_unbiased*xi**2
    records = np.zeros(xi.size, dtype=dtype)
    records["step"] = np.arange(xi.size)
    records["odp_window"] = window
    records["odp_xi"] = xi
    records["odp_cv_raw"][:, 0] = xi
    records["odp_cv_scaled"][:, 0] = xi
    records["odp_cv_perpendicular"][:, 0] = 0.0
    records["odp_perpendicular_norm"] = 0.0
    records["odp_bias_parallel_hartree"] = umbrella
    records["odp_bias_perpendicular_hartree"] = 0.0
    records["odp_bias_hartree"] = umbrella
    records["e_unbiased_pot_hartree"] = unbiased
    records["e_pot_hartree"] = unbiased + umbrella
    records["e_tot_hartree"] = unbiased + umbrella
    provenance = {
        "enabled": True, "window": window, "cv": ["distance(1,2)"],
        "cv_atom_indexing": "1-based", "cv_native_units": ["bohr"],
        "scale": [1.0],
        "reference_r": [0.0], "reference_p": [1.0],
        "scaled_path_length": 1.0,
        "center": center, "k_parallel_hartree": k_umbrella,
        "k_perpendicular_hartree": 0.0,
        "perpendicular_restraint": False,
        "projection": "signed_scaled_dot_over_path_norm_squared",
    }
    restart_identity = json.loads(json.dumps(system_identity))
    restart_identity["system"]["sha256"] = f"restart-window-{window}"
    header = {
        "schema_version": NAMD_TRAJECTORY_SCHEMA_VERSION,
        "nstate": 1, "natom": 1, "ncv": 1,
        "record_bytes": dtype.itemsize,
        "signature": json.dumps(restart_identity, sort_keys=True),
        "wham_system_identity": system_identity["system"],
        "ensemble": "NVT", "odp": provenance,
        "independent_controls": {
            "thermostat": {
                "ensemble": "nvt", "type": "langevin",
                "temperature_kelvin": temperature,
                "friction_ps_inverse": 1.0,
            },
        },
    }
    encoded = json.dumps(header, sort_keys=True).encode("utf-8")
    path = os.path.join(root, f"window-{window}.namd.trj")
    with open(path, "wb") as stream:
        stream.write(NAMD_TRAJECTORY_MAGIC)
        stream.write(struct.pack("<Q", len(encoded)))
        stream.write(encoded)
        stream.write(records.tobytes(order="C"))
    paths.append(path)

output = os.path.join(root, "wham.npz")
result = odp_wham(paths, temperature, bins=100, tolerance=1.0e-11, output=output)
with open(paths[0], "rb") as stream:
    segment_magic = stream.read(8)
    segment_header_size = struct.unpack("<Q", stream.read(8))[0]
    segment_header = stream.read(segment_header_size)
    segment_payload = stream.read()
segment_paths = []
segment_record_count = 6000
for name, payload in (
    ("segment-a.namd.trj", segment_payload[:segment_record_count*dtype.itemsize]),
    ("segment-b.namd.trj", segment_payload[segment_record_count*dtype.itemsize:]),
):
    segment_path = os.path.join(root, name)
    with open(segment_path, "wb") as stream:
        stream.write(segment_magic)
        stream.write(struct.pack("<Q", segment_header_size))
        stream.write(segment_header)
        stream.write(payload)
    segment_paths.append(segment_path)
segment_result = odp_wham(
    segment_paths, temperature, bins=100, discard=100, stride=3)
segment_steps = segment_result["sample_step"]
overlap_path = os.path.join(root, "copied-and-extended-window.namd.trj")
shutil.copyfile(paths[0], overlap_path)
extra_record = np.frombuffer(
    segment_payload[-dtype.itemsize:], dtype=dtype, count=1).copy()
extra_record["step"] = 12000
with open(overlap_path, "ab") as stream:
    stream.write(extra_record.tobytes(order="C"))
try:
    odp_wham([paths[0], overlap_path], temperature, bins=100)
except ValueError as error:
    overlapping_segments_rejected = "non-overlapping steps" in str(error)
else:
    overlapping_segments_rejected = False
try:
    odp_wham([paths[0]], temperature, bins=[-0.1, 0.0, 0.1])
except ValueError as error:
    incomplete_bins_rejected = "cover every selected" in str(error)
else:
    incomplete_bins_rejected = False
extension_input = os.path.join(root, "extension-window.npz")
shutil.copyfile(paths[0], extension_input)
try:
    odp_wham([extension_input], temperature, output=extension_input)
except ValueError as error:
    output_alias_rejected = "must not overwrite" in str(error)
else:
    output_alias_rejected = False
try:
    odp_wham(
        [extension_input], temperature,
        output=os.path.join(root, "extension-window"),
    )
except ValueError as error:
    implicit_extension_alias_rejected = "must not overwrite" in str(error)
else:
    implicit_extension_alias_rejected = False
hardlink_output = os.path.join(root, "hardlink-output.npz")
os.link(paths[0], hardlink_output)
try:
    odp_wham(paths, temperature, output=hardlink_output)
except ValueError as error:
    hardlink_output_rejected = "must not overwrite" in str(error)
else:
    hardlink_output_rejected = False
try:
    odp_wham(
        [paths[0], os.path.join(root, ".", os.path.basename(paths[0]))],
        temperature,
    )
except ValueError as error:
    duplicate_path_rejected = "duplicate trajectory input" in str(error)
else:
    duplicate_path_rejected = False
duplicate_copy = os.path.join(root, "duplicate-copy.namd.trj")
shutil.copyfile(paths[0], duplicate_copy)
try:
    odp_wham([paths[0], duplicate_copy], temperature)
except ValueError as error:
    duplicate_hash_rejected = "duplicate trajectory content" in str(error)
else:
    duplicate_hash_rejected = False
with open(paths[0], "rb") as stream:
    magic = stream.read(8)
    header_size = struct.unpack("<Q", stream.read(8))[0]
    other_header = json.loads(stream.read(header_size).decode("utf-8"))
    other_payload = stream.read()
other_header["wham_system_identity"]["sha256"] = "system-b"
other_encoded = json.dumps(other_header, sort_keys=True).encode("utf-8")
other_system = os.path.join(root, "other-system.namd.trj")
with open(other_system, "wb") as stream:
    stream.write(magic)
    stream.write(struct.pack("<Q", len(other_encoded)))
    stream.write(other_encoded)
    stream.write(other_payload)
try:
    odp_wham([paths[0], other_system], temperature)
except ValueError as error:
    other_system_rejected = "different molecular systems" in str(error)
else:
    other_system_rejected = False
restraint_header = json.loads(json.dumps(other_header))
restraint_header["wham_system_identity"] = system_identity["system"]
restraint_signature = json.loads(restraint_header["signature"])
restraint_signature["independent_controls"] = {
    "droplet": {
        "enabled": True, "center_angstrom": [0.0, 0.0, 0.0],
        "radius_angstrom": 12.0, "buffer_angstrom": 1.0,
        "force_constant_kcal_mol_angstrom2": 10.0,
        "target": "atoms", "atoms_zero_based": "0",
        "water_resnames": [], "max_penetration_angstrom": 5.0,
        "group_count": 1,
    },
    "solute_com": {"enabled": False},
}
restraint_header["signature"] = json.dumps(
    restraint_signature, sort_keys=True)
restraint_encoded = json.dumps(restraint_header, sort_keys=True).encode("utf-8")
other_restraint = os.path.join(root, "other-restraint.namd.trj")
with open(other_restraint, "wb") as stream:
    stream.write(magic)
    stream.write(struct.pack("<Q", len(restraint_encoded)))
    stream.write(restraint_encoded)
    stream.write(other_payload)
try:
    odp_wham([paths[0], other_restraint], temperature)
except ValueError as error:
    other_restraint_rejected = "different molecular systems" in str(error)
else:
    other_restraint_rejected = False
pcm_header = json.loads(json.dumps(other_header))
pcm_header["wham_system_identity"] = system_identity["system"]
pcm_signature = json.loads(pcm_header["signature"])
pcm_signature["pcm"]["epsilon"] = 40.0
pcm_header["signature"] = json.dumps(pcm_signature, sort_keys=True)
pcm_encoded = json.dumps(pcm_header, sort_keys=True).encode("utf-8")
other_pcm = os.path.join(root, "other-pcm.namd.trj")
with open(other_pcm, "wb") as stream:
    stream.write(magic)
    stream.write(struct.pack("<Q", len(pcm_encoded)))
    stream.write(pcm_encoded)
    stream.write(other_payload)
try:
    odp_wham([paths[0], other_pcm], temperature)
except ValueError as error:
    other_pcm_rejected = "different molecular systems" in str(error)
else:
    other_pcm_rejected = False
d4_header = json.loads(json.dumps(other_header))
d4_header["wham_system_identity"] = system_identity["system"]
d4_signature = json.loads(d4_header["signature"])
d4_signature["d4"] = True
d4_header["signature"] = json.dumps(d4_signature, sort_keys=True)
d4_encoded = json.dumps(d4_header, sort_keys=True).encode("utf-8")
other_d4 = os.path.join(root, "other-d4.namd.trj")
with open(other_d4, "wb") as stream:
    stream.write(magic)
    stream.write(struct.pack("<Q", len(d4_encoded)))
    stream.write(d4_encoded)
    stream.write(other_payload)
try:
    odp_wham([paths[0], other_d4], temperature)
except ValueError as error:
    other_d4_rejected = "different molecular systems" in str(error)
else:
    other_d4_rejected = False
race_path = os.path.join(root, "append-race.namd.trj")
shutil.copyfile(paths[0], race_path)
race_snapshot_bytes = os.path.getsize(race_path)
with open(race_path, "rb") as stream:
    stream.seek(-dtype.itemsize, os.SEEK_END)
    appended_record = stream.read(dtype.itemsize)
original_reader = namd_module.read_odp_wham_series
def read_then_append(path):
    series = original_reader(path)
    with open(path, "ab") as stream:
        stream.write(appended_record)
    return series
namd_module.read_odp_wham_series = read_then_append
try:
    race_result = odp_wham([race_path], temperature, bins=100)
finally:
    namd_module.read_odp_wham_series = original_reader
with open(race_path, "rb") as stream:
    race_expected_hash = hashlib.sha256(
        stream.read(race_snapshot_bytes)).hexdigest()
mean = np.sum(result["sample_weights"]*result["sample_xi"])
variance = np.sum(result["sample_weights"]*(result["sample_xi"] - mean)**2)
with np.load(output, allow_pickle=False) as saved:
    metadata = json.loads(str(saved["metadata_json"][0]))
    saved_centers = saved["window_centers"].tolist()
print("ODP_WHAM=" + json.dumps({
    "mean": mean,
    "variance": variance,
    "expected_variance": 1.0/(beta*k_unbiased),
    "iterations": result["iterations"],
    "residual": result["residual"],
    "probability_sum": float(result["probability"].sum()),
    "overlap_row_sums": result["window_overlap"].sum(axis=1).tolist(),
    "effective_sample_size": result["effective_sample_size"],
    "ensemble_warning": result["ensemble_warning"],
    "metadata_converged": metadata["converged"],
    "metadata_hashes": metadata["trajectory_sha256"],
    "metadata_snapshot_bytes": metadata["trajectory_snapshot_bytes"],
    "saved_centers": saved_centers,
    "duplicate_path_rejected": duplicate_path_rejected,
    "duplicate_hash_rejected": duplicate_hash_rejected,
    "output_alias_rejected": output_alias_rejected,
    "implicit_extension_alias_rejected": implicit_extension_alias_rejected,
    "hardlink_output_rejected": hardlink_output_rejected,
    "other_system_rejected": other_system_rejected,
    "other_restraint_rejected": other_restraint_rejected,
    "other_pcm_rejected": other_pcm_rejected,
    "other_d4_rejected": other_d4_rejected,
    "race_snapshot_bytes": race_result["trajectory_snapshot_bytes"],
    "race_initial_bytes": race_snapshot_bytes,
    "race_record_bytes": dtype.itemsize,
    "race_snapshot_hash": race_result["trajectory_sha256"],
    "race_expected_hash": race_expected_hash,
    "race_sample_count": int(race_result["sample_xi"].size),
    "race_final_bytes": os.path.getsize(race_path),
    "system_identity": result["system_identity"],
    "segment_sample_count": int(segment_steps.size),
    "segment_first_step": int(segment_steps[0]),
    "segment_first_step_after_boundary": int(
        segment_steps[np.searchsorted(segment_steps, segment_record_count)]),
    "incomplete_bins_rejected": incomplete_bins_rejected,
    "overlapping_segments_rejected": overlapping_segments_rejected,
}))
'''
    env = child_env()
    env["OQP_ODP_WHAM_ROOT"] = str(tmp_path)
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(
        (line for line in result.stdout.splitlines()
         if line.startswith("ODP_WHAM=")), None)
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("ODP_WHAM="))
    assert abs(values["mean"]) < 0.025
    assert values["variance"] == pytest.approx(
        values["expected_variance"], rel=0.04)
    assert values["residual"] < 1.0e-11
    assert values["probability_sum"] == pytest.approx(1.0, abs=1.0e-12)
    assert values["overlap_row_sums"] == pytest.approx([1.0]*5, abs=1.0e-12)
    assert values["effective_sample_size"] > 10000
    assert values["ensemble_warning"] is None
    assert values["metadata_converged"] is True
    example_output = tmp_path / "wham-example.npz"
    example_paths = sorted(tmp_path.glob("window-*.namd.trj"))
    example = subprocess.run(
        [
            sys.executable, str(ROOT / "examples" / "ODP" / "odp_wham.py"),
            "--temperature", "300", "--bins", "100",
            "--output", str(example_output),
            *(str(path) for path in example_paths),
        ],
        cwd=ROOT, env=env, capture_output=True, text=True, check=False,
    )
    assert example.returncode == 0, example.stdout + example.stderr
    assert "WHAM converged" in example.stdout
    with np.load(example_output, allow_pickle=False) as saved:
        example_metadata = json.loads(str(saved["metadata_json"][0]))
        assert example_metadata["converged"] is True
        assert len(example_metadata["trajectory_sha256"]) == 5
    assert all(len(value) == 64 for value in values["metadata_hashes"])
    assert all(value > 0 for value in values["metadata_snapshot_bytes"])
    assert values["saved_centers"] == pytest.approx(
        [-1.2, -0.6, 0.0, 0.6, 1.2])
    assert values["duplicate_path_rejected"] is True
    assert values["duplicate_hash_rejected"] is True
    assert values["output_alias_rejected"] is True
    assert values["implicit_extension_alias_rejected"] is True
    assert values["hardlink_output_rejected"] is True
    assert values["other_system_rejected"] is True
    assert values["other_restraint_rejected"] is True
    assert values["other_pcm_rejected"] is True
    assert values["other_d4_rejected"] is True
    assert values["race_snapshot_bytes"] == [values["race_initial_bytes"]]
    assert values["race_final_bytes"] == (
        values["race_initial_bytes"] + values["race_record_bytes"])
    assert values["race_snapshot_hash"] == [values["race_expected_hash"]]
    assert values["race_sample_count"] == 12000
    assert values["system_identity"]["system"]["sha256"] == "system-a"
    assert values["segment_sample_count"] == len(range(100, 12000, 3))
    assert values["segment_first_step"] == 100
    assert values["segment_first_step_after_boundary"] == 6001
    assert values["incomplete_bins_rejected"] is True
    assert values["overlapping_segments_rejected"] is True


def test_odp_restart_example_runs_generated_manifest(tmp_path):
    env = child_env()
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples" / "ODP" / "odp_restart.py"),
            "--workdir", str(tmp_path),
        ],
        cwd=ROOT, env=env, capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout + result.stderr)
    assert "ODP restart completed through step 2" in result.stdout
    assert (tmp_path / "odp-restart.namd.restart.npz").is_file()
    assert (tmp_path / "odp-restart.namd.trj").is_file()


def _run_openqp_deck(deck_path, workdir, threads="4"):
    env = child_env()
    env["OMP_NUM_THREADS"] = threads
    return subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", "--nompi", str(deck_path)],
        cwd=workdir, env=env, capture_output=True, text=True, check=False,
    )


def _trajectory_values(workdir, trajectory, extra=""):
    script = r"""
import json, sys, warnings
import numpy as np
from oqp.library.namd import read_namd_trajectory
from oqp.library.odp import odp_wham
header, records = read_namd_trajectory(sys.argv[1])
adjusted = (records["e_tot_hartree"]
            - records["thermostat_exchange_cumulative_hartree"])
values = {
    "ensemble": header["ensemble"],
    "thermostat": header["independent_controls"]["thermostat"],
    "window": header["odp"]["window"],
    "steps": records["step"].tolist(),
    "xi": records["odp_xi"].tolist(),
    "cv_raw": records["odp_cv_raw"].tolist(),
    "coordinates": records["coordinates_bohr"].tolist(),
    "bias_split": float(np.max(np.abs(
        records["e_unbiased_pot_hartree"] + records["odp_bias_hartree"]
        - records["e_pot_hartree"]))),
    "exchange": float(np.max(np.abs(
        records["thermostat_exchange_cumulative_hartree"]))),
    "adjusted_drift": float(np.max(np.abs(adjusted - adjusted[0]))),
}
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    result = odp_wham([sys.argv[1]], 300.0, bins=4)
values["wham_warning"] = result["ensemble_warning"]
values["wham_runtime_warnings"] = [
    str(item.message) for item in caught
    if issubclass(item.category, RuntimeWarning)]
try:
    odp_wham([sys.argv[1]], 310.0, bins=4)
    values["wham_310"] = None
except ValueError as exc:
    values["wham_310"] = str(exc)
print("TRJ=" + json.dumps(values))
"""
    result = subprocess.run(
        [sys.executable, "-c", script, trajectory], cwd=workdir,
        env=child_env(), capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        pytest.fail(result.stdout + result.stderr)
    marker = next(line for line in result.stdout.splitlines()
                  if line.startswith("TRJ="))
    return json.loads(marker.removeprefix("TRJ="))


def test_odp_langevin_nvt_namd_run_feeds_canonical_wham(tmp_path):
    name = "H2CO_BHHLYP-MRSF-NAMD-ODP-NVT"
    workdir = tmp_path / "ODP"
    workdir.mkdir()
    (tmp_path / "geometries").mkdir()
    shutil.copy(ROOT / "examples" / "ODP" / f"{name}.oqp", workdir)
    shutil.copy(ROOT / "examples" / "geometries" / "CH2O-2bc62dda4b8a.xyz",
                tmp_path / "geometries")
    result = _run_openqp_deck(workdir / f"{name}.oqp", workdir)
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout[-4000:] + result.stderr[-4000:])
    values = _trajectory_values(workdir, f"{name}.namd.trj")
    assert values["ensemble"] == "NVT"
    assert values["thermostat"]["type"] == "langevin"
    assert values["thermostat"]["temperature_kelvin"] == 300.0
    assert values["steps"] == [0, 1, 2, 3]
    assert values["bias_split"] < 1.0e-12
    # the thermostat really acted, and the bias stayed conservative: the
    # energy minus the cumulative thermostat exchange is conserved
    assert values["exchange"] > 1.0e-6
    assert values["adjusted_drift"] < 1.0e-5
    assert values["wham_warning"] is None
    assert values["wham_runtime_warnings"] == []
    assert "sampled at 300 K" in values["wham_310"]


PME_ODP_DECK = """\
mrsf(nstate=2)/bhhlyp/6-31g
namd(S0,scheme=Overlap,nstep=3,ensemble=nvt,thermostat=langevin,thermostat_friction=5.0,seed=11,restart_interval=3)
scf(conv=1e-09)
qmmm(forcefield_files="amber14-all.xml amber14/tip3p.xml",qm_atoms="8,9,16,17,18",cutoff=PME,ewald_tol=1e-06,lj_switch=true,h_lj=true,mm_charge_width=0.7)
odp(enabled=true,cv="distance(9,10); distance(9,17)",scale="1.0,1.0",reference_r="2.2,2.5",reference_p="2.5,2.8",center=0.4,k_parallel=0.05,k_perpendicular=0.01,window=2)
geom="ala_box.pdb"
"""


def test_odp_periodic_pme_qmmm_namd_run_measures_minimum_image_cvs(tmp_path):
    pytest.importorskip("openmm")
    shutil.copy(ROOT / "examples" / "QMMM" / "ala_box.pdb", tmp_path / "ala_box.pdb")
    deck = tmp_path / "pme.oqp"
    deck.write_text(PME_ODP_DECK)
    result = _run_openqp_deck(deck, tmp_path, threads="8")
    if result.returncode != 0:
        if "cannot load" in result.stderr or "No module named" in result.stderr:
            pytest.skip("compiled OpenQP runtime is not available")
        pytest.fail(result.stdout[-4000:] + result.stderr[-4000:])
    values = _trajectory_values(tmp_path, "pme.namd.trj")
    assert values["ensemble"] == "NVT"
    assert values["window"] == 2
    assert values["bias_split"] < 1.0e-12
    box = 16.0/0.52917721092
    for xyz, cvs in zip(values["coordinates"], values["cv_raw"]):
        xyz = np.asarray(xyz)
        for (a, b), cv in zip(((8, 9), (8, 16)), cvs):
            d = xyz[b] - xyz[a]
            d -= box*np.round(d/box)
            assert abs(cv - float(np.linalg.norm(d))) < 1.0e-10
    assert values["wham_warning"] is None
