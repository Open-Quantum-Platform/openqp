"""Native ODP umbrella equations, invariants, restart, and provenance gates."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


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
    sys.path.insert(0, str(ROOT / "pyoqp"))
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
    sys.path.insert(0, str(ROOT / "pyoqp"))
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
    env = os.environ.copy()
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
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


def test_odp_packed_trj_and_restart_preserve_wham_provenance(tmp_path):
    script = r'''
import json
import os
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
        "namd(S1,nstep=2,dt=0.5)\n"
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
    def get_state_tracking():
        return None
    def put_data(self, data):
        self.loaded = data

d = NAMD.__new__(NAMD)
d.mol = Mol(); d.nstate = 2; d.dt_fs = 0.5; d.dt_adaptive = False
d._t_fs = 0.0; d.active = 1; d.coef = np.array([1.+0j, 0.+0j])
d.vel = np.zeros((2, 3)); d.trajectory_interval = 1
d.trajectory_file = os.path.join(root, "odp.namd.trj")
d.nacme_audit_file = os.path.join(root, "odp.nacme.tsv")
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
d2.mol = Mol(); d2.nstate = 2; d2.dt_fs = 0.5
d2.seed = 1; d2.rng_stream = 0; d2.restart_requested = True
d2._restart_system_identity = d._restart_system_identity
d2._wham_system_identity = d._wham_system_identity
d2.restart_file = d.restart_file; d2.trajectory_file = d.trajectory_file
d2.nacme_audit_file = d.nacme_audit_file
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
    env = os.environ.copy()
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
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
try:
    odp_wham(paths, temperature, output=paths[0])
except ValueError as error:
    output_alias_rejected = "must not overwrite" in str(error)
else:
    output_alias_rejected = False
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
    "hardlink_output_rejected": hardlink_output_rejected,
    "other_system_rejected": other_system_rejected,
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
}))
'''
    env = os.environ.copy()
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
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
    assert values["hardlink_output_rejected"] is True
    assert values["other_system_rejected"] is True
    assert values["other_pcm_rejected"] is True
    assert values["other_d4_rejected"] is True
    assert values["race_snapshot_bytes"] == [values["race_initial_bytes"]]
    assert values["race_final_bytes"] == (
        values["race_initial_bytes"] + values["race_record_bytes"])
    assert values["race_snapshot_hash"] == [values["race_expected_hash"]]
    assert values["race_sample_count"] == 12000
    assert values["system_identity"]["system"]["sha256"] == "system-a"
