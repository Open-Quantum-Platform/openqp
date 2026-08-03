"""Finite-droplet, COM-restraint, and NVT native-kernel regression gates."""

import json
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_droplet_production_path_is_native_and_has_a_c_abi():
    source = (ROOT / "source/modules/namd.F90").read_text()
    header = (ROOT / "include/oqp.h").read_text()
    wrapper = (ROOT / "pyoqp/oqp/__init__.py").read_text()
    driver = (ROOT / "pyoqp/oqp/library/namd.py").read_text()

    assert 'bind(C, name="oqp_namd_droplet_boundary")' in source
    assert 'bind(C, name="oqp_namd_com_restraint")' in source
    assert 'bind(C, name="oqp_namd_langevin_thermostat")' in source
    assert "int oqp_namd_droplet_boundary(" in header
    assert "int oqp_namd_com_restraint(" in header
    assert "int oqp_namd_langevin_thermostat(" in header
    for name in (
        "oqp_namd_droplet_boundary",
        "oqp_namd_com_restraint",
        "oqp_namd_langevin_thermostat",
    ):
        assert repr(name) in wrapper
        assert f"oqp.{name}(" in driver
    assert "ODP" not in source


def test_built_droplet_force_matches_central_finite_difference_and_fails_safe():
    script = r'''
import json
import numpy as np
import oqp

def droplet(xyz, mass, groups, center=(0, 0, 0), radius=3.0, buffer=0.5,
            k=0.2, limit=20.0):
    xyz = np.ascontiguousarray(xyz, dtype=np.float64).reshape((-1, 3))
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    groups = np.ascontiguousarray(groups, dtype=np.int64)
    center = np.ascontiguousarray(center, dtype=np.float64)
    force = np.zeros_like(xyz); energy = np.zeros(1); penetration = np.zeros(1)
    active = np.zeros(1, dtype=np.int64)
    status = oqp.oqp_namd_droplet_boundary(
        len(mass), int(groups.max()),
        oqp.ffi.cast('double *', xyz.ctypes.data),
        oqp.ffi.cast('double *', mass.ctypes.data),
        oqp.ffi.cast('int64_t *', groups.ctypes.data),
        oqp.ffi.cast('double *', center.ctypes.data),
        radius, buffer, k, limit,
        oqp.ffi.cast('double *', energy.ctypes.data),
        oqp.ffi.cast('double *', force.ctypes.data),
        oqp.ffi.cast('double *', penetration.ctypes.data),
        oqp.ffi.cast('int64_t *', active.ctypes.data))
    return int(status), float(energy[0]), force, float(penetration[0]), int(active[0])

inside = droplet([[2.5, 0, 0]], [16.0], [1])
outside = droplet([[4.0, 0, 0]], [16.0], [1])

# Molecular-COM target: compare every Cartesian force component with -dU/dx.
xyz = np.array([[4.10, 0.20, -0.10], [3.80, -0.10, 0.15]])
mass = np.array([16.0, 1.0]); groups = np.array([1, 1], dtype=np.int64)
status, energy, force, penetration, active = droplet(xyz, mass, groups)
h = 1.0e-5
fd = np.zeros_like(xyz)
for atom in range(2):
    for axis in range(3):
        plus = xyz.copy(); minus = xyz.copy()
        plus[atom, axis] += h; minus[atom, axis] -= h
        ep = droplet(plus, mass, groups)[1]
        em = droplet(minus, mass, groups)[1]
        fd[atom, axis] = -(ep-em)/(2*h)

# Repeat inside the smoothing buffer rather than only in the harmonic tail.
smooth_xyz = np.array([[3.25, 0.10, 0.05], [3.10, -0.05, -0.02]])
smooth_force = droplet(smooth_xyz, mass, groups)[2]
smooth_fd = np.zeros_like(smooth_xyz)
for atom in range(2):
    for axis in range(3):
        plus = smooth_xyz.copy(); minus = smooth_xyz.copy()
        plus[atom, axis] += h; minus[atom, axis] -= h
        ep = droplet(plus, mass, groups)[1]
        em = droplet(minus, mass, groups)[1]
        smooth_fd[atom, axis] = -(ep-em)/(2*h)

# Continuity on both sides of the onset and smooth-buffer/harmonic join.
eps = 1.0e-7
onset_left = droplet([[3.0-eps, 0, 0]], [1.0], [1])
onset_right = droplet([[3.0+eps, 0, 0]], [1.0], [1])
join_left = droplet([[3.5-eps, 0, 0]], [1.0], [1])
join_right = droplet([[3.5+eps, 0, 0]], [1.0], [1])

# The penetration cap returns a recoverable status before huge-polynomial work.
huge = droplet([[1.0e100, 0, 0]], [1.0], [1], limit=10.0)

def com_restraint(xyz, mass, selected, center=(0.2, -0.1, 0.3), k=0.4):
    xyz = np.ascontiguousarray(xyz, dtype=np.float64).reshape((-1, 3))
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    selected = np.ascontiguousarray(selected, dtype=np.int64)
    center = np.ascontiguousarray(center, dtype=np.float64)
    force = np.zeros_like(xyz); energy = np.zeros(1); displacement = np.zeros(1)
    status = oqp.oqp_namd_com_restraint(
        len(mass), oqp.ffi.cast('double *', xyz.ctypes.data),
        oqp.ffi.cast('double *', mass.ctypes.data),
        oqp.ffi.cast('int64_t *', selected.ctypes.data),
        oqp.ffi.cast('double *', center.ctypes.data), k,
        oqp.ffi.cast('double *', energy.ctypes.data),
        oqp.ffi.cast('double *', force.ctypes.data),
        oqp.ffi.cast('double *', displacement.ctypes.data))
    return int(status), float(energy[0]), force, float(displacement[0])

selected = np.array([1, 1], dtype=np.int64)
cs, ce, cf, cd = com_restraint(xyz, mass, selected)
cfd = np.zeros_like(xyz)
for atom in range(2):
    for axis in range(3):
        plus = xyz.copy(); minus = xyz.copy()
        plus[atom, axis] += h; minus[atom, axis] -= h
        ep = com_restraint(plus, mass, selected)[1]
        em = com_restraint(minus, mass, selected)[1]
        cfd[atom, axis] = -(ep-em)/(2*h)

print('DROPLET=' + json.dumps({
    'inside': [inside[0], inside[1], inside[2].tolist(), inside[4]],
    'outside': [outside[0], outside[1], outside[2].tolist(), outside[4]],
    'fd_error': float(np.max(np.abs(force-fd))),
    'smooth_fd_error': float(np.max(np.abs(smooth_force-smooth_fd))),
    'com_fd_error': float(np.max(np.abs(cf-cfd))),
    'onset': [onset_left[1], onset_right[1],
              onset_left[2][0,0], onset_right[2][0,0]],
    'join_jump': [abs(join_left[1]-join_right[1]),
                  abs(join_left[2][0,0]-join_right[2][0,0])],
    'huge': [huge[0], huge[1], huge[3], huge[4]],
    'finite': bool(np.isfinite(force).all() and np.isfinite(cf).all()),
}))
'''
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        if ("cannot load" in result.stderr or "Cannot locate" in result.stderr
                or "undefined symbol" in result.stderr):
            pytest.skip("compiled OpenQP runtime with droplet ABI is not available")
        pytest.fail(result.stdout + result.stderr)
    marker = next(line for line in result.stdout.splitlines()
                  if line.startswith("DROPLET="))
    values = json.loads(marker.removeprefix("DROPLET="))
    assert values["inside"] == [0, 0.0, [[0.0, 0.0, 0.0]], 0]
    assert values["outside"][0] == 0
    assert values["outside"][1] > 0.0
    assert values["outside"][2][0][0] < 0.0
    assert values["outside"][3] == 1
    assert values["fd_error"] < 2.0e-9
    assert values["smooth_fd_error"] < 2.0e-9
    assert values["com_fd_error"] < 2.0e-9
    assert max(abs(item) for item in values["onset"]) < 1.0e-12
    assert values["join_jump"][0] < 3.0e-8
    assert values["join_jump"][1] < 3.0e-7
    assert values["huge"][0] == 1
    assert values["huge"][1] == 0.0
    assert values["huge"][2] > 1.0e90
    assert values["huge"][3] == 1
    assert values["finite"]


def test_built_langevin_step_is_reproducible_and_reports_exact_heat():
    script = r'''
import json
import numpy as np
import oqp

mass = np.ascontiguousarray([1837.0, 29156.0], dtype=np.float64)
initial = np.ascontiguousarray([[0.001, -0.002, 0.0005],
                                [0.0002, 0.0001, -0.0003]], dtype=np.float64)

def run(step):
    vel = initial.copy(); heat = np.zeros(1)
    before = 0.5*np.sum(mass[:,None]*vel**2)
    status = oqp.oqp_namd_langevin_thermostat(
        len(mass), 20.0, 300.0, 1.0e-4, 20260803, 7, step,
        oqp.ffi.cast('double *', mass.ctypes.data),
        oqp.ffi.cast('double *', vel.ctypes.data),
        oqp.ffi.cast('double *', heat.ctypes.data))
    after = 0.5*np.sum(mass[:,None]*vel**2)
    return int(status), vel, float(heat[0]), float(after-before)

a = run(9); b = run(9); c = run(10)
print('THERMOSTAT=' + json.dumps({
    'status': [a[0], b[0], c[0]],
    'same': bool(np.array_equal(a[1], b[1])),
    'different': bool(not np.array_equal(a[1], c[1])),
    'heat_error': abs(a[2]-a[3]),
    'finite': bool(np.isfinite(a[1]).all() and np.isfinite(a[2])),
}))
'''
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        pytest.skip("compiled OpenQP runtime with Langevin ABI is not available")
    marker = next(line for line in result.stdout.splitlines()
                  if line.startswith("THERMOSTAT="))
    values = json.loads(marker.removeprefix("THERMOSTAT="))
    assert values == {
        "status": [0, 0, 0], "same": True, "different": True,
        "heat_error": pytest.approx(0.0, abs=2.0e-15), "finite": True,
    }


def test_small_water_droplet_nve_nvt_and_restart_smoke():
    """Propagate three rigid water COM groups using the production kernels."""
    script = r'''
import json
import numpy as np
import oqp

mass = np.ascontiguousarray([29156.0, 1837.0, 1837.0]*3, dtype=np.float64)
groups = np.ascontiguousarray([1,1,1,2,2,2,3,3,3], dtype=np.int64)
center = np.zeros(3, dtype=np.float64)
base = np.array([[0.0,0.0,0.0], [0.12,0.0,0.0], [-0.03,0.11,0.0]])
xyz0 = np.vstack((base+[3.4,0.0,0.0],
                  base+[-3.3,0.2,0.0],
                  base+[0.0,3.5,-0.1])).astype(np.float64)
vel0 = np.zeros_like(xyz0)
vel0[0:3,1] = 1.0e-4; vel0[3:6,2] = -8.0e-5; vel0[6:9,0] = 9.0e-5
dt = 0.1

def boundary(xyz):
    xyz = np.ascontiguousarray(xyz, dtype=np.float64)
    force = np.zeros_like(xyz); energy=np.zeros(1); penetration=np.zeros(1)
    active=np.zeros(1,dtype=np.int64)
    status=oqp.oqp_namd_droplet_boundary(
        len(mass), 3, oqp.ffi.cast('double *',xyz.ctypes.data),
        oqp.ffi.cast('double *',mass.ctypes.data),
        oqp.ffi.cast('int64_t *',groups.ctypes.data),
        oqp.ffi.cast('double *',center.ctypes.data), 3.0,0.4,0.08,5.0,
        oqp.ffi.cast('double *',energy.ctypes.data),
        oqp.ffi.cast('double *',force.ctypes.data),
        oqp.ffi.cast('double *',penetration.ctypes.data),
        oqp.ffi.cast('int64_t *',active.ctypes.data))
    if status != 0: raise RuntimeError(status)
    return force,float(energy[0]),float(penetration[0]),int(active[0])

def thermostat(vel,step):
    vel=np.ascontiguousarray(vel,dtype=np.float64); heat=np.zeros(1)
    status=oqp.oqp_namd_langevin_thermostat(
        len(mass),dt,300.0,2.0e-4,77,4,step,
        oqp.ffi.cast('double *',mass.ctypes.data),
        oqp.ffi.cast('double *',vel.ctypes.data),
        oqp.ffi.cast('double *',heat.ctypes.data))
    if status != 0: raise RuntimeError(status)
    return vel,float(heat[0])

def integrate(nstep,nvt=False,state=None,start=0):
    if state is None:
        xyz=xyz0.copy(); vel=vel0.copy(); force,epot,pmax,nactive=boundary(xyz)
    else:
        xyz,vel,force,epot,pmax,nactive,qcum=state
        xyz=xyz.copy(); vel=vel.copy(); force=force.copy()
    qcum = 0.0 if state is None else float(state[6])
    e0=0.5*np.sum(mass[:,None]*vel**2)+epot-qcum
    max_drift=0.0
    for step in range(start+1,start+nstep+1):
        accel=force/mass[:,None]
        xyz=xyz+vel*dt+0.5*accel*dt*dt
        force_new,epot,pmax,nactive=boundary(xyz)
        vel=vel+0.5*(force+force_new)/mass[:,None]*dt
        force=force_new
        if nvt:
            vel,heat=thermostat(vel,step); qcum+=heat
        accounted=0.5*np.sum(mass[:,None]*vel**2)+epot-qcum
        max_drift=max(max_drift,abs(accounted-e0))
    return (xyz,vel,force,epot,pmax,nactive,qcum),max_drift

nve,nve_drift=integrate(400,nvt=False)
nvt,nvt_drift=integrate(40,nvt=True)
half,_=integrate(20,nvt=True)
split,split_drift=integrate(20,nvt=True,state=half,start=20)
print('SMOKE='+json.dumps({
    'nve_drift':nve_drift,'nvt_adjusted_drift':nvt_drift,
    'nvt_heat':nvt[6],'active':nve[5],'penetration':nve[4],
    'restart_xyz':bool(np.array_equal(nvt[0],split[0])),
    'restart_vel':bool(np.array_equal(nvt[1],split[1])),
    'restart_heat':bool(nvt[6]==split[6]),
    'finite':bool(all(np.isfinite(x).all() if hasattr(x,'shape') else np.isfinite(x)
                      for x in nvt)),
}))
'''
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=ROOT, env=env,
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        pytest.skip("compiled OpenQP runtime with droplet ABI is not available")
    marker = next(line for line in result.stdout.splitlines()
                  if line.startswith("SMOKE="))
    values = json.loads(marker.removeprefix("SMOKE="))
    assert values["nve_drift"] < 2.0e-12
    assert values["nvt_adjusted_drift"] < 2.0e-10
    assert abs(values["nvt_heat"]) > 1.0e-12
    assert values["active"] >= 1
    assert values["penetration"] >= 0.0
    assert values["restart_xyz"]
    assert values["restart_vel"]
    assert values["restart_heat"]
    assert values["finite"]


@pytest.mark.skipif(
    importlib.util.find_spec("openmm") is None,
    reason="optional OpenMM backend is not installed",
)
def test_openmm_five_water_droplet_nve_nvt_smoke(tmp_path):
    """Exercise the complete QM/MM path with five real TIP3P waters."""
    example = ROOT / "examples/QMMM"
    geometry = ROOT / "examples/geometries/CH2O-2bc62dda4b8a.xyz"
    common = "\n".join((
        "mrsf(nstate=2)/bhhlyp/6-31g*",
        "{driver}",
        "guess(type=huckel)",
        (
            f'qmmm(pdb_file="{example / "formaldehyde_water.pdb"}",'
            f'forcefield_files="{example / "formaldehyde.xml"} '
            f'{example / "tip3p.xml"}",qm_atoms="0-3",cutoff=NoCutoff,'
            'embedding=electrostatic)'
        ),
        (
            'droplet(enabled=true,center="0,0,0",radius=3.3,buffer=1.0,'
            'force_constant=10.0,target=water_com,max_penetration=3.0)'
        ),
        ('solute_com(enabled=true,center="0,0,0",force_constant=2.0,'
         'atoms="0-3")'),
        f'geom="{geometry}"',
        "",
    ))
    drivers = {
        "nve": (
            "namd(T0,nstep=1,dt=0.5,substep=50,init_temp=300,"
            "velocity=maxwell,seed=8128,rng_stream=7,first_hop_step=1,"
            "nacme_check=off,nacme_gate=off,nve_gate=warn,ensemble=nve,"
            "thermostat=off,trajectory_interval=1,restart_interval=1)"
        ),
        "nvt": (
            "namd(T0,nstep=1,dt=0.5,substep=50,init_temp=300,"
            "velocity=maxwell,seed=8128,rng_stream=7,first_hop_step=1,"
            "nacme_check=off,nacme_gate=off,nve_gate=off,ensemble=nvt,"
            "thermostat=langevin,thermostat_temperature=300,"
            "thermostat_friction=1,trajectory_interval=1,restart_interval=1)"
        ),
    }
    env = os.environ.copy()
    env["OPENQP_ROOT"] = str(ROOT)
    env["OMP_NUM_THREADS"] = "2"
    env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get(
        "PYTHONPATH", "")

    for ensemble, driver in drivers.items():
        input_file = tmp_path / f"{ensemble}.oqp"
        input_file.write_text(common.format(driver=driver), encoding="utf-8")
        result = subprocess.run(
            [sys.executable, "-m", "oqp.pyoqp", "--silent", "--nompi",
             str(input_file)],
            cwd=tmp_path, env=env, capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            if ("cannot load" in result.stderr or "Cannot locate" in result.stderr
                    or "undefined symbol" in result.stderr):
                pytest.skip("compiled OpenQP runtime is not available")
            pytest.fail(result.stdout + result.stderr)

        inspection = subprocess.run(
            [sys.executable, "-c", (
                "import json,sys\n"
                "from oqp.library.namd import read_namd_trajectory\n"
                "header,records=read_namd_trajectory(sys.argv[1])\n"
                "last=records[-1]\n"
                "print(json.dumps({"
                "'groups':header['independent_controls']['droplet']['group_count'],"
                "'length':len(records),"
                "'active':int(last['droplet_active_count']),"
                "'energy':float(last['droplet_energy_hartree']),"
                "'nve_verdict':int(last['nve_verdict'])}))\n"
            ), str(tmp_path / f"{ensemble}.namd.trj")],
            cwd=tmp_path, env=env, capture_output=True, text=True, check=False,
        )
        if inspection.returncode != 0:
            pytest.fail(inspection.stdout + inspection.stderr)
        trajectory = json.loads(inspection.stdout)
        with np.load(
                tmp_path / f"{ensemble}.namd.restart.npz",
                allow_pickle=False) as checkpoint:
            assert float(checkpoint["droplet_energy"][0]) > 0.0
            assert float(checkpoint["droplet_max_penetration"][0]) > 0.0
            assert int(checkpoint["droplet_active_count"][0]) == 5
            thermostat_exchange = float(
                checkpoint["thermostat_exchange_cumulative"][0])
        assert trajectory["groups"] == 5
        assert trajectory["length"] == 2
        assert trajectory["active"] == 5
        assert trajectory["energy"] > 0.0
        log = (tmp_path / f"{ensemble}.log").read_text(encoding="utf-8")
        assert "drop_n=5" in log
        assert "droplet(enabled=true" in (
            tmp_path / f"{ensemble}.namd.restart.oqp").read_text(
                encoding="utf-8")
        if ensemble == "nve":
            assert trajectory["nve_verdict"] == 1
            assert thermostat_exchange == 0.0
        else:
            assert trajectory["nve_verdict"] == -1
            assert abs(thermostat_exchange) > 1.0e-12

    assert (tmp_path / "nve.namd.restart.oqp").exists()
    assert (tmp_path / "nvt.namd.restart.oqp").exists()
