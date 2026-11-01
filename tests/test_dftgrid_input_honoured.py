"""An explicit [dftgrid] rad_npts / ang_npts request must reach the grid builder.

Until 2026-10 the ``pruned`` option defaulted to a fixed preset (SG2, before
that SG1).  A preset defines its own radial and angular grids, so
``dft_set_options`` (source/dftlib/dft.F90) never read ``rad_npts`` /
``ang_npts`` for the atoms it covers.  An input such as
``dftgrid(rad_npts=90,ang_npts=302)`` -- or the legacy ``[dftgrid]`` section
with the same two keys -- therefore ran on SG2 (75 radial points), for every
runtype, with the requested values still echoed in the printed configuration.
Measured on H2O BHHLYP/6-31G*: the "90/302" energy was bit-identical to SG2 and
9.5e-6 Eh away from the real unpruned 90/302 grid.

``pruned`` now defaults to ``auto`` (oqp.molecule.oqpdata.resolve_dftgrid_pruned):
SG2 while both sizes keep their defaults, the unpruned rad_npts x ang_npts grid
as soon as either is changed.  These tests run the real input routes (concise
.oqp, legacy INI, energy and MD runtypes) through the installed ``openqp``
entry point and read the grid that the native builder announces.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from oqp.molecule.oqpdata import resolve_dftgrid_pruned
from oqp.utils.input_checker import check_input_values


H2O_XYZ = """3
water
O 0.0 0.0 -0.041061554
H -0.533194329 0.533194329 -0.614469223
H 0.533194329 -0.533194329 -0.614469223
"""

H2O_INLINE = """
 8 0.0 0.0 -0.041061554
 1 -0.533194329 0.533194329 -0.614469223
 1 0.533194329 -0.533194329 -0.614469223"""

METHOD = "rks/bhhlyp/6-31g* scf(conv=1e-9)"

LEGACY = """[input]
system={inline}
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=hf

[scf]
type=rhf
conv=1e-9

[dftgrid]
{grid}
"""

SG2_RE = re.compile(r"Standard Grid 2 \(SG2\).*?NRAD\s*=\s*(\d+)", re.S)
LEBEDEV_RE = re.compile(
    r"Lebedev grid-based DFT options.*?NRAD\s*=\s*(\d+)\s+NLEB\s*=\s*(\d+)", re.S)
ENERGY_RE = re.compile(r"Final RHF energy is\s+(-?\d+\.\d+)")


# --------------------------------------------------------------------------
# Resolution rule
# --------------------------------------------------------------------------

@pytest.mark.parametrize("section, expected", [
    ({}, "SG2"),
    ({"pruned": "auto"}, "SG2"),
    ({"pruned": "auto", "rad_npts": 96, "ang_npts": 302}, "SG2"),
    ({"pruned": "auto", "rad_npts": 90}, "none"),
    ({"pruned": "auto", "ang_npts": 590}, "none"),
    ({"rad_npts": "90", "ang_npts": "302"}, "none"),        # untyped strings
    ({"pruned": "SG2", "rad_npts": 90}, "SG2"),             # explicit preset wins
    ({"pruned": "sg1", "rad_npts": 200}, "SG1"),
    ({"pruned": "SG3"}, "SG3"),
    ({"pruned": "", "rad_npts": 96}, "none"),
    ({"pruned": "None"}, "none"),
    ({"pruned": "off"}, "none"),
])
def test_resolution_rule(section, expected):
    assert resolve_dftgrid_pruned(section) == expected


def _checker_messages(dftgrid):
    config = {"input": {"method": "hf", "functional": "bhhlyp",
                        "system": H2O_INLINE, "runtype": "energy"},
              "dftgrid": dftgrid}
    report = check_input_values(config, raise_error=False, emit=False)
    return [(d.severity, d.message) for d in report.diagnostics
            if d.path == "dftgrid.pruned"]


def test_checker_rejects_unknown_preset():
    messages = _checker_messages({"pruned": "SG9"})
    assert any(sev == "ERROR" for sev, _ in messages), messages


def test_checker_warns_when_explicit_preset_ignores_sizes():
    messages = _checker_messages({"pruned": "SG2", "rad_npts": 90, "ang_npts": 302})
    assert any(sev == "WARNING" and "rad_npts=90" in msg for sev, msg in messages), messages
    assert _checker_messages({"pruned": "auto", "rad_npts": 90}) == []
    assert _checker_messages({"pruned": "none", "rad_npts": 90}) == []


# --------------------------------------------------------------------------
# Native state after configuration
# --------------------------------------------------------------------------

def _configured_dft(tmp_path, grid_lines):
    from oqp.molecule import Molecule

    inp = tmp_path / "state.inp"
    inp.write_text(LEGACY.format(inline=H2O_INLINE, grid=grid_lines), encoding="utf-8")
    mol = Molecule("state", str(inp), str(tmp_path / "state.log"), silent=1)
    mol.load_config(str(inp))
    return mol


def _native_grid(mol):
    dft = mol.data._data.dft
    name = b"".join(dft.grid_pruned_name).decode("ascii").strip("\x00 ")
    return bool(dft.grid_pruned), name, int(dft.grid_rad_size), int(dft.grid_ang_size)


def test_native_state_follows_the_request(tmp_path):
    mol = _configured_dft(tmp_path, "rad_npts=90\nang_npts=302")
    assert _native_grid(mol) == (False, "", 90, 302)
    assert mol.config["dftgrid"]["pruned"] == "none"

    mol = _configured_dft(tmp_path, "")
    assert _native_grid(mol)[:2] == (True, "SG2")
    assert mol.config["dftgrid"]["pruned"] == "SG2"


def test_reconfiguring_to_unpruned_clears_an_earlier_preset(tmp_path):
    """MD/NAMD/API set() re-apply a config to a live Molecule; no stale preset."""
    mol = _configured_dft(tmp_path, "pruned=SG2")
    assert _native_grid(mol)[:2] == (True, "SG2")
    inp = tmp_path / "state2.inp"
    inp.write_text(LEGACY.format(inline=H2O_INLINE, grid="pruned=none\nrad_npts=90"),
                   encoding="utf-8")
    mol.load_config(str(inp))
    assert _native_grid(mol) == (False, "", 90, 302)


# --------------------------------------------------------------------------
# End to end through the openqp entry point
# --------------------------------------------------------------------------

def _run_openqp(tmp_path, name, text):
    workdir = tmp_path / name
    workdir.mkdir()
    (workdir / "h2o.xyz").write_text(H2O_XYZ, encoding="utf-8")
    suffix = ".inp" if text.lstrip().startswith("[") else ".oqp"
    deck = workdir / f"{name}{suffix}"
    deck.write_text(text, encoding="utf-8")
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "2")
    proc = subprocess.run(
        [sys.executable, "-c", "import sys; from oqp.pyoqp import main; sys.exit(main())",
         deck.name],
        cwd=workdir, env=env, capture_output=True, text=True, timeout=900)
    log = workdir / f"{name}.log"
    assert proc.returncode == 0, (proc.stdout[-3000:], proc.stderr[-3000:])
    assert log.exists(), proc.stdout[-3000:]
    return log.read_text(encoding="utf-8", errors="replace")


def _grid_and_energy(log):
    lebedev = LEBEDEV_RE.search(log)
    sg2 = SG2_RE.search(log)
    energies = ENERGY_RE.findall(log)
    assert energies, "no SCF energy in log"
    if lebedev:
        grid = ("unpruned", int(lebedev.group(1)), int(lebedev.group(2)))
    elif sg2:
        grid = ("SG2", int(sg2.group(1)))
    else:
        raise AssertionError("no grid announcement in log")
    return grid, energies


@pytest.fixture(scope="module")
def reference_runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("dftgrid_ref")
    geom = 'geom="h2o.xyz"'
    return {
        "none_90_302": _grid_and_energy(_run_openqp(
            tmp, "none", f"{METHOD} dftgrid(pruned=none,rad_npts=90,ang_npts=302)\n{geom}\n")),
        "default": _grid_and_energy(_run_openqp(tmp, "default", f"{METHOD}\n{geom}\n")),
    }


def test_default_grid_is_still_sg2(reference_runs):
    grid, _ = reference_runs["default"]
    assert grid == ("SG2", 75)


def test_concise_energy_runs_the_requested_grid(tmp_path, reference_runs):
    log = _run_openqp(tmp_path, "concise",
                      f'{METHOD} dftgrid(rad_npts=90,ang_npts=302)\ngeom="h2o.xyz"\n')
    grid, energies = _grid_and_energy(log)
    ref_grid, ref_energies = reference_runs["none_90_302"]
    assert grid == ("unpruned", 90, 302) == ref_grid
    assert energies[-1] == ref_energies[-1]
    # and it is not the preset: the two grids differ by ~1e-5 Eh here
    assert abs(float(energies[-1]) - float(reference_runs["default"][1][-1])) > 1e-6


def test_legacy_section_runs_the_requested_grid(tmp_path, reference_runs):
    log = _run_openqp(tmp_path, "legacy",
                      LEGACY.format(inline=H2O_INLINE, grid="rad_npts=90\nang_npts=302"))
    grid, energies = _grid_and_energy(log)
    assert grid == ("unpruned", 90, 302)
    assert energies[-1] == reference_runs["none_90_302"][1][-1]


def test_md_runs_the_requested_grid(tmp_path, reference_runs):
    log = _run_openqp(
        tmp_path, "md",
        f"{METHOD} md(nstep=1,dt=0.5,velocity=zero,seed=1) "
        f'dftgrid(rad_npts=90,ang_npts=302)\ngeom="h2o.xyz"\n')
    grid, energies = _grid_and_energy(log)
    assert grid == ("unpruned", 90, 302)
    assert "Standard Grid" not in log
    # step 0 is the same geometry as the single point
    assert energies[0] == reference_runs["none_90_302"][1][-1]
