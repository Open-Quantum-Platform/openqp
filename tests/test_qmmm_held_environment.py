"""Held-atom regressions without an electronic-structure backend.

Execute the production Python methods directly; OpenMM system construction is
optional, and the WHAM identity tests need only NumPy. No molecular calculation
or replacement implementation of the selection/identity algorithms is used.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / 'pyoqp' / 'oqp' / 'library'


def _production_class(filename, class_name, methods, namespace, functions=(), constants=()):
    tree = ast.parse((LIBRARY / filename).read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name)
    body = [n for n in cls.body if
            (isinstance(n, ast.FunctionDef) and n.name in methods) or
            (isinstance(n, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id in constants for t in n.targets))]
    selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in functions]
    selected.append(ast.ClassDef(name=class_name, bases=[], keywords=[], body=body, decorator_list=[]))
    module = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    exec(compile(module, str(LIBRARY / filename), 'exec'), namespace)
    return namespace[class_name]


@pytest.fixture
def wham_driver():
    cls = _production_class('namd.py', 'NAMD_QMMM',
        {'_qmmm_wham_system_identity', '_qmmm_identity_config'},
        {'np': np, 'json': json, 'hashlib': hashlib, 'struct': struct},
        functions={'_restart_identity_digest'},
        constants={'_QMMM_OPTIONAL_HAMILTONIAN_KEYS', '_QMMM_SELECTION_KEYS'})
    driver = cls()
    residue = SimpleNamespace(name='MOL', index=0, chain=SimpleNamespace(id='A'))
    atoms = [SimpleNamespace(name=f'H{i}', element=SimpleNamespace(atomic_number=1),
                             residue=residue, index=i) for i in range(4)]
    driver.pdb = SimpleNamespace(topology=SimpleNamespace(
        atoms=lambda: iter(atoms), bonds=lambda: iter([(atoms[0], atoms[1])])))
    driver._mm = SimpleNamespace(XmlSerializer=SimpleNamespace(serialize=lambda _: '<System/>'))
    driver.natom_all = 4
    driver.qm_atoms = np.array([0, 1])
    driver.m_all = np.ones(4)
    driver.r_all = np.arange(12, dtype=float).reshape(4, 3)
    driver._held_atoms = {3}
    return driver


def _identity(driver, **selection):
    return driver._qmmm_wham_system_identity(None, {'embedding': 'electrostatic', **selection})


def test_wham_distinguishes_resolved_held_sets_with_identical_input(wham_driver):
    before = _identity(wham_driver, active_from_pdb=True)
    wham_driver._held_atoms = {2}
    assert _identity(wham_driver, active_from_pdb=True) != before


def test_wham_distinguishes_fixed_coordinates(wham_driver):
    before = _identity(wham_driver)
    wham_driver.r_all[3, 0] += 0.5
    assert _identity(wham_driver) != before


def test_wham_allows_different_mobile_coordinates_and_equivalent_selections(wham_driver):
    before = _identity(wham_driver, active_atoms='0-2')
    wham_driver.r_all[:3] += 1.0
    assert _identity(wham_driver, frozen_atoms='3') == before
    # Repeated calls own no mutable coordinate views or stale digest cache.
    wham_driver.r_all[3, 1] += 0.5
    assert _identity(wham_driver) != before
    wham_driver.r_all[3, 1] -= 0.5
    assert _identity(wham_driver) == before


def test_wham_unconstrained_identity_is_backward_compatible(wham_driver):
    wham_driver._held_atoms = set()
    before = _identity(wham_driver)
    assert before['sha256'] == '1bd49fbc7a4f62bc0dd1bb145f57a40d6f530cc522780b353dbb9dece8ae4558'
    wham_driver.r_all += 2.0
    assert _identity(wham_driver, active_atoms='0-3') == before
    # Legacy callers with no selection state retain the all-moving identity.
    del wham_driver._held_atoms
    assert _identity(wham_driver) == before


@pytest.fixture
def md_driver():
    mm = pytest.importorskip('openmm')
    from openmm import app, unit
    spec = importlib.util.spec_from_file_location('selection_under_test', LIBRARY / 'qmmm_active.py')
    selection = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(selection)
    namespace = {'np': np, 'mm': mm, 'app': app, 'unit': unit}
    for name in ('selection_requested', 'resolve_active_set', 'freeze_constrained_partners', 'held_atoms'):
        namespace[name] = getattr(selection, name)
    cls = _production_class('qmmm_md.py', 'QMMM_MD',
        {'_resolve_frozen_atoms', '_build_md_system'}, namespace,
        functions={'_copy_virtual_site', '_rigid_water_constraints', '_to_kJmol'})
    driver = cls()
    pdb_path = ROOT / 'examples' / 'QMMM' / 'formaldehyde_water_active.pdb'
    driver.pdb = app.PDBFile(str(pdb_path))
    driver._pdb_path = str(pdb_path)
    driver.qm_atoms = np.array([0, 1, 2, 3])
    driver.rigidwater = False
    driver.ensemble = 'npt'
    driver._selection_cfg = {}
    driver.system_md = None
    driver.pressure = 1.0 * unit.bar
    driver.temperature = 300.0 * unit.kelvin
    driver.barostat_interval = 25
    count = driver.pdb.topology.getNumAtoms()
    system = mm.System()
    for _ in range(count):
        system.addParticle(1.0)
    driver.mm_systems = {'sys0': system}
    driver.oqp_driver = SimpleNamespace(
        _box_lengths_bohr=lambda: None,
        compute_force=Mock(return_value=(0.0, np.zeros((count, 3)))))
    return driver


@pytest.mark.parametrize('selection', [
    {'frozen_atoms': '7-9'}, {'active_atoms': '0-6'},
    {'active_from_pdb': True}, {'active_radius': 0.1},
])
def test_npt_rejects_resolved_held_atoms_before_force_evaluation(md_driver, selection):
    md_driver._selection_cfg = selection
    with pytest.raises(ValueError, match='NPT.*held atoms'):
        md_driver._build_md_system()
    md_driver.oqp_driver.compute_force.assert_not_called()
    assert md_driver.system_md is None
    # A rejected initialization can be retried with a supported selection.
    md_driver._selection_cfg = {}
    md_driver._build_md_system()
    assert not md_driver.frozen_atoms


@pytest.mark.parametrize('ensemble', ['nve', 'nvt'])
def test_fixed_volume_ensembles_still_hold_atoms(md_driver, ensemble):
    from openmm import unit
    md_driver.ensemble = ensemble
    md_driver._selection_cfg = {'active_atoms': '0-6'}
    md_driver._build_md_system()
    assert md_driver.frozen_atoms == set(range(7, 19))
    for i in range(19):
        mass = md_driver.system_md.getParticleMass(i).value_in_unit(unit.dalton)
        assert mass == (1.0 if i < 7 else 0.0)


@pytest.mark.parametrize('selection', [{}, {'active_atoms': '0-18'}])
def test_npt_without_held_atoms_keeps_the_barostat(md_driver, selection):
    import openmm as mm
    md_driver._selection_cfg = selection
    md_driver._build_md_system()
    assert not md_driver.frozen_atoms
    assert any(isinstance(force, mm.MonteCarloBarostat) for force in md_driver.system_md.getForces())
