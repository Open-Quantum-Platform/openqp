"""Regression tests for the periodic/embedding QM/MM controls added with the
link-atom + Ewald work ([qmmm] ewald_tol, lj_switch, h_lj, mm_charge_width) and
for the review fixes around them:

* lj_switch / h_lj are boolean schema keys (checked from the schema source, no
  runtime needed);
* the Ewald branch is entered only for genuinely periodic OpenMM nonbonded
  methods (PME / Ewald / LJPME / CutoffPeriodic), never for NoCutoff or
  CutoffNonPeriodic even when the topology carries box vectors;
* mm_charge_width / ewald_tol must be finite and positive;
* h_lj assigns Lennard-Jones parameters to MM hydrogens only, never to QM atoms;
* the split-embedding NAMD path folds link-atom charges onto their QM hosts
  (total QM charge conserved) before the MM electrostatics;
* the WHAM/restart system identity covers the new Hamiltonian options.

Everything below the schema test needs OpenMM (and the oqp package importable).
"""
import ast
import importlib.util
import os
import types
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_EXAMPLES = os.path.join(_ROOT, 'examples', 'QMMM')
_SCHEMA = os.path.join(_ROOT, 'pyoqp', 'oqp', 'molecule', 'oqpdata.py')

_HAVE_OPENMM = importlib.util.find_spec('openmm') is not None
try:
    from oqp.library.qmmm_driver import OpenQpQMMM
    from oqp.library.namd import NAMD_QMMM
    _HAVE_OQP = True
except Exception:  # pragma: no cover - no compiled runtime / no OpenMM
    _HAVE_OQP = False


def _schema_qmmm_types():
    tree = ast.parse(open(_SCHEMA).read())
    node = next(n.value for n in ast.walk(tree)
                if isinstance(n, ast.Assign)
                and any(getattr(t, 'id', '') == 'OQP_CONFIG_SCHEMA' for t in n.targets))
    for sk, sv in zip(node.keys, node.values):
        if ast.literal_eval(sk) != 'qmmm':
            continue
        out = {}
        for ok, ov in zip(sv.keys, sv.values):
            spec = {ast.literal_eval(k): v for k, v in zip(ov.keys, ov.values)}
            out[ast.literal_eval(ok)] = (spec['type'].id, ast.literal_eval(spec['default']))
        return out
    raise AssertionError('no [qmmm] section in the schema')


class TestSchema(unittest.TestCase):
    def test_periodic_options_types(self):
        t = _schema_qmmm_types()
        self.assertEqual(t['lj_switch'], ('bool', 'False'))
        self.assertEqual(t['h_lj'], ('bool', 'False'))
        self.assertEqual(t['ewald_tol'][0], 'str')
        self.assertEqual(t['mm_charge_width'][0], 'str')

    def test_examples_exercise_the_options(self):
        deck = open(os.path.join(_EXAMPLES, 'ala-box_BHHLYP-MRSF-NAMD-QMMM-PME.inp')).read()
        for key in ('ewald_tol', 'lj_switch', 'h_lj', 'mm_charge_width'):
            self.assertRegex(deck, r'(?m)^\s*%s\s*=\s*\S' % key)
        self.assertRegex(deck, r'(?m)^\s*cutoff\s*=\s*PME')
        self.assertTrue(os.path.exists(os.path.join(_EXAMPLES, 'ala_box.pdb')))


@unittest.skipUnless(_HAVE_OPENMM and _HAVE_OQP, 'needs OpenMM and the oqp package')
class TestDriverGates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import openmm.app as app
        cls.app = app
        cls.box = app.PDBFile(os.path.join(_EXAMPLES, 'ala_box.pdb'))
        cls.cluster = app.PDBFile(os.path.join(_EXAMPLES, 'formaldehyde_water.pdb'))

    def _bare(self, cutoff, topology):
        d = object.__new__(OpenQpQMMM)
        d.Cutoff = cutoff
        d.topology = topology
        return d

    def test_periodic_branch_only_for_periodic_methods(self):
        app = self.app
        for cut in (app.NoCutoff, app.CutoffNonPeriodic):
            d = self._bare(cut, self.box.topology)
            self.assertFalse(d._is_periodic())
            self.assertIsNone(d._box_lengths_bohr())     # box vectors present, still a cluster
        for cut in (app.PME, app.Ewald, app.CutoffPeriodic, app.LJPME):
            d = self._bare(cut, self.box.topology)
            self.assertTrue(d._is_periodic())
            np.testing.assert_allclose(d._box_lengths_bohr(), 16.0 * 1.8897259886 * np.ones(3), rtol=1e-9)

    def test_periodic_method_without_box_raises(self):
        top = self.app.Topology()                    # no CRYST1 -> no box vectors
        self.assertIsNone(top.getPeriodicBoxVectors())
        d = self._bare(self.app.PME, top)
        with self.assertRaisesRegex(ValueError, 'box vectors'):
            d._box_lengths_bohr()
        # a cluster PDB that happens to carry a CRYST1 record is still a cluster
        d = self._bare(self.app.NoCutoff, self.cluster.topology)
        self.assertIsNotNone(self.cluster.topology.getPeriodicBoxVectors())
        self.assertIsNone(d._box_lengths_bohr())

    def test_option_validation(self):
        common = dict(positions=None, topology=None, forcefield=None, qm_atoms=[0])
        for bad in (-0.5, 0.0 + float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'mm_charge_width'):
                OpenQpQMMM(mm_charge_width=bad, **common)
        for bad in (0.0, -1e-6, float('nan')):
            with self.assertRaisesRegex(ValueError, 'ewald_tol'):
                OpenQpQMMM(ewald_tol=bad, **common)
        # accepted values fall through to the usual "oqp_cfg or mol" check
        with self.assertRaisesRegex(ValueError, 'oqp_cfg'):
            OpenQpQMMM(mm_charge_width=0.7, ewald_tol=1e-6, **common)
        with self.assertRaisesRegex(ValueError, 'oqp_cfg'):
            OpenQpQMMM(mm_charge_width=0, ewald_tol=None, **common)

    def test_h_lj_leaves_qm_hydrogens_alone(self):
        import openmm as mm
        import openmm.unit as unit
        app = self.app
        ff = app.ForceField(os.path.join(_EXAMPLES, 'formaldehyde.xml'),
                            os.path.join(_EXAMPLES, 'tip3p.xml'))
        pdb = self.cluster
        d = object.__new__(OpenQpQMMM)
        d.positions, d.topology, d.forcefield = pdb.positions, pdb.topology, ff
        d.qm_atoms = np.array([0, 1, 2, 3])          # H2CO: QM hydrogens are atoms 2, 3
        d.Cutoff = app.NoCutoff
        d.ewald_tol, d.lj_switch, d.h_lj = None, False, True
        d.Embedding, d.espf_full = 'electrostatic', True
        ref = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                              constraints=None, rigidWater=False)
        nb_ref = next(f for f in ref.getForces() if isinstance(f, mm.NonbondedForce))
        systems = d.prepare_mm()
        nb = next(f for f in systems['sys0'].getForces() if isinstance(f, mm.NonbondedForce))
        n_mm_h_set = 0
        for atom in pdb.topology.atoms():
            if atom.element is None or atom.element.atomic_number != 1:
                continue
            _, s0, e0 = nb_ref.getParticleParameters(atom.index)
            _, s1, e1 = nb.getParticleParameters(atom.index)
            if atom.index in (2, 3):
                self.assertEqual(s1, s0)
                self.assertEqual(e1, e0)             # QM hydrogen untouched
            elif e0.value_in_unit(unit.kilojoule_per_mole) == 0.0:
                self.assertGreater(e1.value_in_unit(unit.kilojoule_per_mole), 0.0)
                n_mm_h_set += 1
        self.assertEqual(n_mm_h_set, 10)             # 5 TIP3P waters x 2 H


@unittest.skipUnless(_HAVE_OPENMM and _HAVE_OQP, 'needs OpenMM and the oqp package')
class TestNamdLinkAtoms(unittest.TestCase):
    def test_fold_link_charges_conserves_total_charge(self):
        n = object.__new__(NAMD_QMMM)
        n.nqm = 3
        n.link_atoms = [types.SimpleNamespace(host_row=1, g=0.7, qm_index=9, mm_index=8),
                        types.SimpleNamespace(host_row=2, g=0.7, qm_index=10, mm_index=12)]
        pchg = np.array([0.10, -0.20, 0.30, 0.05, -0.07])
        q = n._fold_link_charges(pchg)
        np.testing.assert_allclose(q, [0.10, -0.15, 0.23])
        self.assertAlmostEqual(q.sum(), pchg.sum(), places=14)
        n.link_atoms = []
        np.testing.assert_array_equal(n._fold_link_charges(pchg[:3]), pchg[:3])

    def test_restart_identity_tracks_periodic_options(self):
        import openmm as mm
        import openmm.app as app
        pdb = app.PDBFile(os.path.join(_EXAMPLES, 'formaldehyde_water.pdb'))
        ff = app.ForceField(os.path.join(_EXAMPLES, 'formaldehyde.xml'),
                            os.path.join(_EXAMPLES, 'tip3p.xml'))
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
        n = object.__new__(NAMD_QMMM)
        n.pdb, n._mm = pdb, mm
        n.qm_atoms = np.array([0, 1, 2, 3])
        n.natom_all = pdb.topology.getNumAtoms()
        n.m_all = np.ones(n.natom_all)
        base = {'embedding': 'electrostatic', 'cutoff': 'NoCutoff'}
        ref = n._qmmm_wham_system_identity(system, dict(base))['sha256']
        self.assertEqual(n._qmmm_wham_system_identity(system, dict(base))['sha256'], ref)
        for key, value in (('mm_charge_width', '0.7'), ('h_lj', True),
                           ('lj_switch', True), ('ewald_tol', '1e-6')):
            cfg = dict(base); cfg[key] = value
            self.assertNotEqual(n._qmmm_wham_system_identity(system, cfg)['sha256'], ref, key)


if __name__ == '__main__':
    unittest.main()
