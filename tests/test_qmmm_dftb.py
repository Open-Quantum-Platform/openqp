"""Finite-difference consistency of the DFTB QM/MM (ESPF/OpenMM) driver.

These tests are THE arbiter that the ``method=dftb`` energy/force bookkeeping in
``oqp.library.qmmm_driver.OpenQpQMMM`` is internally consistent: no double
counting and no missing nuclear/coupling term. For a chosen atom (QM or MM) we
displace it by a small step along a random direction and check that the
projected total QM/MM force equals minus the central finite difference of the
total QM/MM energy.

The DFTB electrostatic-embedding contract (see oqp.library.openqp_dftb): setting
``mol.dftb_external_potential`` folds the per-atom MM potential straight into the
SCC Hamiltonian, so the returned state energy already contains the full
net-charge coupling ``E_ext = sum_A q_A phi_A`` (no separate nuclear
``sum_A Z_A phi_A`` term is added for dftb, unlike the native ESPF path) and the
returned gradient is ``d(E_embedded)/dR`` at fixed potential (the Pulay coupling
term is inside it, so ``grad_esp_qmmm`` / ``OQP::ESPF_GRAD`` are skipped). The
classical dphi/dR forces remain the driver's job via ``_coupling_forces``, which
reads the ``OQP::partial_charges`` the adapter publishes after the gradient call.

Both a whole-molecule QM region (H2 solute in TIP3P waters) and a covalent
QM/MM boundary with a hydrogen link atom (ethane split QM=CH3 / MM=CH3) are
exercised, including atoms on both sides of the boundary and the link-atom host.

The system topologies and force fields are built programmatically (a minimal
custom ForceField XML plus TIP3P), so the driver reads MM charges from the
NonbondedForce exactly as in a force-field-typed run without depending on
amber14 residue templates for the non-standard QM residues.

Requires OpenMM, the compiled OpenQP runtime (OPENQP_ROOT), the standalone
openqp-dftb shared library (libopenqp_dftb_c), and the minimal DFTB parameter
fixtures; the tests skip cleanly when any of these is absent.
"""

import os
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

# Minimal DFTB parameter fixtures shipped with the openqp-dftb repo. Overridable
# so the test is not pinned to one checkout; skipped when they are missing.
_HH_PARAM = os.environ.get(
    "OPENQP_DFTB_HH_PARAM",
    "/Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-dftb/tests/data/minimal_hh.opdftb",
)
_CH_SKF_DIR = os.environ.get(
    "OPENQP_DFTB_CH_SKF_DIR",
    "/Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-dftb/tests/data/minimal_skf_dir",
)

try:
    import openmm as mm  # noqa: F401
    import openmm.app as app
    from openmm import unit as u
    from openmm.app import element
    _HAVE_OPENMM = True
except Exception:  # pragma: no cover - optional dependency
    _HAVE_OPENMM = False

try:
    from oqp.library.qmmm_driver import OpenQpQMMM, _normalize_embedding
    _HAVE_OQP = True
except Exception:  # pragma: no cover - uncompiled backend / missing OPENQP_ROOT
    _HAVE_OQP = False


_ANG2NM = 0.1


def _e_f(driver, topology, qm_atoms, positions):
    """Total QM/MM energy (kJ/mol) and force array (kJ/mol/nm) at *positions*.

    Mirrors the QMMM_MD/NAMD caller contract: the MM context positions are
    synced before compute_force (the driver reads the MM-MM energy from the
    OpenMM context, whereas the QM geometry and embedding come from the passed
    positions)."""
    driver.mm_systems["sim0"].context.setPositions(positions)
    e, f = driver.compute_force(positions, topology, driver.mm_systems, qm_atoms)
    ev = e.value_in_unit(u.kilojoule_per_mole) if hasattr(e, "value_in_unit") else float(e)
    fv = (f.value_in_unit(u.kilojoule_per_mole / u.nanometer)
          if hasattr(f, "value_in_unit") else np.asarray(f))
    return ev, np.asarray(fv)


def _fd_residuals(driver, topology, qm_atoms, positions, atoms, *, step_ang=1e-4, seed=0):
    """Directional finite-difference residual per tested atom.

    Displaces each atom by ``step_ang`` Angstrom along a random unit vector and
    returns a list of (atom, |F_analytic|, residual) where residual =
    |(-dE/ds) - F.dot(u)|, i.e. the mismatch between the projected analytic
    force and the central finite difference of the total energy."""
    rng = np.random.default_rng(seed)
    h = step_ang * _ANG2NM
    p0 = np.array(positions.value_in_unit(u.nanometer))
    _, f0 = _e_f(driver, topology, qm_atoms, positions)
    out = []
    for a in atoms:
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        pp = p0.copy(); pp[a] += h * d
        pm = p0.copy(); pm[a] -= h * d
        ep, _ = _e_f(driver, topology, qm_atoms, pp * u.nanometer)
        em, _ = _e_f(driver, topology, qm_atoms, pm * u.nanometer)
        fd = -(ep - em) / (2.0 * h)
        analytic = float(np.dot(f0[a], d))
        out.append((a, float(np.linalg.norm(f0[a])), abs(fd - analytic)))
    return out


# --- custom force fields / geometries --------------------------------------

_H2_XML = """<ForceField>
 <AtomTypes>
  <Type name="QMH" class="QMH" element="H" mass="1.008"/>
 </AtomTypes>
 <Residues>
  <Residue name="HHQ">
   <Atom name="H1" type="QMH"/>
   <Atom name="H2" type="QMH"/>
   <Bond atomName1="H1" atomName2="H2"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="QMH" charge="0.0" sigma="0.25" epsilon="0.1"/>
 </NonbondedForce>
</ForceField>
"""

_TIP3P_XML = """<ForceField>
 <AtomTypes>
  <Type name="tip3p-O" class="OW" element="O" mass="15.99943"/>
  <Type name="tip3p-H" class="HW" element="H" mass="1.007947"/>
 </AtomTypes>
 <Residues>
  <Residue name="HOH">
   <Atom name="O" type="tip3p-O"/>
   <Atom name="H1" type="tip3p-H"/>
   <Atom name="H2" type="tip3p-H"/>
   <Bond atomName1="O" atomName2="H1"/>
   <Bond atomName1="O" atomName2="H2"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond class1="OW" class2="HW" length="0.09572" k="462750.4"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle class1="HW" class2="OW" class3="HW" angle="1.82421813418" k="836.8"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="tip3p-O" charge="-0.834" sigma="0.31507524065751241" epsilon="0.635968"/>
  <Atom type="tip3p-H" charge="0.417" sigma="1" epsilon="0"/>
 </NonbondedForce>
</ForceField>
"""

# Ethane residue with real MM point charges so the boundary carries a field.
_ETH_XML = """<ForceField>
 <AtomTypes>
  <Type name="ec" class="ec" element="C" mass="12.011"/>
  <Type name="eh" class="eh" element="H" mass="1.008"/>
 </AtomTypes>
 <Residues>
  <Residue name="ETH">
   <Atom name="C1" type="ec"/><Atom name="H11" type="eh"/><Atom name="H12" type="eh"/><Atom name="H13" type="eh"/>
   <Atom name="C2" type="ec"/><Atom name="H21" type="eh"/><Atom name="H22" type="eh"/><Atom name="H23" type="eh"/>
   <Bond atomName1="C1" atomName2="C2"/>
   <Bond atomName1="C1" atomName2="H11"/><Bond atomName1="C1" atomName2="H12"/><Bond atomName1="C1" atomName2="H13"/>
   <Bond atomName1="C2" atomName2="H21"/><Bond atomName1="C2" atomName2="H22"/><Bond atomName1="C2" atomName2="H23"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond class1="ec" class2="ec" length="0.1535" k="253634.0"/>
  <Bond class1="ec" class2="eh" length="0.1092" k="284512.0"/>
 </HarmonicBondForce>
 <NonbondedForce coulomb14scale="0.8333" lj14scale="0.5">
  <Atom type="ec" charge="-0.18" sigma="0.339967" epsilon="0.457730"/>
  <Atom type="eh" charge="0.06"  sigma="0.264953" epsilon="0.065710"/>
 </NonbondedForce>
</ForceField>
"""


def _write(tmpdir, name, text):
    path = os.path.join(tmpdir, name)
    with open(path, "w", encoding="ascii") as handle:
        handle.write(text)
    return path


def _add_molecule(top, resname, names, elems, bonds):
    chain = top.addChain()
    res = top.addResidue(resname, chain)
    atoms = [top.addAtom(n, e, res) for n, e in zip(names, elems)]
    for i, j in bonds:
        top.addBond(atoms[i], atoms[j])
    return atoms


def _h2_water_system():
    """QM H2 solute in a small cluster of TIP3P waters (non-periodic)."""
    top = app.Topology()
    coords = []
    # QM H2 (indices 0,1)
    _add_molecule(top, "HHQ", ["H1", "H2"], [element.hydrogen] * 2, [(0, 1)])
    coords += [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]
    # two TIP3P waters at asymmetric offsets
    for base in ([2.6, 0.3, 0.2], [-2.2, 1.8, -1.1]):
        _add_molecule(top, "HOH", ["O", "H1", "H2"],
                      [element.oxygen, element.hydrogen, element.hydrogen],
                      [(0, 1), (0, 2)])
        bx, by, bz = base
        coords += [[bx, by, bz], [bx + 0.5, by - 0.7, bz + 0.4],
                   [bx + 0.45, by + 0.75, bz - 0.3]]
    positions = (np.array(coords) * _ANG2NM) * u.nanometer
    return top, positions, [0, 1]


def _ethane_link_system():
    """Ethane split QM=CH3(+link H) / MM=CH3 across the C-C bond.

    The C-C bond length is chosen so the hydrogen link atom lands at the DFTB
    C-H equilibrium (~1.09 A), i.e. the QM fragment is a symmetric CH4 that the
    minimal C/H SK tables converge cleanly."""
    c1 = np.array([0.0, 0.0, 0.0])
    c2 = np.array([0.0, 0.0, 1.4854])   # link atom lands at 1.09 A (g_C-C=0.7338)
    h1 = [c1 + np.array([1.028, 0.0, -0.363]),
          c1 + np.array([-0.514, 0.890, -0.363]),
          c1 + np.array([-0.514, -0.890, -0.363])]
    h2 = [c2 + np.array([1.028, 0.0, 0.363]),
          c2 + np.array([-0.514, 0.890, 0.363]),
          c2 + np.array([-0.514, -0.890, 0.363])]
    coords = np.array([c1, h1[0], h1[1], h1[2], c2, h2[0], h2[1], h2[2]])
    elems = [element.carbon] + [element.hydrogen] * 3 + [element.carbon] + [element.hydrogen] * 3
    names = ["C1", "H11", "H12", "H13", "C2", "H21", "H22", "H23"]
    bonds = [(0, 4), (0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)]
    top = app.Topology()
    _add_molecule(top, "ETH", names, elems, bonds)
    positions = (coords * _ANG2NM) * u.nanometer
    return top, positions, [0, 1, 2, 3]   # QM = C1 + its three H


def _dftb_cfg(parameter_path, *, excited=False):
    cfg = {
        "input.method": "dftb", "input.runtype": "grad", "input.basis": "6-31g",
        "input.functional": "", "input.charge": 0,
        "dftb.parameter_path": parameter_path, "dftb.scc_tolerance": 1e-10,
    }
    if excited:
        cfg.update({
            "properties.grad": "1", "dftb.type": "mrsf",
            "tdhf.type": "mrsf", "tdhf.nstate": 2,
            "scf.type": "rohf", "scf.multiplicity": 3,
        })
    else:
        cfg.update({"properties.grad": "0", "dftb.type": "ground", "scf.type": "rhf"})
    return cfg


@unittest.skipUnless(_HAVE_OQP, "compiled OpenQP backend unavailable")
class TestEmbeddingNormalization(unittest.TestCase):
    """`_normalize_embedding` is a pure function -- no OpenMM, no SCF.

    Kept out of the OpenMM-gated class below so these actually run: gating a
    string-validation test behind a simulation backend means it never executes
    anywhere, which is how the typo-tolerance it guards got in.
    """

    def test_unknown_embedding_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown QM/MM embedding"):
            _normalize_embedding("electrostatc")

    def test_embedding_is_case_normalized(self):
        self.assertEqual(_normalize_embedding(" Mechanical "), "mechanical")

    def test_every_dispatched_embedding_is_accepted(self):
        # The dispatch in compute_force/espf_full branches on these exact
        # spellings; if one were missing from _VALID_EMBEDDINGS the driver
        # would reject a mode it can actually run.
        for embedding in ("mechanical", "electrostatic", "espf",
                          "espf_full", "split"):
            self.assertEqual(_normalize_embedding(embedding), embedding)


@unittest.skipUnless(_HAVE_OPENMM and _HAVE_OQP,
                     "OpenMM or compiled OpenQP backend unavailable")
class TestQMMMDFTBEmbeddingDispatch(unittest.TestCase):
    def test_mechanical_embedding_passes_none_to_dftb(self):
        driver = object.__new__(OpenQpQMMM)
        driver.use_mol = True
        driver.Embedding = "mechanical"
        driver.mol = SimpleNamespace(
            config={"input": {"method": "dftb"}}
        )
        driver._update_mol_positions = mock.Mock()
        expected = object()
        driver._forces_qm_dftb = mock.Mock(return_value=expected)

        result = driver.forces_qm_openqp(
            potmm=np.zeros(2), potqm=np.zeros((2, 2))
        )

        self.assertIs(result, expected)
        driver._forces_qm_dftb.assert_called_once_with(driver.mol, None)

    def test_electrostatic_embedding_forwards_the_field_to_dftb(self):
        # Counterpart to the mechanical case: without this, a dispatch that
        # passed None unconditionally would still satisfy the test above.
        driver = object.__new__(OpenQpQMMM)
        driver.use_mol = True
        driver.Embedding = "electrostatic"
        driver.mol = SimpleNamespace(config={"input": {"method": "dftb"}})
        driver._update_mol_positions = mock.Mock()
        driver._forces_qm_dftb = mock.Mock(return_value=object())
        potmm = np.arange(2, dtype=float)

        driver.forces_qm_openqp(potmm=potmm, potqm=np.zeros((2, 2)))

        forwarded = driver._forces_qm_dftb.call_args[0][1]
        self.assertIsNotNone(forwarded)
        np.testing.assert_allclose(forwarded, potmm)

    def test_zero_embedding_matches_the_electrostatic_array_shapes(self):
        # _zero_embedding must be sized like the arrays the ESPF path builds
        # (nqm + nlink), or the einsum against espf_op_corr silently mismatches.
        driver = object.__new__(OpenQpQMMM)
        driver.qm_atoms = np.array([0, 1, 2])
        driver.link_atoms = [object(), object()]

        potmm, potqm = driver._zero_embedding()

        nqm_c = len(driver.qm_atoms) + len(driver.link_atoms)
        self.assertEqual(potmm.shape, (nqm_c,))
        self.assertEqual(potqm.shape, (nqm_c, nqm_c))
        self.assertEqual(potmm.dtype, np.float64)
        self.assertEqual(potqm.dtype, np.float64)
        # "Zero field" is the whole correctness argument -- assert it.
        self.assertFalse(potmm.any())
        self.assertFalse(potqm.any())


@unittest.skipUnless(_HAVE_OPENMM and _HAVE_OQP, "OpenMM or compiled OpenQP backend unavailable")
@unittest.skipUnless(os.path.exists(_HH_PARAM), f"missing DFTB fixture {_HH_PARAM}")
@unittest.skipUnless(os.path.exists(_CH_SKF_DIR), f"missing DFTB SKF dir {_CH_SKF_DIR}")
class TestQMMMDFTB(unittest.TestCase):
    # FD consistency tolerance: a component passes when the residual is below a
    # small absolute floor OR a relative fraction of the analytic force. The
    # residual is dominated by SCC-charge finite-difference noise and, along
    # stiff bond coordinates, central-difference truncation -- both far below
    # any bookkeeping error (a missing coupling/nuclear term yields a residual
    # ~100% of the force).
    ABS_FLOOR = 5.0     # kJ/mol/nm
    REL = 0.02          # 2% of |F_atom|

    def _assert_consistent(self, residuals, label):
        for atom, fmag, resid in residuals:
            tol = max(self.ABS_FLOOR, self.REL * fmag)
            self.assertLess(
                resid, tol,
                f"{label}: atom {atom} FD residual {resid:.4g} kJ/mol/nm "
                f"exceeds tol {tol:.4g} (|F|={fmag:.4g})",
            )

    def _driver(self, tmpdir, xmls, top, positions, qm_atoms, cfg):
        ff = app.ForceField(*[_write(tmpdir, f"ff{k}.xml", x) for k, x in enumerate(xmls)])
        return OpenQpQMMM(positions, top, ff, qm_atoms, oqp_cfg=cfg,
                          Cutoff=app.NoCutoff, Embedding="electrostatic")

    def test_h2_in_water_ground_state(self):
        """Whole-molecule QM region: ground-state DFTB energy/force consistency
        for QM and MM atoms (the primary bookkeeping arbiter)."""
        import tempfile
        top, positions, qm = _h2_water_system()
        with tempfile.TemporaryDirectory() as tmp:
            drv = self._driver(tmp, [_H2_XML, _TIP3P_XML], top, positions, qm,
                               _dftb_cfg(_HH_PARAM))
            self.assertEqual(len(drv.link_atoms), 0)
            # QM H (0), MM water O (2) and MM water H (6)
            res = _fd_residuals(drv, top, qm, positions, [0, 1, 2, 6])
            self._assert_consistent(res, "H2/water ground")

    def test_ethane_link_atom_ground_state(self):
        """Covalent QM/MM boundary: ground-state DFTB energy/force consistency
        including the QM host (C1), a QM hydrogen, the MM host (C2) across the
        boundary, and an MM hydrogen -- exercising the hydrogen link atom and
        its gradient/coupling-force projection onto the real host atoms."""
        import tempfile
        top, positions, qm = _ethane_link_system()
        with tempfile.TemporaryDirectory() as tmp:
            drv = self._driver(tmp, [_ETH_XML], top, positions, qm,
                               _dftb_cfg(_CH_SKF_DIR))
            self.assertEqual(len(drv.link_atoms), 1)
            self.assertEqual(drv.link_atoms[0].qm_index, 0)   # QM frontier C1
            self.assertEqual(drv.link_atoms[0].mm_index, 4)   # MM host C2
            # C1 (QM host, index 0), QM H (1), MM host C2 (4), MM H (5)
            res = _fd_residuals(drv, top, qm, positions, [0, 1, 4, 5])
            self._assert_consistent(res, "ethane link ground")

    def test_h2_in_water_excited_state(self):
        """Excited-state (MRSF-TDDFTB S1) QM/MM consistency on the H2/water box.

        The coupling force uses the relaxed excited-state charges the adapter
        publishes in OQP::partial_charges. Excited-state Z-vector charge
        relaxation is being finished in a parallel change, so a modestly larger
        residual than the ground state is expected here; the check uses a looser
        tolerance and is informational rather than a hard bookkeeping arbiter.
        """
        import tempfile
        top, positions, qm = _h2_water_system()
        with tempfile.TemporaryDirectory() as tmp:
            drv = self._driver(tmp, [_H2_XML, _TIP3P_XML], top, positions, qm,
                               _dftb_cfg(_HH_PARAM, excited=True))
            res = _fd_residuals(drv, top, qm, positions, [0, 1, 2, 6])
            worst = max(r for _, _, r in res)
            # Loose informational bound (see docstring); tighten once excited
            # Z-vector relaxed charges land.
            for atom, fmag, resid in res:
                tol = max(15.0, 0.05 * fmag)
                self.assertLess(
                    resid, tol,
                    f"excited: atom {atom} residual {resid:.4g} exceeds "
                    f"loose tol {tol:.4g} (|F|={fmag:.4g})",
                )
            self.assertTrue(np.isfinite(worst))


if __name__ == "__main__":
    unittest.main()
