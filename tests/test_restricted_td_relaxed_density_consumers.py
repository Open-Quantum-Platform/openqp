"""Regression guard: restricted TDHF/TDDFT (RPA, TDA) relaxed-density consumers.

The restricted Z-vector solver (tdhf_z_vector) stores OQP::td_p with ONE
column, the alpha part of a closed-shell relaxed difference density; its spin
sum is 2*td_p(:,1), the density tdhf_gradient uses (dmat_a + 2*p(:,1)).  SF and
MRSF store alpha and beta columns.  The excited-state consumers all added
td_p(:,1) + td_p(:,2) and required OQP::td_abxc, which only SF/MRSF create:

* gas phase, ``[properties] td_prop=true`` with RPA/TDA: the relaxed dipole read
  past the end of td_p, then mulliken_excited stopped with
  "Record OQP::td_abxc not found";
* QM/MM (ESPF) RPA/TDA excited-state gradient: grad_esp_qmmm_excited stopped
  with the same message, and the restricted Z-vector never formed the relaxed
  ESPF charges the MM forces use (the MRSF Z-vector does).

Tests:

* ``GasPhaseTdPropTests``: RPA and TDA gradients with td_prop finish, and the
  reported relaxed excited-state dipole equals an independent recomputation
  from DM_A + 2*td_p and the AO dipole integrals (OQP::td_dip_ao of an MRSF run
  on the same geometry and basis).  A factor 1 instead of 2 moves the dipole by
  far more than the tolerance (asserted).
* ``QMMMRestrictedExcitedForceTests`` (needs OpenMM): RPA and TDA excited-state
  QM/MM forces on a TIP3P water dimer match central finite differences of the
  excited-state energy for a QM and an MM coordinate.  The same check with the
  ground-state ESPF density (ESPF_ROHF=1, i.e. no relaxed charges, which is what
  the restricted path used before) is violated by much more than the
  tolerance, so the test sees missing relaxed charges.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
OPENMM_REQUIRED = os.getenv("OQP_REQUIRE_OPENMM_TESTS") == "1"
_RUNTIME_ERROR = ""


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()
try:
    import openmm  # noqa: F401
    HAVE_OPENMM = True
except Exception:  # optional dependency
    HAVE_OPENMM = False

DEBYE_PER_AU = 2.541746473

H2O = """  8    0.000000000    0.000000000   -0.041061554
  1   -0.533194329    0.533194329   -0.614469223
  1    0.533194329   -0.533194329   -0.614469223"""

RESTRICTED_GRAD_INPUT = """[input]
system=
{geom}
charge=0
method=tdhf
basis=6-31g
runtype=grad
functional=
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rhf
multiplicity=1
conv=1.0e-10
save_molden=False

[properties]
grad=1
td_prop=true

[tdhf]
type={td}
multiplicity=1
nstate=2
conv=1.0e-10
zvconv=1.0e-10
"""

# Same geometry and basis, MRSF energy: only used for OQP::td_dip_ao (AO dipole
# integrals at the centre of mass), which depend on geometry and basis alone.
MRSF_ENERGY_INPUT = """[input]
system=
{geom}
charge=0
method=tdhf
basis=6-31g
runtype=energy
functional=
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rohf
multiplicity=3
conv=1.0e-8
save_molden=False

[tdhf]
type=mrsf
multiplicity=1
nstate=2
conv=1.0e-8
"""

# Runs both inputs in one interpreter and prints the tags the check needs.
_GAS_DRIVER = r"""
import json, sys
import numpy as np
from oqp.pyoqp import Runner
grad_inp, mrsf_inp = sys.argv[1], sys.argv[2]
r = Runner(project="grad", input_file=grad_inp, log="grad.log", silent=1, usempi=False)
r.run()
dma = np.asarray(r.mol.data["OQP::DM_A"], dtype=float).ravel()
tdp = np.asarray(r.mol.data["OQP::td_p"], dtype=float)
xyz = np.asarray(r.mol.get_atoms2("coords"), dtype=float).reshape(-1, 3)
zn = np.asarray(r.mol.get_atoms(), dtype=float).ravel()
mass = np.asarray(r.mol.get_atoms2("mass"), dtype=float).ravel()
m = Runner(project="mrsf", input_file=mrsf_inp, log="mrsf.log", silent=1, usempi=False)
m.run()
dip_ao = np.asarray(m.mol.data["OQP::td_dip_ao"], dtype=float)
print("RESULT " + json.dumps({"dm_a": dma.tolist(), "td_p": tdp.tolist(),
      "td_p_shape": list(tdp.shape), "xyz": xyz.tolist(), "zn": zn.tolist(),
      "mass": mass.tolist(), "dip_ao": dip_ao.tolist(),
      "dip_ao_shape": list(dip_ao.shape)}))
"""


def _packed_trace(m_packed, d_packed, nbf):
    """traceprod_sym_packed: sum_ij M_ij D_ij for packed L-triangles."""
    import numpy as np
    m = np.zeros((nbf, nbf))
    d = np.zeros((nbf, nbf))
    il = np.tril_indices(nbf)
    m[il] = m_packed
    d[il] = d_packed
    m = m + m.T - np.diag(np.diag(m))
    d = d + d.T - np.diag(np.diag(d))
    return float(np.sum(m * d))


def _excited_dipole_from_log(text):
    """Excited-state dipole (Debye) printed by electric_moments_excited."""
    i = text.rfind("electric dipole (Debye) - excited state")
    if i < 0:
        return None
    nums = re.findall(r"-?\d+\.\d+", text[i:i + 600])
    return [float(x) for x in nums[:3]]


class RestrictedTdBuildGate(unittest.TestCase):
    @unittest.skipUnless(NATIVE_TESTS_REQUIRED,
                         "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP")
    def test_required_runtime(self):
        self.assertTrue(RUNTIME_AVAILABLE, f"compiled OpenQP runtime is required: {_RUNTIME_ERROR}")

    @unittest.skipUnless(OPENMM_REQUIRED, "set OQP_REQUIRE_OPENMM_TESTS=1 where OpenMM is installed")
    def test_required_openmm(self):
        self.assertTrue(HAVE_OPENMM, "OpenMM is required (OQP_REQUIRE_OPENMM_TESTS=1)")


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class GasPhaseTdPropTests(unittest.TestCase):
    def _check(self, td):
        import numpy as np
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        work = tmp.name
        grad_inp = os.path.join(work, "grad.inp")
        mrsf_inp = os.path.join(work, "mrsf.inp")
        with open(grad_inp, "w") as fh:
            fh.write(RESTRICTED_GRAD_INPUT.format(geom=H2O, td=td))
        with open(mrsf_inp, "w") as fh:
            fh.write(MRSF_ENERGY_INPUT.format(geom=H2O))
        env = dict(os.environ)
        env.setdefault("OMP_NUM_THREADS", "1")
        proc = subprocess.run([sys.executable, "-c", _GAS_DRIVER, grad_inp, mrsf_inp],
                              cwd=work, env=env, capture_output=True, text=True, timeout=900)
        detail = (f"returncode={proc.returncode}\n--- stdout ---\n{proc.stdout[-3000:]}"
                  f"\n--- stderr ---\n{proc.stderr[-3000:]}")
        self.assertEqual(proc.returncode, 0, detail)
        lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
        self.assertEqual(len(lines), 1, detail)
        d = json.loads(lines[0][7:])
        with open(os.path.join(work, "grad.log"), errors="replace") as fh:
            reported = _excited_dipole_from_log(fh.read())
        self.assertIsNotNone(reported, "excited-state dipole missing from the log")

        tdp = np.asarray(d["td_p"]).reshape(d["td_p_shape"])
        self.assertEqual(tdp.ndim, 2)
        self.assertEqual(min(tdp.shape), 1, "restricted td_p must have one column")
        tdp = tdp.reshape(-1)
        dma = np.asarray(d["dm_a"])
        nbf2 = dma.size
        nbf = int(round((np.sqrt(8 * nbf2 + 1) - 1) / 2))
        # OQP::td_dip_ao is a Fortran (nbf2, 3) array; the bridge hands back
        # its column-major buffer, so column k is the k-th contiguous block.
        dip_ao = np.asarray(d["dip_ao"]).reshape(-1)
        self.assertEqual(dip_ao.size, 3 * nbf2)
        dip_ao = dip_ao.reshape(3, nbf2)

        xyz = np.asarray(d["xyz"])
        zn = np.asarray(d["zn"])
        mass = np.asarray(d["mass"])
        com = (xyz * mass[:, None]).sum(0) / mass.sum()

        def dipole(scale):
            dens = dma + scale * tdp
            nuc = ((xyz - com) * zn[:, None]).sum(0)
            return np.array([nuc[k] - _packed_trace(dip_ao[k], dens, nbf) for k in range(3)])

        # Orientation guard: this H2O geometry is C2v with the axis along z,
        # so the ground-state dipole has no x or y component.
        ground = dipole(0.0) * DEBYE_PER_AU
        self.assertLess(np.max(np.abs(ground[:2])), 1.0e-8,
                        f"dipole integrals misread (ground dipole {ground})")
        ref = dipole(2.0) * DEBYE_PER_AU
        wrong = dipole(1.0) * DEBYE_PER_AU
        got = np.asarray(reported)
        # The log prints 6 decimals in Debye.
        self.assertLess(np.max(np.abs(got - ref)), 2.0e-5,
                        f"{td}: reported {got} vs DM_A + 2*td_p {ref}")
        self.assertGreater(np.max(np.abs(wrong - ref)), 1.0e-2,
                           f"{td}: test cannot tell factor 1 from 2 ({wrong} vs {ref})")

    def test_rpa_td_prop_relaxed_dipole(self):
        self._check("rpa")

    def test_tda_td_prop_relaxed_dipole(self):
        self._check("tda")


def _water_dimer():
    """TIP3P water dimer as an in-memory OpenMM topology (positions in nm)."""
    import openmm.app as app
    import openmm.unit as unit
    from openmm import Vec3
    top = app.Topology()
    chain = top.addChain()
    pos = []
    for shift in (0.0, 0.30):
        res = top.addResidue("HOH", chain)
        o = top.addAtom("O", app.element.oxygen, res)
        h1 = top.addAtom("H1", app.element.hydrogen, res)
        h2 = top.addAtom("H2", app.element.hydrogen, res)
        top.addBond(o, h1)
        top.addBond(o, h2)
        pos += [Vec3(0.0, 0.0, shift), Vec3(0.09572, 0.0, shift),
                Vec3(-0.02400, 0.09267, shift)]
    return top, unit.Quantity(pos, unit.nanometer)


# Energy and forces of the [properties] grad state at the base geometry and at
# +-h displacements of the requested (atom, xyz) coordinates; one driver per
# geometry so every point starts from the same guess.
_QMMM_FD_DRIVER = r"""
import copy, json, sys
import numpy as np
import openmm.app as app
import openmm.unit as unit
from openmm import Vec3
from oqp.library.qmmm_driver import OpenQpQMMM
exec(sys.argv[1])            # defines _water_dimer()
cfg, coords, h = json.loads(sys.argv[2]), json.loads(sys.argv[3]), float(sys.argv[4])
top, pos0 = _water_dimer()
ff = app.ForceField("tip3p.xml")
base = np.array([[p.x, p.y, p.z] for p in pos0.value_in_unit(unit.nanometer)])

def point(xyz):
    pos = unit.Quantity([Vec3(*r) for r in xyz], unit.nanometer)
    drv = OpenQpQMMM(pos, top, ff, [0, 1, 2], oqp_cfg=copy.deepcopy(cfg),
                     Cutoff=app.NoCutoff, Embedding="electrostatic")
    drv.mm_systems["sim0"].context.setPositions(pos)
    e, f = drv.compute_force(pos, top, drv.mm_systems, [0, 1, 2])
    e = float(e.value_in_unit(unit.kilojoule_per_mole) if hasattr(e, "value_in_unit") else e)
    f = (f.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
         if hasattr(f, "value_in_unit") else f)
    return e, np.asarray(f, dtype=float).reshape(-1, 3)

e0, f0 = point(base)
fd = []
for a, k in coords:
    xp = base.copy(); xp[a, k] += h
    xm = base.copy(); xm[a, k] -= h
    ep, _ = point(xp)
    em, _ = point(xm)
    fd.append({"atom": a, "xyz": k, "analytic": float(f0[a, k]),
               "fd": float(-(ep - em) / (2 * h))})
print("RESULT " + json.dumps({"energy": e0, "fd": fd}))
"""

# input.qmmm_flag=True as in every QM/MM input: without it the Z-vector solvers
# do not form the relaxed ESPF charges at all.
QMMM_RESTRICTED_CFG = {
    "input.qmmm_flag": "True", "input.method": "tdhf", "input.runtype": "grad", "input.basis": "sto-3g",
    "input.functional": "", "input.charge": 0, "input.d4": "False",
    "guess.type": "huckel", "scf.type": "rhf", "scf.multiplicity": 1,
    "scf.conv": "1.0e-10", "properties.grad": "1",
    "tdhf.multiplicity": 1, "tdhf.nstate": 2,
    "tdhf.conv": "1.0e-10", "tdhf.zvconv": "1.0e-10",
}


@unittest.skipUnless(RUNTIME_AVAILABLE and HAVE_OPENMM,
                     "compiled OpenQP runtime or OpenMM not available")
class QMMMRestrictedExcitedForceTests(unittest.TestCase):
    # MM water O (atom 3) along the QM-MM axis, and the QM O (atom 0).
    COORDS = [[3, 2], [3, 0], [0, 2]]
    H = 1.0e-4          # nm
    # kJ/mol/nm.  The HF ground state shows the same residual on QM atoms
    # (0.019, ESPF grid weights); MM atoms agree to ~5e-4.
    TOL = 0.05

    def _fd(self, td, env_extra=None):
        import inspect
        cfg = dict(QMMM_RESTRICTED_CFG, **{"tdhf.type": td})
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = "1"
        env.pop("ESPF_ROHF", None)
        if env_extra:
            env.update(env_extra)
        proc = subprocess.run(
            [sys.executable, "-c", _QMMM_FD_DRIVER, inspect.getsource(_water_dimer),
             json.dumps(cfg), json.dumps(self.COORDS), repr(self.H)],
            cwd=tmp.name, env=env, capture_output=True, text=True, timeout=1800)
        detail = (f"returncode={proc.returncode}\n--- stdout ---\n{proc.stdout[-3000:]}"
                  f"\n--- stderr ---\n{proc.stderr[-3000:]}")
        self.assertEqual(proc.returncode, 0, detail)
        lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
        self.assertEqual(len(lines), 1, detail)
        return json.loads(lines[0][7:])["fd"]

    def _check(self, td):
        res = self._fd(td)
        for r in res:
            self.assertLess(abs(r["analytic"] - r["fd"]), self.TOL,
                            f"{td} atom {r['atom']} xyz {r['xyz']}: analytic "
                            f"{r['analytic']:.5f} vs FD {r['fd']:.5f} kJ/mol/nm")
        # Control: without the relaxed ESPF density the MM-atom force is wrong
        # by much more than TOL, so the check above can see it.
        ctl = self._fd(td, {"ESPF_ROHF": "1"})
        worst = max(abs(r["analytic"] - r["fd"]) for r in ctl if r["atom"] == 3)
        self.assertGreater(worst, 5 * self.TOL,
                           f"{td}: ground-state ESPF density passes the FD check "
                           f"(max |dF| {worst:.4f}); the test is not discriminating")

    def test_rpa_qmmm_excited_force_matches_finite_difference(self):
        self._check("rpa")

    def test_tda_qmmm_excited_force_matches_finite_difference(self):
        self._check("tda")


if __name__ == "__main__":
    unittest.main()
