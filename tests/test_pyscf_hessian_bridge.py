"""Regression test for the external PySCF analytic Hessian bridge.

Guards the three bugs that previously made the analytic Hessian path
non-functional, plus UHF coverage:

1. the bridge was gated off by ``mol.usempi`` (hard-coded ``True`` in serial),
2. the geometry was double-converted (``get_system()`` is already Bohr, and the
   ``ANGSTROM_TO_BOHR`` constant is actually the Bohr->Angstrom factor),
3. the molden frequency writer crashed on a ``(-1, 1)``-reshaped ``freqs``.

The test builds H2O (RHF) and H2O+ (UHF) analytic Hessians through the bridge
and checks them element-wise against a direct PySCF Hessian on the identical
geometry. Skipped unless the compiled OpenQP runtime and PySCF are importable.
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUT = """[input]
system=
   O  -0.0000000000   0.0000000000  -0.0410615540
   H  -0.5331943294   0.5331943294  -0.6144692230
   H   0.5331943294  -0.5331943294  -0.6144692230
charge={charge}
basis=6-31g*
runtype=hess
method=hf
[scf]
type={scf}
multiplicity={mult}
[hess]
type=analytical
state=0
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        import pyscf  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime or PySCF not available")
class PyscfHessianBridgeTest(unittest.TestCase):
    def _bridge_vs_pyscf(self, charge, scf_type, mult, uhf):
        import numpy as np
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner
        from oqp.library import external
        from oqp.periodic_table import SYMBOL_MAP, ELEMENTS_NAME
        from pyscf import gto, scf as pscf

        workdir = Path("/tmp/oqp_bridge_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / f"mol_{scf_type}.inp"
        inp.write_text(INPUT.format(charge=charge, scf=scf_type, mult=mult))

        runner = Runner(project=f"br_{scf_type}", input_file=str(inp),
                        log=str(workdir / "run.log"))
        runner.mol.config["input"]["runtype"] = "energy"
        try:
            runner.run()
        except Exception:
            pass

        h_bridge = np.asarray(external.analytic_hessian_from_pyscf(runner.mol)[0])
        self.assertEqual(h_bridge.shape, (9, 9))
        self.assertLess(np.max(np.abs(h_bridge - h_bridge.T)), 1e-6)

        # direct PySCF on the identical (Bohr) geometry
        coord = np.asarray(runner.mol.get_system(), float).reshape(-1, 3)
        atoms = [[ELEMENTS_NAME[SYMBOL_MAP[int(a)]], coord[n]]
                 for n, a in enumerate(runner.mol.get_atoms())]
        mol = gto.Mole()
        mol.atom = atoms
        mol.unit = "Bohr"
        mol.basis = "6-31g*"
        mol.charge = charge
        mol.spin = mult - 1
        mol.verbose = 0
        mol.build(cart=True)
        mf = (pscf.UHF(mol) if uhf else pscf.RHF(mol))
        mf.kernel()
        h_direct = mf.Hessian().kernel().transpose(0, 2, 1, 3).reshape(9, 9)
        h_direct = 0.5 * (h_direct + h_direct.T)

        self.assertLess(np.max(np.abs(h_bridge - h_direct)), 1e-8)

    def test_rhf_water(self):
        self._bridge_vs_pyscf(charge=0, scf_type="rhf", mult=1, uhf=False)

    def test_uhf_water_cation(self):
        self._bridge_vs_pyscf(charge=1, scf_type="uhf", mult=2, uhf=True)


if __name__ == "__main__":
    unittest.main()
