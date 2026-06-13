"""Production guards for SOC-NAMD-QMMM option selection and hop bookkeeping.

These tests avoid importing the compiled OpenQP runtime. They pin the Python
dispatch and source-level invariants that make ``soc_basis=mch`` the preferred
SOC-QMMM production path.
"""

from __future__ import annotations

import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNFUNC = ROOT / "pyoqp" / "oqp" / "library" / "runfunc.py"
NAMD = ROOT / "pyoqp" / "oqp" / "library" / "namd.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"


def load_runfunc_with_namd_stubs():
    """Import runfunc.py with fake NAMD classes and no compiled oqp module."""
    calls = []

    class _Runner:
        def __init__(self, mol):
            self.mol = mol

        def run(self):
            calls.append(type(self).__name__)

    namd = types.ModuleType("oqp.library.namd")
    for name in (
        "NAMD",
        "NAMD_QMMM",
        "NAMD_SOC",
        "NAMD_SOC_QMMM",
        "NAMD_SOC_MCH",
        "NAMD_SOC_MCH_QMMM",
    ):
        setattr(namd, name, type(name, (_Runner,), {}))

    oqp = types.ModuleType("oqp")
    for name in ("electric_moments", "mulliken", "lowdin", "resp_charges"):
        setattr(oqp, name, lambda *_args, **_kwargs: None)

    single_point = types.ModuleType("oqp.library.single_point")
    noop = type("_Noop", (), {
        "__init__": lambda self, *_args, **_kwargs: None,
        "energy": lambda self: None,
        "gradient": lambda self: None,
        "hessian": lambda self: None,
        "compute": lambda self, *_args, **_kwargs: None,
        "overlap": lambda self: None,
        "nacme": lambda self: None,
        "nac": lambda self: None,
    })
    for name in ("SinglePoint", "Gradient", "Hessian", "LastStep", "BasisOverlap", "NACME", "NAC"):
        setattr(single_point, name, noop)

    modules = {
        "oqp": oqp,
        "oqp.library": types.ModuleType("oqp.library"),
        "oqp.library.namd": namd,
        "oqp.library.single_point": single_point,
        "oqp.library.libscipy": types.ModuleType("oqp.library.libscipy"),
        "oqp.library.libgeometric": types.ModuleType("oqp.library.libgeometric"),
        "oqp.library.liboqp": types.ModuleType("oqp.library.liboqp"),
    }
    for module_name, names in {
        "oqp.library.libscipy": ("StateSpecificOpt", "MECIOpt", "MECPOpt", "MEP", "QMMMOpt"),
        "oqp.library.libgeometric": (
            "GeometricIRCOpt",
            "GeometricMECIOpt",
            "GeometricMECPOpt",
            "GeometricNEBOpt",
            "GeometricOpt",
            "GeometricTSOpt",
            "GeometricQMMMOpt",
        ),
        "oqp.library.liboqp": (
            "OQPOpt",
            "OQPTSOpt",
            "OQPMECIOpt",
            "OQPMECPOpt",
            "OQPTCIOpt",
            "OQPNEBOpt",
            "OQPIRCOpt",
            "OQPMEPOpt",
        ),
    }.items():
        for name in names:
            setattr(modules[module_name], name, noop)

    saved = {name: sys.modules.get(name) for name in modules}
    sys.modules.update(modules)

    def cleanup():
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value

    spec = importlib.util.spec_from_file_location("runfunc_soc_namd_under_test", RUNFUNC)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, calls, cleanup


class DummyMol:
    def __init__(self, *, qmmm, soc, soc_basis="adiabatic"):
        self.config = {
            "input": {"qmmm_flag": qmmm},
            "md": {"soc": soc, "soc_basis": soc_basis},
        }


class SOCNAMDQMMMProductionTests(unittest.TestCase):
    def test_soc_basis_mch_dispatches_to_mch_qmmm_class(self):
        runfunc, calls, cleanup = load_runfunc_with_namd_stubs()
        try:
            runfunc.compute_namd(DummyMol(qmmm=True, soc=True, soc_basis="mch"))
        finally:
            cleanup()

        self.assertEqual(calls, ["NAMD_SOC_MCH_QMMM"])

    def test_adiabatic_soc_qmmm_dispatch_remains_available(self):
        runfunc, calls, cleanup = load_runfunc_with_namd_stubs()
        try:
            runfunc.compute_namd(DummyMol(qmmm=True, soc=True, soc_basis="adiabatic"))
        finally:
            cleanup()

        self.assertEqual(calls, ["NAMD_SOC_QMMM"])

    def test_soc_basis_schema_documents_mch_and_correction_switches(self):
        src = OQPDATA.read_text()

        self.assertIn("'soc_basis'", src)
        self.assertIn("'adiabatic' (SHARC) | 'mch' (spin-pure exact-gradient)", src)
        self.assertIn("'soc_du_dt_corr'", src)
        self.assertIn("'soc_tdc_grad_corr'", src)

    def test_soc_qmmm_hops_rescale_for_espf_energy_change(self):
        src = NAMD.read_text()

        self.assertGreaterEqual(src.count("de_espf = ((epot_old - epot) +"), 2)
        self.assertIn("eval_ha[self.active - 1] - eval_ha[active_old - 1]", src)
        self.assertIn("e_mch[self.active - 1] - e_mch[active_old - 1]", src)
        self.assertIn("self.v_all *= np.sqrt(max(0.0, 1.0 + de_espf / ekin_all))", src)

    def test_mch_qmmm_recomputes_exact_gradient_after_accepted_hop(self):
        src = NAMD.read_text()
        block = re.search(
            r"class NAMD_SOC_MCH_QMMM\(.*?dump_log\(mol, title='PyOQP: SOC-MCH-QMMM-NAMD trajectory complete'\)",
            src,
            re.S,
        )
        self.assertIsNotNone(block)
        body = block.group(0)

        self.assertIn("hopped = self._mch_propagate_and_hop(h_mch, e_mch)", body)
        self.assertIn("if hopped:", body)
        self.assertIn("g_qm, e_pure, mult, state, pchg = self._mch_exact_gradient_qmmm(self.active)", body)
        self.assertIn("f_all, epot = self._total_force_soc(potmm, g_qm, e_pure, pchg)", body)

    def test_overlap_tracking_and_rng_are_reproducible_contracts(self):
        src = NAMD.read_text()

        self.assertIn("self.rng = np.random.default_rng(self.seed)", src)
        self.assertIn("cfg['properties']['back_door'] = True", src)
        self.assertIn("NAMD_SOC._store_prev(self, self.r_all[self.qm_atoms].reshape((self.natom, 3)), u, eval_ha)", src)
        self.assertIn("BasisOverlap(mol).overlap()", src)
        self.assertIn("s_mch = NAMD_SOC._mch_overlap(self)", src)


if __name__ == "__main__":
    unittest.main()
