"""ROHF status and NMR claim-boundary guards.

ROHF ground-state SCF support is distinct from ROHF NMR response support.  These
checks prevent requested ROHF NMR cases from being silently routed through RHF,
UHF, or a GIAO-to-CGO fallback while the NMR response remains unvalidated.
"""
from __future__ import annotations

import importlib.util
import os
import platform
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNFUNC = ROOT / "pyoqp" / "oqp" / "library" / "runfunc.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_CHECKER = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"
SCF = ROOT / "source" / "scf.F90"
NMR = ROOT / "source" / "modules" / "nmr_shielding.F90"
STATUS = ROOT / "source" / "modules" / "NMR_SHIELDING_STATUS.md"

ROHF_NMR_INPUT = """\
[input]
system=
  H 0.0 0.0 0.0
charge=0
runtype=energy
basis=sto-3g
method=hf

[guess]
type=huckel

[scf]
type=rohf
multiplicity=2
maxit=5

[properties]
scf_prop=nmr
nmr_gauge=cgo
"""


class DummyMol:
    def __init__(self, scf_type="rohf", nmr_gauge="cgo"):
        self.config = {
            "scf": {"type": scf_type},
            "properties": {"scf_prop": ["nmr"], "nmr_gauge": nmr_gauge},
        }


def load_runfunc_with_stubs():
    """Import runfunc.py without requiring the compiled oqp extension."""
    calls = {"nmr": 0}

    oqp = types.ModuleType("oqp")

    def nmr_shielding(_mol):
        calls["nmr"] += 1

    setattr(oqp, "nmr_shielding", nmr_shielding)
    setattr(oqp, "electric_moments", lambda _mol: None)
    setattr(oqp, "mulliken", lambda _mol: None)
    setattr(oqp, "lowdin", lambda _mol: None)
    setattr(oqp, "resp_charges", lambda _mol: None)

    library = types.ModuleType("oqp.library")
    dftbplus = types.ModuleType("oqp.library.dftbplus")
    setattr(dftbplus, "optimize_openqp_molecule", lambda *args, **kwargs: None)
    setattr(dftbplus, "run_openqp_molecule", lambda *args, **kwargs: None)

    single_point = types.ModuleType("oqp.library.single_point")
    class _Noop:
        def __init__(self, *args, **kwargs):
            pass
        def energy(self):
            pass
        def compute(self, *args, **kwargs):
            pass
        def gradient(self):
            pass
        def hessian(self):
            pass
        def reference(self):
            return None
        def excitation(self, *args, **kwargs):
            pass
        def overlap(self):
            pass
        def nacme(self):
            pass
        def nac(self):
            pass
    for name in ("SinglePoint", "Gradient", "Hessian", "LastStep", "BasisOverlap", "NACME", "NAC"):
        setattr(single_point, name, _Noop)

    libscipy = types.ModuleType("oqp.library.libscipy")
    for name in ("StateSpecificOpt", "MECIOpt", "MECPOpt", "MEP"):
        setattr(libscipy, name, _Noop)

    libdlfind = types.ModuleType("oqp.library.libdlfind")
    for name in ("DLFindMin", "DLFindTS", "DLFindMECI"):
        setattr(libdlfind, name, _Noop)

    libgeometric = types.ModuleType("oqp.library.libgeometric")
    for name in (
        "GeometricIRCOpt",
        "GeometricMECIOpt",
        "GeometricMECPOpt",
        "GeometricNEBOpt",
        "GeometricOpt",
        "GeometricTSOpt",
    ):
        setattr(libgeometric, name, _Noop)

    saved = {name: sys.modules.get(name) for name in (
        "oqp",
        "oqp.library",
        "oqp.library.dftbplus",
        "oqp.library.single_point",
        "oqp.library.libscipy",
        "oqp.library.libdlfind",
        "oqp.library.libgeometric",
    )}
    sys.modules.update({
        "oqp": oqp,
        "oqp.library": library,
        "oqp.library.dftbplus": dftbplus,
        "oqp.library.single_point": single_point,
        "oqp.library.libscipy": libscipy,
        "oqp.library.libdlfind": libdlfind,
        "oqp.library.libgeometric": libgeometric,
    })
    try:
        spec = importlib.util.spec_from_file_location("runfunc_under_test", RUNFUNC)
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module, calls
    finally:
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value


def load_input_checker_with_minimal_stubs():
    """Load the real input checker with only its lightweight oqp stubs."""
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    setattr(mpi_utils, "MPIManager", MPIManager)
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    spec = importlib.util.spec_from_file_location(
        "input_checker_rohf_nmr_under_test", INPUT_CHECKER
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ROHFStatusAndInterfaceTests(unittest.TestCase):
    def test_rohf_is_recognized_in_input_interface(self):
        oqpdata = OQPDATA.read_text()
        checker = INPUT_CHECKER.read_text()
        self.assertIn('"rohf": 3', oqpdata)
        self.assertIn("SCF_TYPES", checker)
        self.assertIn('"rohf"', checker)
        self.assertIn("SF/MRSF requires an ROHF reference", checker)

    def test_real_input_checker_leaves_rohf_nmr_to_runtime_guard(self):
        input_checker = load_input_checker_with_minimal_stubs()
        config = {
            "input": {
                "method": "hf",
                "runtype": "energy",
                "basis": "sto-3g",
                "system": "\nO 0 0 0\nH 0 0 1\nH 0 1 0",
            },
            "scf": {"type": "rohf", "multiplicity": 3},
            "properties": {"scf_prop": ["nmr"], "nmr_gauge": "cgo"},
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())
        self.assertEqual(config["scf"]["type"], "rohf")

    def test_rohf_has_native_scf_driver_scaffold(self):
        scf = SCF.read_text()
        for required in (
            "scf_rohf",
            "get_scf_name",
            "form_rohf_fock",
            "Guest-Saunders ROHF Fock",
            "mo_b = mo_a",
            "dmat_a = pdmat(:,1) - pdmat(:,2)",
        ):
            self.assertIn(required, scf)

    def test_rohf_nmr_raises_notimplemented_without_rhf_or_uhf_fallback(self):
        runfunc, calls = load_runfunc_with_stubs()
        mol = DummyMol(scf_type="rohf", nmr_gauge="cgo")
        before = mol.config["scf"]["type"]
        with self.assertRaisesRegex(NotImplementedError, "ROHF NMR shielding requested") as ctx:
            runfunc.compute_scf_prop(mol)
        self.assertIn("ROHF SCF/reference support exists", str(ctx.exception))
        self.assertIn("ROHF NMR response support is separate", str(ctx.exception))
        self.assertEqual(mol.config["scf"]["type"], before)
        self.assertEqual(calls["nmr"], 0, "ROHF NMR must not dispatch to the CGO routine")

    def test_live_rohf_nmr_request_fails_after_rohf_scf_without_cgo_dispatch(self):
        suffix = {
            "Darwin": "dylib",
            "Linux": "so",
            "Windows": "dll",
        }.get(platform.system(), "so")
        runtime_root = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
        if not (runtime_root / "lib" / f"liboqp.{suffix}").exists():
            self.skipTest("OpenQP shared library not built for this platform; skipping live ROHF/NMR guard")
        if importlib.util.find_spec("basis_set_exchange") is None:
            self.skipTest("basis_set_exchange missing for live PyOQP driver")

        with tempfile.TemporaryDirectory() as workdir:
            inp = Path(workdir) / "h_rohf_nmr.inp"
            inp.write_text(ROHF_NMR_INPUT)
            env = dict(os.environ)
            env["OPENQP_ROOT"] = str(runtime_root)
            # Use this checkout's Python wrapper while allowing it to load the
            # selected compiled OpenQP runtime from OPENQP_ROOT.  This keeps the
            # ROHF/NMR guard test from accidentally exercising an older
            # pip-installed wrapper that may not know properties.nmr_gauge.
            root_parent = str(Path(runtime_root).parent)
            source_parent = str(ROOT / "pyoqp")
            env["PYTHONPATH"] = os.pathsep.join(
                item for item in (source_parent, root_parent, env.get("PYTHONPATH", "")) if item
            )
            proc = subprocess.run(
                [sys.executable, "-m", "oqp.pyoqp", str(inp)],
                cwd=workdir,
                env=env,
                capture_output=True,
                text=True,
                timeout=60,
            )
            log = inp.with_suffix(".log").read_text() if inp.with_suffix(".log").exists() else ""

        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("ROHF NMR shielding requested", proc.stderr)
        self.assertIn("ROHF SCF/reference support exists", proc.stderr)
        self.assertIn("ROHF NMR response support is separate", proc.stderr)
        self.assertIn("SCF type = ROHF", log)
        self.assertIn("Final ROHF energy", log)
        self.assertNotIn("Final RHF energy", log)
        self.assertNotIn("Final UHF energy", log)
        self.assertNotIn("NMR shielding", log)

    def test_rohf_nmr_giao_is_not_routed_to_cgo(self):
        runfunc, calls = load_runfunc_with_stubs()
        mol = DummyMol(scf_type="rohf", nmr_gauge="giao")
        with self.assertRaises(NotImplementedError):
            runfunc.compute_scf_prop(mol)
        self.assertEqual(calls["nmr"], 0, "ROHF-GIAO must not use CGO fallback")

    def test_nmr_fortran_still_declares_closed_shell_only(self):
        nmr = NMR.read_text()
        self.assertIn("NMR shielding currently supports closed-shell", nmr)
        self.assertIn("infos%control%scftype == 2 .or. infos%control%scftype == 3", nmr)

    def test_status_ledger_records_conservative_rohf_status(self):
        status = STATUS.read_text()
        for required in (
            "ROHF implementation / audit checkpoint",
            "ROHF ground-state SCF support and ROHF NMR response support are separate tasks",
            "ROHF-CGO NMR response | Not implemented / not validated",
            "ROHF-GIAO NMR response | Not implemented / not validated",
            "ROHF-GIAO shielding | Not validated",
            "RHF/UHF fallback for requested ROHF | Not allowed",
            "ROHF SCF alone is insufficient for open-shell NMR",
            "ABI/build caution",
        ):
            self.assertIn(required, status)

    def test_status_ledger_does_not_advertise_rohf_nmr_as_validated(self):
        status = STATUS.read_text().lower()
        forbidden_claims = (
            "rohf-cgo nmr response | validated",
            "rohf-giao nmr response | validated",
            "rohf-giao shielding | validated",
            "rohf nmr shielding is validated",
        )
        for claim in forbidden_claims:
            self.assertNotIn(claim, status)


if __name__ == "__main__":
    unittest.main()
