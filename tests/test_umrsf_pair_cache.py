"""Signed production-density pin for the quadratic UMRSF spin-pair matrix."""
import ctypes
import os
import json
import sys
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch
import unittest


class UmrsfPairCache(unittest.TestCase):
    def test_signed_sigma_and_repeated_lifetime(self):
        import oqp.runtime as runtime
        root, suffix = runtime.resolve_oqp_root()
        lib = ctypes.CDLL(runtime.library_path(root, suffix))
        check = lib.umrsf_pair_selftest
        check.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int)]
        check.restype = None
        for option in ("1", "0", None):
            with self.subTest(option=option), patch.dict(os.environ):
                os.environ.pop("OQP_UMRSF_PAIR_CACHE", None)
                if option is not None:
                    os.environ["OQP_UMRSF_PAIR_CACHE"] = option
                error, count, failed = ctypes.c_double(), ctypes.c_int(), ctypes.c_int()
                check(ctypes.byref(error), ctypes.byref(count), ctypes.byref(failed))
                self.assertEqual(count.value, 160)
                self.assertEqual(failed.value, 0)
                self.assertLessEqual(error.value, 2e-10)
                print(f"UMRSF_PAIR_PIN cases={count.value} max_abs={error.value:.16g}")


class UmrsfPairLive(unittest.TestCase):
    def test_molecular_dispatch_and_spectrum(self):
        import numpy as np
        source = Path(__file__).resolve().parents[1] / "examples/UMRSF-TDDFT/H2O_UMRSF_PAIR_CACHE.inp"
        cases = (("singlet", 1, "6-31g", ""), ("triplet", 3, "6-31g", ""),
                 ("spc", 1, "6-31g", "spc_coco=0.2\nspc_ovov=0.7\nspc_coov=0.35\n"),
                 ("spherical", 1, "cc-pvdz", ""))
        for name, spin, basis, extra in cases:
            results = []
            for option in ((None, "0", "1") if name == "singlet" else ("0", "1")):
                with self.subTest(case=name, option=option), tempfile.TemporaryDirectory(prefix="umrsf_pair_") as wd:
                    text = source.read_text().replace("basis=6-31g", "basis=" + basis)
                    text = text.replace("multiplicity=1", f"multiplicity={spin}") + extra
                    Path(wd, "run.inp").write_text(text)
                    env = dict(os.environ, OQP_ROUTEC_SIG="0")
                    env.pop("OQP_UMRSF_PAIR_CACHE", None)
                    if option is not None:
                        env["OQP_UMRSF_PAIR_CACHE"] = option
                    proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--live"],
                                          cwd=wd, env=env, text=True, capture_output=True, timeout=180)
                    self.assertEqual(proc.returncode, 0, proc.stdout[-3000:] + proc.stderr[-3000:])
                    result = json.loads(Path(wd, "result.json").read_text())
                    results.append(result)
                    log = Path(wd, "run.log").read_text()
                    self.assertIn("energies converged in", log)
                    self.assertEqual("1 dynamic exchange; quadratic spin-pair matrix" in log, option == "1")
            for result in results[1:]:
                np.testing.assert_allclose(result, results[0], atol=1e-8, rtol=0)


if __name__ == '__main__':
    if sys.argv[1:] == ["--live"]:
        from oqp.pyoqp import Runner
        runner = Runner(project="pair", input_file="run.inp", log="run.log", silent=1, usempi=False)
        runner.run()
        Path("result.json").write_text(json.dumps(runner.mol.data["OQP::td_energies"].tolist()))
    else:
        unittest.main()
