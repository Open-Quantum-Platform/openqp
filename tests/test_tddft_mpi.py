"""TDDFT/TDHF RPA and TDA excitation energies under MPI equal the serial ones.

The XC kernel contraction distributes grid slices over the MPI ranks
(dft_gridint: ``mod(iChunk/chunkSize, pe%size) == pe%rank``) and reduces the
response Fock matrices afterwards.  A buffer missing from that reduction, or a
rank-local quantity used before it, changes the excitation energies with the
rank count.  Covered: the closed-shell singlet and triplet branches, the
CONTIGUOUS MO-pointer consumer (GGA and the meta-GGA ``moG1`` branch), the
spin-antisymmetric triplet contraction, the CAM two-pass integrals, and the
UMRSF quintet operator (HF, global hybrid and CAM) in the complete response space.

Runs only with mpi4py and an MPI launcher.  ``OQP_REQUIRE_MPI_TESTS=1`` (set
in the CI MPI leg) turns the skip into a failure, so a broken MPI environment
cannot pass silently.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

H2O_C1 = ("   O 0.010000 -0.020000 0.117300\n   H 0.080000 0.770200 -0.449200\n"
          "   H -0.050000 -0.731200 -0.509200\n")

CASES = [  # functional, [tdhf] type, multiplicity
    ("b3lypv5", "rpa", 1),
    ("b3lypv5", "rpa", 3),
    ("b3lypv5", "tda", 3),
    ("m06-2x", "rpa", 3),
    ("cam-b3lyp", "rpa", 3),
]
RANKS = (2, 4)
REQUIRED = os.environ.get("OQP_REQUIRE_MPI_TESTS", "") == "1"


def _launcher():
    try:
        import mpi4py  # noqa: F401
        import oqp  # noqa: F401
    except (ImportError, OSError, RuntimeError):
        # not installed, or liboqp/libmpi not loadable (oqp.runtime raises
        # RuntimeError when it cannot load the native library)
        return None
    return shutil.which("mpiexec") or shutil.which("mpirun")


def _deck(functional, td_type, mult):
    return (f"[input]\nsystem=\n{H2O_C1}charge=0\nruntype=energy\nbasis=6-31g*\n"
            f"method=tdhf\nfunctional={functional}\n\n[guess]\ntype=huckel\n\n"
            "[scf]\ntype=rhf\nmultiplicity=1\nmaxit=100\nconv=1e-10\n\n"
            f"[tdhf]\ntype={td_type}\nnstate=4\nmaxit=100\nconv=1e-10\nmultiplicity={mult}\n\n"
            "[dftgrid]\nrad_npts=96\nang_npts=302\n")


def _run(text, nrank, launcher):
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        cmd = [sys.executable, "-m", "oqp.pyoqp", "job.inp"]
        if nrank > 1:
            cmd = [launcher, "-n", str(nrank), "--oversubscribe", "--bind-to", "none"] + cmd
        proc = subprocess.run(cmd, cwd=tmp, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=900)
        log_path = Path(tmp) / "job.log"
        log = log_path.read_text(errors="ignore") if log_path.exists() else ""
        summary = log.split("Summary of the TD-DFT calculation")[-1]
        ev = [float(x) for x in re.findall(r"^\s*0 ->\s*\d+\s+(-?[\d.]+)", summary, re.M)]
        return proc.returncode, ev, proc.stdout.decode(errors="ignore")[-1500:] + log[-2000:]


def _quintet_snapshot(functional, nrank, launcher):
    # Use full-precision eigenpairs, not rounded energies from the summary.
    text = ("[input]\nsystem=\n C 0 0 0.104\n H 0 0.862 -0.468\n"
            " H 0 -0.862 -0.468\nbasis=sto-3g\nruntype=energy\nmethod=tdhf\n"
            f"functional={functional}\n[guess]\ntype=huckel\n"
            "[scf]\ntype=uhf\nmultiplicity=3\nconv=1e-10\nmaxit=200\n"
            "[tdhf]\ntype=umrsf\nmultiplicity=5\nnstate=6\nconv=1e-10\nmaxit=100\n")
    script = """
import oqp
import json
import numpy as np
from pathlib import Path
from oqp.pyoqp import Runner
from oqp.utils.mpi_utils import MPIManager
r = Runner(input_file='job.inp', log='job.log')
r.run()
d = r.mol.data
assert r.mol.mol_energy.Davidson_converged
assert d['nelec_A'] == 5 and d['nelec_B'] == 3
payload = dict(energies=np.array(d['OQP::td_energies']).ravel().tolist(),
               vectors=np.array(d['OQP::td_bvec_mo']).reshape(6, 6).T.tolist(),
               a=np.array(d['OQP::VEC_MO_A']).T.tolist(),
               b=np.array(d['OQP::VEC_MO_B']).T.tolist(),
               overlap=np.array(d['OQP::SM']).ravel().tolist(),
               hamilton=int(d['hamilton']), cam_flag=bool(d['dft'].cam_flag),
               cam=[float(d[k]) for k in ('cam_alpha', 'cam_beta', 'cam_mu')],
               ranks=MPIManager().size)
manager = MPIManager()
all_ranks = manager.comm.allgather(payload) if manager.use_mpi else [payload]
# Native outputs must be valid on every rank, including non-root ranks.
for item in all_ranks:
    for key in ('energies', 'vectors', 'a', 'b', 'overlap'):
        value = np.asarray(item[key])
        assert np.isfinite(value).all()
        np.testing.assert_allclose(value, all_ranks[0][key], rtol=0, atol=1e-9,
                                   equal_nan=False)
if manager.rank == 0:
    Path('quintet.json').write_text(json.dumps(payload, allow_nan=False))
manager.barrier()
"""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root/'job.inp').write_text(text)
        cmd = [sys.executable, '-c', script]
        if nrank > 1:
            cmd = [launcher, '-n', str(nrank), '--oversubscribe', '--bind-to', 'none'] + cmd
        proc = subprocess.run(cmd, cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS='1'),
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
        assert proc.returncode == 0, proc.stdout.decode(errors='replace')[-4000:]
        log = (root/'job.log').read_text()
        if nrank > 1:
            # This banner comes from native pe%size/use_mpi, not mpi4py alone:
            # launching multiple copies of a serial library must not pass.
            assert re.search(rf'MPI Processors:\s+{nrank}\b', log), log[-4000:]
        result = json.loads((root/'quintet.json').read_text())
        assert result['ranks'] == nrank
        return result


class TestTddftUnderMpi(unittest.TestCase):
    def setUp(self):
        self.launcher = _launcher()
        if self.launcher is None:
            if REQUIRED:
                self.fail("OQP_REQUIRE_MPI_TESTS=1 but mpi4py, oqp or an MPI launcher is missing")
            self.skipTest("mpi4py, the compiled OpenQP runtime and an MPI launcher are needed")

    def test_excitation_energies_independent_of_rank_count(self):
        for functional, td_type, mult in CASES:
            text = _deck(functional, td_type, mult)
            rc0, e0, d0 = _run(text, 1, self.launcher)
            self.assertEqual(rc0, 0, d0)
            self.assertEqual(len(e0), 4, d0)
            for nrank in RANKS:
                with self.subTest(functional=functional, td=td_type, mult=mult, ranks=nrank):
                    rc, e, d = _run(text, nrank, self.launcher)
                    self.assertEqual(rc, 0, d)
                    self.assertEqual(len(e), 4, d)
                    for a, b in zip(e0, e):
                        self.assertAlmostEqual(a, b, delta=1e-6)

    def test_umrsf_quintet_operator_independent_of_rank_count(self):
        import numpy as np
        for functional in ('', 'bhhlyp', 'cam-b3lyp'):
            matrices = []
            energies = []
            orbitals = []
            for nrank in (1, *RANKS):
                with self.subTest(functional=functional or 'HF', ranks=nrank):
                    result = _quintet_snapshot(functional, nrank, self.launcher)
                    if functional == 'cam-b3lyp':
                        self.assertEqual(result['hamilton'], 20)
                        self.assertTrue(result['cam_flag'])
                        np.testing.assert_allclose(result['cam'], [0.19, 0.46, 0.33],
                                                   rtol=0, atol=1e-12, equal_nan=False)
                    e = np.asarray(result['energies'])
                    x = np.asarray(result['vectors'])
                    self.assertEqual(e.shape, (6,))
                    self.assertEqual(x.shape, (6, 6))
                    self.assertTrue(np.isfinite(e).all() and np.isfinite(x).all())
                    np.testing.assert_allclose(x.T@x, np.eye(6), rtol=0, atol=1e-10,
                                               equal_nan=False)
                    a = np.asarray(result['a'])
                    b = np.asarray(result['b'])
                    overlap = np.zeros((7, 7))
                    overlap[np.tril_indices(7)] = result['overlap']
                    overlap += np.tril(overlap, -1).T
                    self.assertTrue(all(np.isfinite(v).all() for v in (a, b, overlap)))
                    orbitals.append((a, b))
                    a0, b0 = orbitals[0]
                    uc = b0[:, :3].T@overlap@b[:, :3]
                    uv = a0[:, 5:].T@overlap@a[:, 5:]
                    np.testing.assert_allclose(uc.T@uc, np.eye(3), rtol=0, atol=1e-9,
                                               equal_nan=False)
                    np.testing.assert_allclose(uv.T@uv, np.eye(2), rtol=0, atol=1e-9,
                                               equal_nan=False)
                    transform = np.kron(uv, uc)
                    # Express each response matrix in the serial MO gauge.
                    matrix = transform@((x*e)@x.T)@transform.T
                    self.assertTrue(np.isfinite(matrix).all())
                    energies.append(e)
                    matrices.append(matrix)
                    # Signed matrix elements retain the operator's signs while
                    # removing the arbitrary phase of individual Ritz vectors.
                    np.testing.assert_allclose(e, energies[0], rtol=0, atol=1e-9,
                                               equal_nan=False)
                    np.testing.assert_allclose(matrix, matrices[0], rtol=0, atol=1e-9,
                                               equal_nan=False)


class TestLauncherProbe(unittest.TestCase):
    """The probe treats an unloadable runtime as "MPI unavailable"; only
    OQP_REQUIRE_MPI_TESTS turns that into a failure (setUp above)."""

    def test_runtime_loader_error_means_unavailable(self):
        import builtins
        from unittest import mock
        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "mpi4py":
                return mock.MagicMock()
            if name == "oqp":
                raise RuntimeError("cannot load liboqp")
            return real_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=fake_import):
            self.assertIsNone(_launcher())


class TestTripletGradientUnderMpi(unittest.TestCase):
    """The triplet gradient (utddft_fxc, utddft_gxc, utddft_xc_gradient with
    the f_xc moving-grid term) at 2 ranks equals the serial gradient."""

    def setUp(self):
        self.launcher = _launcher()
        if self.launcher is None:
            if REQUIRED:
                self.fail("OQP_REQUIRE_MPI_TESTS=1 but mpi4py, oqp or an MPI launcher is missing")
            self.skipTest("mpi4py, the compiled OpenQP runtime and an MPI launcher are needed")

    def _grad(self, td_type, nrank):
        text = _deck("b3lypv5", td_type, 3)
        text = text.replace("runtype=energy", "runtype=grad").replace(
            "[guess]\ntype=huckel\n", "[guess]\ntype=huckel\nsave_mol=true\n")
        text += "\n[properties]\ngrad=1\n"
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "job.inp").write_text(text)
            env = dict(os.environ, OMP_NUM_THREADS="1")
            cmd = [sys.executable, "-m", "oqp.pyoqp", "job.inp"]
            if nrank > 1:
                cmd = [self.launcher, "-n", str(nrank), "--oversubscribe", "--bind-to", "none"] + cmd
            proc = subprocess.run(cmd, cwd=tmp, env=env, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, timeout=900)
            res = Path(tmp) / "job.json"
            data = json.loads(res.read_text()) if res.exists() else {}
            return proc.returncode, data.get("grad"), proc.stdout.decode(errors="ignore")[-1500:]

    def test_triplet_gradient_independent_of_rank_count(self):
        for td_type in ("tda", "rpa"):
            with self.subTest(td=td_type):
                rc0, g0, d0 = self._grad(td_type, 1)
                rc2, g2, d2 = self._grad(td_type, 2)
                self.assertEqual(rc0, 0, d0)
                self.assertEqual(rc2, 0, d2)
                self.assertIsNotNone(g0, d0)
                self.assertIsNotNone(g2, d2)
                for a, b in zip(g0, g2):
                    self.assertAlmostEqual(a, b, delta=1e-7)


if __name__ == "__main__":
    unittest.main()
