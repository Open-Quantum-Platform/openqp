"""Under mpiexec only world rank 0 may append to the shared log.

``dump_log`` carries the ``mpi_dump`` guard, but the Runner writes the
performance-settings block directly.  Without its own guard every rank
appended the block, with its section heading, before rank 0 had written the
banner, which the CI step "Check MPI banner single-writer ordering" rejects.
"""
import os
import tempfile
import types
import unittest
from unittest import mock

try:
    import oqp.pyoqp as pyoqp_module
    _HAVE = True
except Exception:  # pragma: no cover - uncompiled backend
    _HAVE = False


@unittest.skipUnless(_HAVE, "compiled OpenQP backend unavailable")
class PerformanceSettingsSingleWriter(unittest.TestCase):
    def _write(self, world_rank):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "run.log")
            runner = types.SimpleNamespace(mol=types.SimpleNamespace(
                log=log, usempi=True, silent=True, perf_level=1,
                perf_report={"scf.xc_c2f": "off"}, perf_warns=[]))
            manager = types.SimpleNamespace(world_rank=world_rank)
            with mock.patch.object(pyoqp_module, "MPIManager", return_value=manager), \
                    mock.patch("oqp.utils.perf_levels.format_report",
                               return_value="Performance settings (perf = 1)"):
                pyoqp_module.Runner._log_perf_settings(runner)
            return open(log).read() if os.path.exists(log) else ""

    def test_only_world_rank_zero_writes_the_block(self):
        self.assertIn("PyOQP LOG | INPUT AND REFERENCE", self._write(0))
        self.assertEqual(self._write(1), "")


if __name__ == "__main__":
    unittest.main()
