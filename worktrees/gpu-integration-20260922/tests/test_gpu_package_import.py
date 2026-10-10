"""Importing planners must not start the historical DF preparation script."""
import os
from pathlib import Path
import subprocess
import sys
import unittest

class PackageImportTests(unittest.TestCase):
    def test_import_without_pyscf_or_cli_arguments(self):
        env = os.environ.copy()
        env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1] / 'python')
        result = subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['pyscf'] = None
sys.argv = ['import-only']
from openqp_gpu import build_df_tensor, solve_rhf
from openqp_gpu.gpu_workspace import GpuWorkspaceManager
assert callable(build_df_tensor) and callable(solve_rhf)
assert GpuWorkspaceManager() is not None
'''], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
