"""Bounded XC cache: missing blocks, exact replay and content invalidation."""
import ctypes
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

class ResponseCacheTest(unittest.TestCase):
    def test_storage_with_bounds_checks(self):
        compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
        if compiler is None:
            self.skipTest('GNU Fortran required for bounds checks')
        root = Path(__file__).resolve().parents[1]
        source = (root/'tests/fortran/response_cache_selftest.F90').read_text()
        source = source.split('subroutine response_cache_compare', 1)[0]
        source += '''
program check_storage
  use iso_c_binding, only: c_int
  implicit none
  interface
    subroutine response_cache_storage_selftest(failed) bind(C)
      import c_int
      integer(c_int), intent(out) :: failed
    end subroutine
  end interface
  integer(c_int) :: failed
  call response_cache_storage_selftest(failed)
  if(failed/=0) error stop 'response cache storage failed'
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            check = Path(tmp)/'check.f90'
            check.write_text(source)
            exe = Path(tmp)/'check'
            subprocess.run([compiler, '-fdefault-integer-8', '-fopenmp',
                            '-fcheck=all', '-finit-real=snan', '-ffpe-trap=invalid,zero,overflow',
                            str(root/'source/precision.F90'),
                            str(root/'source/dftlib/dft_gridint_response_cache.F90'),
                            str(check), '-o', str(exe)], cwd=tmp, check=True,
                           capture_output=True, text=True)
            subprocess.run([str(exe)], check=True, capture_output=True, text=True)

    def test_budget_and_missing_blocks(self):
        import oqp
        from oqp.runtime import library_path
        lib = ctypes.CDLL(str(library_path(oqp.oqp_root, oqp.suffix)))
        func = lib.response_cache_storage_selftest
        func.argtypes = [ctypes.POINTER(ctypes.c_int)]
        func.restype = None
        failed = ctypes.c_int()
        func(ctypes.byref(failed))
        self.assertEqual(failed.value, 0)

    def test_response_replay_and_invalidation(self):
        import oqp
        from oqp.pyoqp import Runner
        oqp.ffi.cdef('void response_cache_compare(oqp_handle_t *, double *, int *);', override=True)
        for basis, functional in (('6-31g*', 'bhhlyp'), ('cc-pvdz', 'bhhlyp'),
                                  ('6-31g*', 'slater')):
            with self.subTest(basis=basis, functional=functional), tempfile.TemporaryDirectory() as tmp:
                cfg = {'input': {'system': 'O 0 0 0\nH 0 0 .97\nH .94 0 -.24',
                       'basis': basis, 'method': 'hf', 'functional': functional,
                       'runtype': 'energy', 'd4': 'False'},
                       'scf': {'conv': '1e-10', 'save_molden': 'False'},
                       'guess': {'save_mol': 'False'}}
                r = Runner(project='cache', input_file=None, input_dict=cfg,
                           log=str(Path(tmp)/'run.log'), silent=1, usempi=False)
                r.run(test_mod=True)
                error = oqp.ffi.new('double *'); failed = oqp.ffi.new('int *')
                oqp.lib.response_cache_compare(r.mol.data._data, error, failed)
                print(f'XC_CACHE_COMPARE basis={basis} functional={functional} error={error[0]:.3e} failed={failed[0]}', flush=True)
                self.assertEqual(failed[0], 0)
                self.assertLess(error[0], 1e-10)

    def test_gradient_ao_replay(self):
        import oqp
        from oqp.pyoqp import Runner
        oqp.ffi.cdef('void gradient_cache_compare(oqp_handle_t *, double *, int *);', override=True)
        for basis, functional in (('6-31g*', 'bhhlyp'), ('cc-pvdz', 'slater')):
            with self.subTest(basis=basis, functional=functional), tempfile.TemporaryDirectory() as tmp:
                cfg = {'input': {'system': 'O 0 0 0\nH 0 0 .97\nH .94 0 -.24',
                       'basis': basis, 'method': 'hf', 'functional': functional,
                       'runtype': 'energy', 'd4': 'False'},
                       'scf': {'conv': '1e-10', 'save_molden': 'False'},
                       'guess': {'save_mol': 'False'}}
                r = Runner(project='gradient_cache', input_file=None, input_dict=cfg,
                           log=str(Path(tmp)/'run.log'), silent=1, usempi=False)
                r.run(test_mod=True)
                error = oqp.ffi.new('double *'); failed = oqp.ffi.new('int *')
                oqp.lib.gradient_cache_compare(r.mol.data._data, error, failed)
                print(f'GRADIENT_CACHE_COMPARE basis={basis} functional={functional} error={error[0]:.3e} failed={failed[0]}', flush=True)
                self.assertEqual(failed[0], 0)
                self.assertLess(error[0], 1e-10)
