"""Same-process functional owners, settings and grid-size lifetime regression."""
import tempfile
import unittest
from pathlib import Path


class FunctionalOwnerCache(unittest.TestCase):
    def test_signed_xc_fock_and_owner_lifetime(self):
        import oqp
        from oqp.pyoqp import Runner

        # An absent compiled regression is a failure, including in CI.
        compare = oqp.lib.functional_owner_cache_compare
        for basis in ('6-31g', 'cc-pvdz'):
            with self.subTest(basis=basis), tempfile.TemporaryDirectory() as tmp:
                owners = []
                for index, system in enumerate((
                    'H 0 0 0\nH .13 .07 .74',
                    'O .03 .11 0\nH .21 0 .97\nH .94 .08 -.24',
                )):
                    cfg = {
                        'input': {'system': system, 'basis': basis, 'method': 'hf',
                                  'runtype': 'energy', 'd4': 'False'},
                        'scf': {'conv': '1e-10', 'save_molden': 'False'},
                        'guess': {'save_mol': 'False'},
                    }
                    runner = Runner(project=f'owner{index}', input_file=None,
                                    input_dict=cfg, log=str(Path(tmp)/f'{index}.log'),
                                    silent=1, usempi=False)
                    runner.run(test_mod=True)
                    owners.append(runner)
                error = oqp.ffi.new('double *')
                failed = oqp.ffi.new('int *')
                compare(owners[0].mol.data._data, owners[1].mol.data._data,
                        error, failed)
                print(f'FUNCTIONAL_OWNER basis={basis} error={error[0]:.3e} '
                      f'failed={failed[0]}', flush=True)
                self.assertEqual(failed[0], 0)
                self.assertLessEqual(error[0], 1e-11)
