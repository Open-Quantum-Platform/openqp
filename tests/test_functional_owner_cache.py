"""Same-process functional owners, settings and grid-size lifetime regression."""
import tempfile
import unittest
from pathlib import Path


class FunctionalOwnerCache(unittest.TestCase):
    def test_production_repeated_nmr_reuses_owner(self):
        import numpy as np
        import oqp
        from oqp.library.runfunc import compute_scf_prop
        from oqp.pyoqp import Runner

        revision = oqp.lib.functional_owner_cache_revision
        cleanup = oqp.lib.functional_owner_cache_clean
        with tempfile.TemporaryDirectory() as tmp:
            runner = Runner(project='nmr_owner', input_file=None, input_dict={
                'input': {'system': 'O .03 .11 0\nH .21 0 .97\nH .94 .08 -.24',
                          'basis': 'sto-3g', 'functional': 'pbe0', 'method': 'hf',
                          'runtype': 'energy', 'd4': 'False'},
                'scf': {'conv': '1e-10', 'save_molden': 'False'},
                'guess': {'save_mol': 'False'},
                'dftgrid': {'rad_npts': '24', 'ang_npts': '50'},
                'properties': {'scf_prop': 'nmr', 'nmr_gauge': 'giao'},
            }, log=str(Path(tmp)/'nmr.log'), silent=1, usempi=False)
            try:
                runner.run(test_mod=True)
                owner_revision = revision(runner.mol.data._data)
                before = {key: np.array(runner.mol.data[key], copy=True) for key in (
                    'OQP::nmr_shielding', 'OQP::nmr_pdens')}
                # Repeated property evaluation uses the ordinary production dispatch.
                # Copy native views before the next call replaces their backing data.
                for _ in range(2):
                    compute_scf_prop(runner.mol)
                    self.assertEqual(revision(runner.mol.data._data), owner_revision)
                    for key, expected in before.items():
                        current = np.array(runner.mol.data[key], copy=True)
                        self.assertTrue(np.isfinite(current).all(), key)
                        self.assertLessEqual(np.max(np.abs(current-expected)), 1e-10, key)
            finally:
                cleanup(runner.mol.data._data)

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
