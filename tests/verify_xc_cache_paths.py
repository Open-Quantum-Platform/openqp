"""Compare caller-owned XC caches across response and gradient paths.

Run with the installed OpenQP Python and an explicit task-owned output directory.
"""
import json
import os
import re
import sys
import time
from pathlib import Path
import oqp
import numpy as np
from oqp.pyoqp import Runner

out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('OQP_XC_TIMING', '1')
water = 'O 0 0 -.041061554\nH -.533194329 .533194329 -.614469223\nH .533194329 -.533194329 -.614469223'
cases = [
    ('tda', 'rhf', 1, water, 'cc-pvdz'),
    ('rpa', 'rhf', 1, water, '6-31g*'),
    ('mrsf', 'rohf', 3, water, '6-31g*'),
    ('sf', 'rohf', 3, water, '6-31g*'),
    ('polar_rhf', 'rhf', 1, water, '6-31g*'),
    ('polar_uhf', 'uhf', 2, 'H 0 0 0', 'cc-pvdz'),
    ('polar_rohf', 'rohf', 3, water, '6-31g*'),
    ('adjoint_rohf', 'rohf', 3, water, '6-31g*'),
]
oqp.ffi.cdef('void cphf_static_polarizability(oqp_handle_t *, double *);', override=True)
rows = []
for name, scf, mult, geom, basis in cases:
    if len(sys.argv)>2 and name not in sys.argv[2:]:
        continue
    ref = None
    for mb in (0, 1, 256):
        os.environ['OQP_XC_RESPONSE_CACHE_MB'] = str(mb)
        folder = out / f'{name}-{mb}'
        folder.mkdir(exist_ok=True)
        os.chdir(folder)
        adjoint = name == 'adjoint_rohf'
        polar = name.startswith('polar') or adjoint
        cfg = {'input': {'system': geom, 'basis': basis, 'functional': os.environ.get('OQP_TEST_FUNCTIONAL', 'bhhlyp'),
                          'method': 'hf' if polar else 'tdhf',
                          'runtype': 'energy' if polar else 'grad', 'd4': 'False'},
               'scf': {'type': scf, 'multiplicity': str(mult), 'conv': '1e-11',
                       'maxit': '150', 'save_molden': 'False'},
               'guess': {'type': 'huckel', 'save_mol': 'False'},
               'symmetry': {'enabled': 'False'}}
        if adjoint:
            cfg['tdhf'] = {'type': 'mrsf', 'nstate': '6', 'maxit': '100'}
        if not polar:
            cfg['tdhf'] = {'type': name, 'multiplicity': '1', 'nstate': '6',
                           'conv': '1e-10', 'zvconv': '1e-10', 'maxit': '100',
                           'zv_warmstart': 'False'}
            if name == 'mrsf':
                cfg['tdhf']['z_solver'] = '1'
            cfg['properties'] = {'grad': '1'}
        (folder/'input.json').write_text(json.dumps(cfg, indent=2))
        runner = Runner(project=name, input_file=None, input_dict=cfg,
                        log=str(folder/'run.log'), silent=1, usempi=False)
        t = time.perf_counter()
        runner.run(test_mod=True)
        if polar:
            if adjoint:
                # Exercise the batched MINRES operator through its public
                # one-RHS adapter, with a reproducible nonzero source.
                data = runner.mol.data
                nbf = np.asarray(data['OQP::VEC_MO_A']).shape[0]
                na = int(np.asarray(data['nelec_A']).reshape(-1)[0])
                nb = int(np.asarray(data['nelec_B']).reshape(-1)[0])
                dim = nb*(nbf-nb)+(na-nb)*(nbf-na)
                data['OQP::nac_rohf_rhs'] = 0.01*np.sin(np.arange(1, dim+1))
                oqp.lib.mrsf_nac_rohf_zvector(data._data)
                values = np.asarray(data['OQP::nac_rohf_solution']).ravel().tolist()
            else:
                alpha = oqp.ffi.new('double[9]')
                oqp.lib.cphf_static_polarizability(runner.mol.data._data, alpha)
                values = list(alpha)
            # A direct CPHF call may open Fortran unit 6 after Runner closed it.
            # Keep that diagnostic file inside this case's private directory.
            diagnostic = ''.join(path.read_text() for path in
                                 (folder/'fort.6', folder/'run.log') if path.exists())
            assert 'did not converge' not in diagnostic, diagnostic
            if mb and os.environ['OQP_XC_TIMING'] == '1':
                hits = re.findall(r'kernel_hits=(\d+)', diagnostic)
                assert any(int(hit)>0 for hit in hits), (name, 'no XC cache reuse')
        else:
            values = list(runner.mol.energies) + np.asarray(runner.mol.grads[1]).ravel().tolist()
        elapsed = time.perf_counter() - t
        arr = np.asarray(values)
        assert np.isfinite(arr).all(), (name, mb, values)
        if ref is None:
            ref = arr.copy()
        delta = float(np.max(np.abs(ref-arr)))
        row = {'case': name, 'cache_mib': mb, 'seconds': elapsed,
               'max_difference': delta, 'values': values}
        rows.append(row)
        (out/'results.json').write_text(json.dumps(rows, indent=2))
        print(f'XC_PATH {name} mb={mb} seconds={elapsed:.5f} max_difference={delta:.3e}', flush=True)
        assert delta < 1e-9, row
print('XC_CACHE_PATHS_PASS', flush=True)
