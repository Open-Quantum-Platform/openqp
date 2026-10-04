"""Run finite differences and repeated response solves in a task-owned directory."""
import json
import sys
import time
from pathlib import Path
import oqp  # establish the native runtime before numpy
import numpy as np
from oqp.pyoqp import Runner

out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
water = [[8, 0., 0., -.041061554], [1, -.533194329, .533194329, -.614469223],
         [1, .533194329, -.533194329, -.614469223]]
records = []

def run(name, geom, *, basis='cc-pvdz', scf='rhf', mult=1, excited=False,
        solver=1, runtype='grad', functional=None):
    cfg = {'input': {'system': '\n'.join(str(int(r[0]))+' '+' '.join(map(str, r[1:])) for r in geom),
                     'basis': basis, 'method': 'tdhf' if excited else 'hf',
                     'runtype': runtype, 'd4': 'False'},
           'scf': {'type': scf, 'multiplicity': str(mult), 'conv': '1e-11', 'maxit': '150'},
           'guess': {'type': 'huckel', 'save_mol': 'False'},
           'properties': {'grad': '1' if excited else '0'},
           'symmetry': {'enabled': 'false'}}
    if functional:
        cfg['input']['functional'] = functional
    if excited:
        cfg['tdhf'] = {'type': 'mrsf', 'multiplicity': '1', 'nstate': '6',
                       'conv': '1e-10', 'zvconv': '1e-10', 'z_solver': str(solver),
                       'zv_warmstart': 'False', 'maxit': '100'}
    folder = out / name
    folder.mkdir(exist_ok=True)
    (folder / 'input.json').write_text(json.dumps(cfg, indent=2))
    start = time.perf_counter()
    runner = Runner(project=name, input_file=None, input_dict=cfg,
                    log=str(folder / 'run.log'), silent=1, usempi=False)
    runner.run(test_mod=True)
    energy = float(runner.mol.energies[1 if excited else 0])
    grad = None if runtype == 'energy' else np.asarray(runner.mol.grads[1 if excited else 0])
    record = {'name': name, 'energy': energy, 'seconds': time.perf_counter()-start,
              'gradient': None if grad is None else grad.tolist()}
    records.append(record)
    (out / 'results.json').write_text(json.dumps(records, indent=2))
    print(json.dumps(record), flush=True)
    return energy, grad

for name, geom, kwargs, atom, axis in [
    ('rhf_f', water, dict(basis='cc-pvtz'), 1, 0),
    ('uhf_d', [[8,0.,0.,0.], [1,0.,0.,.97]], dict(basis='cc-pvdz',scf='uhf',mult=2), 1, 2),
]:
    e, g = run(name, geom, **kwargs)
    h = 3e-4
    plus, minus = np.array(geom), np.array(geom)
    plus[atom,axis+1] += h; minus[atom,axis+1] -= h
    ep, _ = run(name+'_plus', plus.tolist(), runtype='energy', **kwargs)
    em, _ = run(name+'_minus', minus.tolist(), runtype='energy', **kwargs)
    fd = (ep-em)/(2*h)*0.529177210903
    err = abs(fd-g[atom,axis])
    print(f'GRADIENT_FD {name} error={err:.3e}', flush=True)
    assert err < 1e-6, (name, fd, g[atom,axis])

reference = None
for solver in (1, 2, 0, 1):
    name = f'mrsf_z{solver}_call{len(records)}'
    e, g = run(name, water, basis='6-31g*', scf='rohf', mult=3,
               excited=True, solver=solver, functional='bhhlyp')
    if reference is None:
        reference = g.copy()
    err = float(np.max(np.abs(g-reference)))
    print(f'ZVECTOR_COMPARE solver={solver} max_gradient_difference={err:.3e}', flush=True)
    assert err < 1e-6, (solver, err)
print('ERI_RESPONSE_MEMORY_VERIFICATION_PASS', flush=True)
