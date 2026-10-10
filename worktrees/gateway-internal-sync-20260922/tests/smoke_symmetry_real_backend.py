"""End-to-end smoke test of symmetry labeling against a real OpenQP backend.

Not collected by pytest (no test_ prefix). Run with a Python environment
that has a built `oqp` package installed:

    <venv>/bin/python3 tests/smoke_symmetry_real_backend.py [workdir]

The script grafts this branch's symmetry modules and Molecule labeling
methods onto the installed oqp package, runs real RHF and MRSF water
calculations, and checks the labels. This validates the real data paths
(get_basis layout, triangular OQP::SM, VEC_MO orientation, td_bvec_mo
layout) that the stub-based unit tests cannot.
"""

import importlib.util
import os
import sys
from pathlib import Path

BRANCH_ROOT = Path(__file__).resolve().parents[1]

RHF_INPUT = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
method=hf

[guess]
type=huckel

[scf]
multiplicity=1
type=rhf
"""

MRSF_INPUT = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf

[tdhf]
type=mrsf
nstate=3
"""

GRAFTED_METHODS = [
    '_parse_bool_like',
    '_parse_enabled_mode',
    'initialize_symmetry_metadata',
    '_detect_symmetry_metadata',
    '_symmetry_labeling_inputs',
    '_mo_coefficients',
    'label_molecular_orbitals',
    '_label_scf_state',
    'label_excited_states',
    'label_normal_modes',
    '_dump_mo_labels_log',
    '_dump_state_labels_log',
]


def load_branch_module(qualname, relpath, package=None):
    spec = importlib.util.spec_from_file_location(qualname, str(BRANCH_ROOT / relpath))
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    if package:
        module.__package__ = package
    sys.modules[qualname] = module
    spec.loader.exec_module(module)
    return module


def graft_branch_methods():
    """Install branch symmetry modules + Molecule methods into installed oqp."""
    import oqp.library  # ensure the package exists before overriding members

    detect = load_branch_module(
        'oqp.library.symmetry_detect', 'pyoqp/oqp/library/symmetry_detect.py')
    symmetry = load_branch_module(
        'oqp.library.symmetry', 'pyoqp/oqp/library/symmetry.py')
    oqp.library.symmetry_detect = detect
    oqp.library.symmetry = symmetry

    branch_molecule = load_branch_module(
        'oqp_branch_molecule_under_test', 'pyoqp/oqp/molecule/molecule.py',
        package='oqp.molecule')

    from oqp.molecule.molecule import Molecule as InstalledMolecule
    for name in GRAFTED_METHODS:
        setattr(InstalledMolecule, name, vars(branch_molecule.Molecule)[name])
    return InstalledMolecule


def run_case(workdir, name, content):
    from oqp.pyoqp import Runner

    case_dir = Path(workdir) / name
    case_dir.mkdir(parents=True, exist_ok=True)
    input_file = case_dir / f'{name}.inp'
    input_file.write_text(content)
    log = case_dir / f'{name}.log'

    runner = Runner(project=name, input_file=str(input_file), log=str(log),
                    usempi=False)
    runner.run()
    return runner.mol


def check(condition, message, failures):
    status = 'ok  ' if condition else 'FAIL'
    print(f'  [{status}] {message}')
    if not condition:
        failures.append(message)


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else '/tmp/oqp_symmetry_smoke'
    failures = []

    graft_branch_methods()

    print('== RHF water (6-31G*, diagonal input orientation) ==')
    mol = run_case(workdir, 'rhf_water', RHF_INPUT)
    mol.config['symmetry'] = {'enabled': 'true'}
    meta = mol.initialize_symmetry_metadata()
    check(meta.get('detected_point_group') == 'c2v',
          f"detected point group c2v (got {meta.get('detected_point_group')})", failures)
    check('detection_error' not in meta,
          f"no detection error (got {meta.get('detection_error')})", failures)

    mo = mol.label_molecular_orbitals()
    check(mo is not None and mo.get('status') == 'ok',
          f"MO labeling status ok (got {mo and mo.get('status')}, "
          f"error={mo and mo.get('error')})", failures)
    if mo and mo.get('status') == 'ok':
        labels = mo['alpha']['labels']
        nocc = 5
        print(f'  occupied alpha labels: {labels[:nocc]}')
        print(f"  max character deviation: {mo['alpha']['max_deviation']:.3e}")
        check('mixed' not in labels[:nocc], 'no mixed occupied orbitals', failures)
        check(labels[:nocc].count('a1') == 3,
              f"three a1 occupied orbitals (got {labels[:nocc].count('a1')})", failures)
        check(mo['alpha']['max_deviation'] < 1.0e-4,
              'occupied character deviation < 1e-4', failures)
        state = mo.get('scf_state')
        check(state is not None and state.get('term') == '1A1',
              f"SCF ground state 1A1 (got {state and state.get('term')})", failures)

    print('== MRSF water (BHHLYP/6-31G*, triplet ROHF reference) ==')
    mol = run_case(workdir, 'mrsf_water', MRSF_INPUT)
    mol.config['symmetry'] = {'enabled': 'true'}
    mol.initialize_symmetry_metadata()
    mo = mol.label_molecular_orbitals()
    check(mo is not None and mo.get('status') == 'ok',
          f"ROHF MO labeling status ok (got {mo and mo.get('status')})", failures)

    states = mol.label_excited_states()
    check(states is not None and states.get('status') == 'ok',
          f"state labeling status ok (got "
          f"{mol.symmetry_metadata.get('state_labels')})", failures)
    if states and states.get('status') == 'ok':
        print(f"  state labels: {states['labels']}")
        print(f"  transition labels: {states['transition_labels']}")
        print(f"  reference (SOMO) labels: {states['reference_labels']}")
        print(f"  max character deviation: {states['max_deviation']:.3e}")
        check(states['labels'][0] == 'a1',
              f"S0 is totally symmetric a1 (got {states['labels'][0]})", failures)
        check('mixed' not in states['labels'],
              'no mixed states among the 3 singlets', failures)

    print()
    if failures:
        print(f'SMOKE TEST FAILED ({len(failures)} failures):')
        for message in failures:
            print(f'  - {message}')
        return 1
    print('SMOKE TEST PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
