"""Phase IV validation: irrep-blocked Davidson vs unblocked excitation energies.

Not collected by pytest. Run with the branch's built venv:

    .venv-symmetry/bin/python3 tests/smoke_response_symmetry_validation.py [workdir]

Gate: excitation energies agree with the unblocked solver to <= 1e-8 Ha
and the state labels are unchanged.
"""

import sys
from pathlib import Path

import numpy as np

GATE = 1.0e-8

WATER = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

# (scf block, extra input lines, tdhf block)
CASES = {
    'water_mrsf_3': ('multiplicity=3\ntype=rohf', 'functional=bhhlyp\n',
                     '[tdhf]\ntype=mrsf\nnstate=3\n\n'),
    'water_mrsf_5': ('multiplicity=3\ntype=rohf', 'functional=bhhlyp\n',
                     '[tdhf]\ntype=mrsf\nnstate=5\n\n'),
    'water_tda_5': ('multiplicity=1\ntype=rhf', 'functional=bhhlyp\n',
                    '[tdhf]\ntype=tda\nnstate=5\n\n'),
    'water_rpa_4': ('multiplicity=1\ntype=rhf', 'functional=bhhlyp\n',
                    '[tdhf]\ntype=rpa\nnstate=4\n\n'),
}

INPUT_TEMPLATE = """\
[input]
system=
{system}
charge=0
runtype=energy
basis=6-31g*
{extra}method=tdhf

[guess]
type=huckel

[scf]
{scf}

{tdhf}[symmetry]
{symmetry}
"""


def run_case(workdir, name, content):
    from oqp.pyoqp import Runner

    case_dir = Path(workdir) / name
    case_dir.mkdir(parents=True, exist_ok=True)
    input_file = case_dir / f'{name}.inp'
    input_file.write_text(content)
    runner = Runner(project=name, input_file=str(input_file),
                    log=str(case_dir / f'{name}.log'), usempi=False)
    runner.run()
    return runner.mol


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else '/tmp/oqp_response_sym'
    failures = []

    for name, (scf, extra, tdhf) in CASES.items():
        ref_input = INPUT_TEMPLATE.format(
            system=WATER, scf=scf, extra=extra, tdhf=tdhf,
            symmetry='enabled=false')
        sym_input = INPUT_TEMPLATE.format(
            system=WATER, scf=scf, extra=extra, tdhf=tdhf,
            symmetry='enabled=true\nuse_response_symmetry=true')

        mol_ref = run_case(workdir, f'{name}_ref', ref_input)
        e_ref = np.asarray(mol_ref.energies, dtype=float)

        mol_sym = run_case(workdir, f'{name}_blocked', sym_input)
        e_sym = np.asarray(mol_sym.energies, dtype=float)

        meta = mol_sym.symmetry_metadata
        active = meta.get('response_symmetry', {})
        labels = meta.get('state_labels', {}).get('labels')

        diff = float(np.max(np.abs(e_sym - e_ref))) if e_ref.shape == e_sym.shape \
            else float('inf')

        ok = diff <= GATE and active.get('status') == 'active'
        status = 'ok  ' if ok else 'FAIL'
        print(f'[{status}] {name:14s} response={active.get("status")} '
              f'max|dE|={diff:.3e} labels={labels}')
        if not ok:
            failures.append(name)

    print()
    if failures:
        print(f'RESPONSE VALIDATION FAILED: {failures}')
        return 1
    print('RESPONSE VALIDATION PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
