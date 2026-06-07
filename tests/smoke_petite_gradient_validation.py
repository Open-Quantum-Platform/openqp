"""Phase II validation: petite-list gradients vs C1 references.

Not collected by pytest. Run with the branch's built venv:

    .venv-symmetry/bin/python3 tests/smoke_petite_gradient_validation.py [workdir]

All geometries are pre-rotated to the standard orientation so the C1
reference and the petite run share the frame (XC grids and the gradient
direction are orientation-dependent). Gate: max |dG| <= 1e-9 Hartree/bohr.
"""

import sys
from pathlib import Path

import numpy as np

GATE = 1.0e-9

WATER = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

ETHYLENE = """\
   6   0.000000000   0.000000000   0.665000000
   6   0.000000000   0.000000000  -0.665000000
   1   0.000000000   0.920000000   1.230000000
   1   0.000000000  -0.920000000   1.230000000
   1   0.000000000   0.920000000  -1.230000000
   1   0.000000000  -0.920000000  -1.230000000"""

# (scf block, extra input lines, tdhf block or '')
CASES = {
    'water_rhf_grad': (WATER, 'multiplicity=1\ntype=rhf', '', ''),
    'ethylene_rhf_grad': (ETHYLENE, 'multiplicity=1\ntype=rhf', '', ''),
    'water_dft_grad': (WATER, 'multiplicity=1\ntype=rhf', 'functional=bhhlyp\n', ''),
    'water_mrsf_s1_grad': (WATER, 'multiplicity=3\ntype=rohf',
                           'functional=bhhlyp\n',
                           '[tdhf]\ntype=mrsf\nnstate=3\n\n[properties]\ngrad=1\n'),
}

INPUT_TEMPLATE = """\
[input]
system=
{system}
charge=0
runtype=grad
basis=6-31g*
{extra}method={method}

[guess]
type=huckel

[scf]
{scf}

{tdhf}[symmetry]
{symmetry}
"""


def standardize_orientation(system):
    from oqp.library.symmetry_detect import detect_point_group

    charges, coords = [], []
    for line in system.splitlines():
        parts = line.split()
        charges.append(float(parts[0]))
        coords.append([float(x) for x in parts[1:4]])
    charges = np.array(charges)
    coords = np.array(coords)
    detection = detect_point_group(charges, coords, tolerance=1.0e-6)
    rotation = np.array(detection['orientation'])
    origin = np.array(detection['origin'])
    standard = (coords - origin) @ rotation.T
    return '\n'.join(
        f'   {int(q)} {x:15.9f} {y:15.9f} {z:15.9f}'
        for q, (x, y, z) in zip(charges, standard)
    )


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
    workdir = sys.argv[1] if len(sys.argv) > 1 else '/tmp/oqp_petite_grad'
    failures = []

    for name, (system, scf, extra, tdhf) in CASES.items():
        system = standardize_orientation(system)
        method = 'tdhf' if tdhf else 'hf'
        ref_input = INPUT_TEMPLATE.format(
            system=system, scf=scf, extra=extra, tdhf=tdhf, method=method,
            symmetry='enabled=false')
        sym_input = INPUT_TEMPLATE.format(
            system=system, scf=scf, extra=extra, tdhf=tdhf, method=method,
            symmetry='enabled=true\nuse_integral_symmetry=true')

        mol_ref = run_case(workdir, f'{name}_c1', ref_input)
        g_ref = np.asarray(mol_ref.grads, dtype=float)

        mol_sym = run_case(workdir, f'{name}_petite', sym_input)
        g_sym = np.asarray(mol_sym.grads, dtype=float)

        meta = mol_sym.symmetry_metadata
        active = meta.get('integral_symmetry', {})
        diff = float(np.max(np.abs(g_sym - g_ref))) if g_ref.shape == g_sym.shape \
            else float('inf')
        gnorm = float(np.max(np.abs(g_ref)))

        ok = diff <= GATE and active.get('status') == 'active'
        status = 'ok  ' if ok else 'FAIL'
        print(f'[{status}] {name:22s} petite={active.get("status")} '
              f'max|G_C1|={gnorm:.3e} max|dG|={diff:.3e}')
        if not ok:
            failures.append(name)

    print()
    if failures:
        print(f'GRADIENT VALIDATION FAILED: {failures}')
        return 1
    print('GRADIENT VALIDATION PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
