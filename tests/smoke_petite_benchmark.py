"""Petite-list performance benchmark: C1 vs symmetry-reduced wall times.

Not collected by pytest. Run with an OpenMP-enabled built venv:

    OMP_NUM_THREADS=8 .venv-symmetry/bin/python3 tests/smoke_petite_benchmark.py [workdir]

Reports wall times for the full energy run and checks the energies still
agree (gate 1e-9, looser than the validation suite since geometries here
are idealized and SCF iteration counts may differ by one).
"""

import sys
import time
from pathlib import Path

import numpy as np

GATE = 1.0e-9


def naphthalene():
    """Idealized D2h naphthalene (rings of side 1.40 A, C-H 1.09 A)."""
    a = 1.40
    h = 1.09
    half = a * np.sqrt(3.0) / 2.0  # 1.2124
    carbons = [(0.0, a / 2), (0.0, -a / 2)]
    for sx in (1.0, -1.0):
        carbons += [
            (sx * half, a),
            (sx * half, -a),
            (sx * 2 * half, a / 2),
            (sx * 2 * half, -a / 2),
        ]
    hydrogens = []
    for sx in (1.0, -1.0):
        for sy in (1.0, -1.0):
            hydrogens.append((sx * half, sy * (a + h)))
            ux, uy = np.array([half, a / 2]) / a  # radial unit from ring center
            hydrogens.append((sx * (2 * half + h * ux), sy * (a / 2 + h * uy)))
    lines = [f'   6 {x:13.6f} {y:13.6f}   0.000000' for x, y in carbons]
    lines += [f'   1 {x:13.6f} {y:13.6f}   0.000000' for x, y in hydrogens]
    return '\n'.join(lines)


CASES = {
    # name: (extra input lines, basis)
    'naphthalene_rhf': ('', '6-31g*'),
    'naphthalene_dft': ('functional=bhhlyp\n', '6-31g*'),
}

INPUT_TEMPLATE = """\
[input]
system=
{system}
charge=0
runtype=energy
basis={basis}
{extra}method=hf

[guess]
type=huckel

[scf]
multiplicity=1
type=rhf

[symmetry]
{symmetry}
"""


def run_case(workdir, name, content):
    from oqp.pyoqp import Runner

    case_dir = Path(workdir) / name
    case_dir.mkdir(parents=True, exist_ok=True)
    input_file = case_dir / f'{name}.inp'
    input_file.write_text(content)
    start = time.time()
    runner = Runner(project=name, input_file=str(input_file),
                    log=str(case_dir / f'{name}.log'), usempi=False)
    runner.run()
    elapsed = time.time() - start
    return runner.mol, elapsed


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else '/tmp/oqp_petite_bench'
    system = naphthalene()
    failures = []

    for name, (extra, basis) in CASES.items():
        ref_input = INPUT_TEMPLATE.format(system=system, extra=extra,
                                          basis=basis, symmetry='enabled=false')
        sym_input = INPUT_TEMPLATE.format(
            system=system, extra=extra, basis=basis,
            symmetry='enabled=true\nuse_integral_symmetry=true')

        mol_ref, t_ref = run_case(workdir, f'{name}_c1', ref_input)
        e_ref = float(mol_ref.energies[0])
        mol_sym, t_sym = run_case(workdir, f'{name}_petite', sym_input)
        e_sym = float(mol_sym.energies[0])

        active = mol_sym.symmetry_metadata.get('integral_symmetry', {})
        diff = abs(e_sym - e_ref)
        speedup = t_ref / t_sym if t_sym > 0 else float('nan')
        ok = diff <= GATE and active.get('status') == 'active'
        status = 'ok  ' if ok else 'FAIL'
        print(f'[{status}] {name:18s} nops={active.get("n_operations")} '
              f'E_C1={e_ref:.9f} dE={diff:.2e} '
              f't_C1={t_ref:7.1f}s t_petite={t_sym:7.1f}s speedup={speedup:4.2f}x')
        if not ok:
            failures.append(name)

    print()
    if failures:
        print(f'BENCHMARK CHECK FAILED: {failures}')
        return 1
    print('BENCHMARK DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())
