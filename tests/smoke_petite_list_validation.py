"""Phase I validation: petite-list/skeleton-Fock vs C1 reference energies.

Not collected by pytest. Run with the branch's built venv:

    .venv-symmetry/bin/python3 tests/smoke_petite_list_validation.py [workdir]

Each case runs twice with the same input geometry: once with symmetry
disabled (C1 reference) and once with use_integral_symmetry enabled
(reorientation + petite list + skeleton-Fock symmetrization). Gate:
energy agreement <= 1e-10 Hartree.
"""

import sys
from pathlib import Path

GATE = 1.0e-10

CASES = {
    # name: (expected subgroup, scf block, geometry block)
    'water_rhf_c2v': ('c2v', 'multiplicity=1\ntype=rhf', """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""),
    'ethylene_rhf_d2h': ('d2h', 'multiplicity=1\ntype=rhf', """\
   6   0.000000000   0.000000000   0.665000000
   6   0.000000000   0.000000000  -0.665000000
   1   0.000000000   0.920000000   1.230000000
   1   0.000000000  -0.920000000   1.230000000
   1   0.000000000   0.920000000  -1.230000000
   1   0.000000000  -0.920000000  -1.230000000"""),
    'tfde_rhf_c2h': ('c2h', 'multiplicity=1\ntype=rhf', """\
   6   0.660000000   0.000000000   0.000000000
   6  -0.660000000   0.000000000   0.000000000
   9   1.400000000   1.100000000   0.000000000
   9  -1.400000000  -1.100000000   0.000000000
   1   1.200000000  -1.000000000   0.000000000
   1  -1.200000000   1.000000000   0.000000000"""),
    'h2o2_rhf_c2': ('c2', 'multiplicity=1\ntype=rhf', """\
   8   0.700000000   0.000000000   0.000000000
   8  -0.700000000   0.000000000   0.000000000
   1   1.000000000   0.700000000   0.500000000
   1  -1.000000000  -0.700000000   0.500000000"""),
    'water_rohf_triplet': ('c2v', 'multiplicity=3\ntype=rohf', """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""),
    'water_uhf_triplet': ('c2v', 'multiplicity=3\ntype=uhf', """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""),
}


def benzene_system():
    import numpy as np

    lines = []
    for k in range(6):
        angle = np.deg2rad(60.0 * k)
        lines.append(f'   6 {1.39*np.cos(angle):15.9f} {1.39*np.sin(angle):15.9f}   0.000000000')
        lines.append(f'   1 {2.49*np.cos(angle):15.9f} {2.49*np.sin(angle):15.9f}   0.000000000')
    return '\n'.join(lines)


# Benzene: D6h molecule, d2h petite subgroup (|G|=8) -- the canonical win.
CASES['benzene_rhf_d2h'] = ('d2h', 'multiplicity=1\ntype=rhf', None)

INPUT_TEMPLATE = """\
[input]
system=
{system}
charge=0
runtype=energy
basis=6-31g*
method=hf

[guess]
type=huckel

[scf]
{scf}

[symmetry]
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
    workdir = sys.argv[1] if len(sys.argv) > 1 else '/tmp/oqp_petite_validation'
    failures = []

    import time

    for name, (subgroup, scf, system) in CASES.items():
        if system is None:
            system = benzene_system()
        ref_input = INPUT_TEMPLATE.format(system=system, scf=scf,
                                          symmetry='enabled=false')
        sym_input = INPUT_TEMPLATE.format(
            system=system, scf=scf,
            symmetry='enabled=true\nuse_integral_symmetry=true')

        t0 = time.time()
        mol_ref = run_case(workdir, f'{name}_c1', ref_input)
        t_ref = time.time() - t0
        e_ref = float(mol_ref.energies[0])

        t0 = time.time()
        mol_sym = run_case(workdir, f'{name}_petite', sym_input)
        t_sym = time.time() - t0
        e_sym = float(mol_sym.energies[0])

        meta = mol_sym.symmetry_metadata
        active = meta.get('integral_symmetry', {})
        detected = meta.get('detected_subgroup')
        diff = abs(e_sym - e_ref)

        # 'disabled_symmetry_broken_scf' is the documented fail-safe when
        # the stability check keeps a broken-symmetry solution (energy must
        # then match C1 trivially and exactly).
        ok = (diff <= GATE
              and active.get('status') in ('active', 'disabled_symmetry_broken_scf')
              and detected == subgroup)
        status = 'ok  ' if ok else 'FAIL'
        print(f'[{status}] {name:24s} subgroup={detected} '
              f'(expected {subgroup}) petite={active.get("status")} '
              f'nops={active.get("n_operations")} '
              f'E_C1={e_ref:.12f} dE={diff:.2e} '
              f't_C1={t_ref:.1f}s t_sym={t_sym:.1f}s')
        if not ok:
            failures.append(name)

    print()
    if failures:
        print(f'PETITE VALIDATION FAILED: {failures}')
        return 1
    print('PETITE VALIDATION PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
