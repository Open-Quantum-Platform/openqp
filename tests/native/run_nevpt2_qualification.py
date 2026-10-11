#!/usr/bin/env python3
"""Run named installed-library NEVPT2 regressions and retain their CI receipts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


TESTS = tuple(
    f'{test}[{name}]'
    for test in ('test_final_contraction_matches_numpy',
                 'test_retained_subspace_outputs_survive_later_calls')
    for name in ('Sr', 'Si', 'Sijrs', 'Sijr', 'Srsi', 'Srs', 'Sij', 'Sir')
) + tuple(f'test_fused_rdms_and_f3_match_numpy[{name}]'
          for name in ('cas2', 'cas3', 'cas4')) + (
    'test_retained_fused_outputs_survive_later_calls',
    'test_actual_native_allocation_failures_and_recovery',
    'test_malloc_fault_shim_survives_late_tls_modules',
    'test_whole_energy_matches_backend_unavailable[core1-cas2-virt2]',
    'test_whole_energy_matches_backend_unavailable[core0-cas3-virt2]',
)


def require_named_passes(junit):
    """A success exit/aggregate alone cannot qualify any selected native node."""
    cases = list(ET.parse(junit).getroot().iter('testcase'))
    if not cases:
        raise RuntimeError('No named NEVPT2 test cases recorded')
    for case in cases:
        if any(case.find(tag) is not None for tag in ('skipped', 'failure', 'error')):
            raise RuntimeError(f"NEVPT2 node did not pass: {case.attrib}")
    identities = [(case.get('classname'), case.get('name')) for case in cases]
    if len(set(identities)) != len(identities):
        raise RuntimeError('Duplicate named NEVPT2 results cannot replace missing parameters')
    expected = {('tests.test_nevpt2_subspaces', name) for name in TESTS}
    actual = set(identities)
    if actual != expected:
        raise RuntimeError(f'Missing or unexpected named NEVPT2 results: '
                           f'missing={sorted(expected - actual)}, extra={sorted(actual - expected)}')
    return [case.attrib for case in cases]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, choices=(1, 4), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    # Thread policy is supplied by the invoking CI step before Python/oqp loads.
    # Do not silently overwrite a mismatch or relabel it as a threaded run.
    for key, expected in (('OMP_NUM_THREADS', str(args.threads)),
                          ('OQP_FCI_THREADS', str(args.threads)),
                          ('OMP_DYNAMIC', 'FALSE'), ('OPENBLAS_NUM_THREADS', '1'),
                          ('MKL_NUM_THREADS', '1')):
        if os.environ.get(key) != expected:
            raise RuntimeError(f'{key} must be {expected}, got {os.environ.get(key)!r}')
    import oqp
    from oqp.library import fci, nevpt2_sc
    from oqp.runtime import library_path

    lib, ffi = nevpt2_sc._koopmans_lib()
    if lib is not oqp.lib or ffi is not oqp.ffi:
        raise RuntimeError('NEVPT2 backend differs from the installed oqp library')
    library = library_path(oqp.oqp_root, oqp.suffix).resolve()
    symbols = ['nevpt2_make_f3', 'nevpt2_f3ca_f3ac_checked'] + [
        'nevpt2_' + name for name in ('sr', 'si', 'sijrs', 'sijr', 'srsi', 'srs', 'sij', 'sir')]
    for symbol in symbols:
        getattr(lib, symbol)  # Missing installed exports fail rather than skip.
    if not oqp.OMP_THREADS_FROM_ENV or fci._fci_lib_threads() != args.threads:
        raise RuntimeError('Installed FCI thread request does not match the CI invocation')
    args.output.mkdir(parents=True, exist_ok=False)
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    provenance = {
        'checkout_sha': source, 'python': sys.version, 'executable': sys.executable,
        'library': str(library), 'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(),
        'required_symbols': symbols, 'omp_requested': args.threads,
        'fci_thread_argument': fci._fci_lib_threads(),
        'omp_threads_from_env': oqp.OMP_THREADS_FROM_ENV,
        'observed_openmp_team_size': None,
        'qualification_limit': 'Thread requests are checked; native team sizes and performance are not measured.',
    }
    (args.output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance), flush=True)
    junit = args.output / 'junit.xml'
    nodes = ['tests/test_nevpt2_subspaces.py::' + test for test in TESTS]
    env = os.environ.copy()
    env['OQP_NEVPT2_RECEIPT_DIR'] = str(args.output.resolve())
    env['OQP_NEVPT2_EXPECTED_LIBRARY'] = str(library)
    env['OQP_NEVPT2_EXPECTED_LIBRARY_SHA256'] = provenance['library_sha256']
    env['OQP_NEVPT2_EXPECTED_PACKAGE'] = str(Path(oqp.__file__).resolve().parent)
    result = subprocess.run([sys.executable, '-m', 'pytest', *nodes, '-vv', '-rA',
                             '--tb=short', f'--junitxml={junit}'], env=env)
    if result.returncode:
        return result.returncode
    cases = require_named_passes(junit)
    for filename in ('pytest-runtime.json', 'fault-runtime.json'):
        runtime = json.loads((args.output / filename).read_text())
        if (runtime['library'] != provenance['library']
                or runtime['library_sha256'] != provenance['library_sha256']):
            raise RuntimeError(f'{filename}: actual subprocess library identity differs')
    (args.output / 'named-passes.json').write_text(json.dumps(cases, indent=2) + '\n')
    print(f'NEVPT2 named qualification passed: {len(cases)} cases, thread request {args.threads}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
