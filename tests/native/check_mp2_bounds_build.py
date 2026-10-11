"""Fail unless the fresh installed GNU package contains bounds-checked MP2.

Consumed by the x86 GCC pytest/MPI steps in .github/workflows/CI.yml. This
records compiler/build evidence; production execution is reported by the
separate MP2 suites, not inferred from the compile command.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shlex
import subprocess


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bounds_enabled(arguments):
    # A disabling clause removes our proof even on GNU versions that retain
    # legacy bounds instrumentation after -fcheck=no-bounds. Conservatively
    # rejecting such a command is safe; compiler quirks cannot admit a build.
    enabled = False
    for argument in arguments:
        if argument in ('-fbounds-check', '-fno-bounds-check'):
            enabled = argument == '-fbounds-check'
        if argument.startswith('@'):
            raise ValueError('unresolved compiler response file')
        if argument in ('-fcheck', '-fno-check'):
            raise ValueError('ambiguous GNU checking option')
        for prefix, positive in (('-fcheck=', True), ('-fno-check=', False)):
            if argument.startswith(prefix):
                for check in argument[len(prefix):].split(','):
                    negative = check.startswith('no-')
                    name = check[3:] if negative else check
                    if name not in ('all', 'bounds', 'array-temps', 'do', 'mem',
                                    'pointer', 'recursion'):
                        raise ValueError(f'unknown GNU checking option: {check}')
                    if name in ('all', 'bounds'):
                        enabled = positive != negative
    return enabled


def compile_entry(database, source):
    source = Path(source).resolve()
    matches = [entry for entry in database
               if (Path(entry['directory']) / entry['file']).resolve() == source]
    if len(matches) != 1:
        raise ValueError('expected exactly one production mp2_lib.F90 compilation')
    entry = matches[0]
    args = entry.get('arguments') or shlex.split(entry['command'])
    if not bounds_enabled(args):
        raise ValueError('MP2 compiler options do not establish the required bounds-check proof')
    if not re.fullmatch(r'gfortran(?:-\d+(?:\.\d+)*)?', Path(args[0]).name):
        raise ValueError('MP2 bounds proof requires the actual GNU Fortran compiler')
    outputs = [args[i + 1] for i, arg in enumerate(args[:-1]) if arg == '-o']
    if len(outputs) != 1:
        raise ValueError('missing or ambiguous MP2 object output')
    obj = (Path(entry['directory']) / outputs[0]).resolve()
    return entry, args, obj


def build_id(path):
    result = subprocess.run(['readelf', '-n', str(path)], check=True,
                            capture_output=True, text=True)
    ids = re.findall(r'Build ID:\s*([0-9a-f]+)', result.stdout)
    if len(ids) != 1:
        raise ValueError(f'expected one ELF build identity: {path}')
    return ids[0]


def verify(build, source, installed):
    build, source, installed = map(lambda p: Path(p).resolve(), (build, source, installed))
    start = int((build / 'package-build-start-ns').read_text())
    success = build / 'package-build-success'
    log = build / 'package-build.log'
    if not success.is_file() or not log.is_file() or not log.stat().st_size:
        raise ValueError('missing successful package-build marker or log')
    entry, args, obj = compile_entry(json.loads((build / 'compile_commands.json').read_text()),
                                     source)
    if not obj.is_relative_to(build):
        raise ValueError('MP2 object is outside the retained package build')
    built = build / 'source/liboqp.so'
    for path in (obj, built, installed):
        if not path.is_file() or not path.stat().st_size or path.stat().st_mtime_ns < start:
            raise ValueError(f'missing, empty or stale package artifact: {path}')
    if success.stat().st_mtime_ns < max(obj.stat().st_mtime_ns, built.stat().st_mtime_ns):
        raise ValueError('package success marker predates compilation/linking')
    identity = build_id(built)
    if build_id(installed) != identity:
        raise ValueError('installed library differs from this package build')
    compiler = subprocess.check_output([args[0], '--version'], text=True)
    version = subprocess.check_output([args[0], '-dumpfullversion'], text=True).strip()
    return {'evidence': 'compiler command and successful package build; execution separate',
            'compiler': compiler, 'compiler_version': version, 'compile_entry': entry,
            'source': str(source), 'source_sha256': sha256(source),
            'object': str(obj), 'object_sha256': sha256(obj),
            'built_library': str(built), 'built_sha256': sha256(built),
            'installed_library': str(installed), 'installed_sha256': sha256(installed),
            'elf_build_id': identity, 'build_log_sha256': sha256(log),
            'checkout_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                                    cwd=source.parents[1], text=True).strip()}


def assert_runtime(record_path, library):
    record = json.loads(Path(record_path).read_text())
    library = Path(library).resolve()
    if str(library) != record['installed_library'] or sha256(library) != record['installed_sha256']:
        raise ValueError('executed MP2 library differs from the bounds-checked package')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('record', type=Path)
    args = parser.parse_args()
    # Do not import oqp or initialize MPI in the launcher.
    dist = importlib.metadata.distribution('OpenQP')
    libs = [dist.locate_file(f) for f in dist.files if str(f) == 'oqp/lib/liboqp.so']
    if len(libs) != 1:
        raise ValueError('expected one installed package liboqp.so')
    source = Path(__file__).resolve().parents[2] / 'source/mp2_lib.F90'
    record = verify(args.build, source, libs[0])
    args.record.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
