"""Negative regressions for the CI's package instrumentation evidence guard."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest


spec = importlib.util.spec_from_file_location(
    'mp2_bounds_guard', Path(__file__).with_name('native') / 'check_mp2_bounds_build.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


@pytest.mark.parametrize('flags,enabled', [
    ([], False), (['-fcheck=bounds'], True), (['-fcheck=all'], True),
    (['-fcheck=bounds', '-fcheck=pointer'], True),
    (['-fcheck=bounds', '-fcheck=no-bounds'], False),
    (['-fcheck=all,no-bounds'], False), (['-fcheck=no-bounds,all'], True),
    (['-fcheck=bounds', '-fno-check=bounds'], False),
    (['-fcheck=bounds', '-fno-check=all'], False),
    (['-fcheck=all', '-fcheck=no-all'], False),
    (['-fno-check=all', '-fcheck=bounds'], True),
    (['-fcheck=bounds', '-fno-bounds-check'], False),
    (['-fno-bounds-check', '-fcheck=bounds'], True),
    (['-fbounds-check'], True),
    (['-fbounds-check', '-fno-bounds-check'], False),
])
def test_effective_ordered_flags(flags, enabled):
    assert guard.bounds_enabled(flags) is enabled


@pytest.mark.parametrize('flags', [['@flags.rsp'], ['-fcheck'], ['-fcheck=unknown']])
def test_ambiguous_instrumentation_fails(flags):
    with pytest.raises(ValueError):
        guard.bounds_enabled(flags)


@pytest.mark.parametrize('flags', [
    ['-fcheck=bounds'], ['-fcheck=bounds', '-fcheck=pointer'],
    ['-fcheck=bounds', '-fcheck=no-bounds'], ['-fcheck=all,no-bounds'],
    ['-fcheck=no-bounds,all'], ['-fcheck=all', '-fcheck=no-all'],
    ['-fcheck=no-all', '-fcheck=bounds'],
    ['-fcheck=bounds', '-fno-bounds-check'],
    ['-fno-bounds-check', '-fcheck=bounds'], ['-fbounds-check'],
    ['-fbounds-check', '-fno-bounds-check'],
])
def test_guard_never_accepts_uninstrumented_gnu_options(tmp_path, flags):
    compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
    assert compiler, 'GNU Fortran is required to check actual ordered flag semantics'
    source = tmp_path / 'bounds.f90'
    source.write_text('''subroutine access(a,n,i,out)
 implicit none
 integer, intent(in) :: n,i
 real(8), intent(in) :: a(n)
 real(8), intent(out) :: out
 out=a(i)
end subroutine
''')
    assembly = tmp_path / 'bounds.s'
    result = subprocess.run([compiler, '-S', '-O0', *flags, str(source), '-o', str(assembly)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    # -fcheck=all also emits recursion diagnostics through runtime_error_at;
    # require the actual lower/upper bounds messages, not only that symbol.
    # GNU 15.3 can retain checks after no-bounds; the guard conservatively
    # rejects those commands. Rejection is not a claim of absent diagnostics.
    text = assembly.read_text()
    instrumented = ('_gfortran_runtime_error_at' in text and
                    'below lower bound' in text and 'above upper bound' in text)
    if guard.bounds_enabled(flags):
        assert instrumented


@pytest.fixture
def entry(tmp_path):
    source = tmp_path / 'source/mp2_lib.F90'
    return source, {'directory': str(tmp_path), 'file': str(source),
                    'command': f'gfortran-14 -fcheck=bounds -c {source} -o mp2.o'}


def test_exact_source_and_object(entry):
    source, row = entry
    _, _, obj = guard.compile_entry([row], source)
    assert obj == source.parents[1] / 'mp2.o'


@pytest.mark.parametrize('change', ['wrong-source', 'duplicate', 'disabled', 'no-object', 'non-gnu'])
def test_invalid_compile_proof_fails(entry, change):
    source, row = entry
    rows = [row]
    if change == 'wrong-source':
        row['file'] = 'elsewhere/mp2_lib.F90'
    elif change == 'duplicate':
        rows.append(row.copy())
    elif change == 'disabled':
        row['command'] += ' -fcheck=no-bounds'
    elif change == 'no-object':
        row['command'] = row['command'].split(' -o ')[0]
    else:
        row['command'] = row['command'].replace('gfortran-14', 'ifx')
    with pytest.raises(ValueError):
        guard.compile_entry(rows, source)


def test_runtime_cannot_use_another_library(tmp_path):
    library = tmp_path / 'liboqp.so'
    library.write_bytes(b'candidate')
    record = tmp_path / 'proof.json'
    record.write_text(json.dumps({'installed_library': str(library),
                                 'installed_sha256': guard.sha256(library)}))
    guard.assert_runtime(record, library)
    library.write_bytes(b'another-build')
    with pytest.raises(ValueError, match='differs'):
        guard.assert_runtime(record, library)


@pytest.fixture
def package_proof(tmp_path, monkeypatch):
    source = tmp_path / 'repo/source/mp2_lib.F90'
    source.parent.mkdir(parents=True)
    source.write_text('production')
    build = tmp_path / 'build'
    (build / 'source').mkdir(parents=True)
    (build / 'package-build-start-ns').write_text('1')
    (build / 'package-build.log').write_text('successful package build')
    obj = build / 'mp2.o'
    obj.write_bytes(b'object')
    built = build / 'source/liboqp.so'
    built.write_bytes(b'linked')
    installed = tmp_path / 'installed.so'
    installed.write_bytes(b'installed')
    (build / 'compile_commands.json').write_text(json.dumps([
        {'directory': str(build), 'file': str(source),
         'arguments': ['gfortran-14', '-fcheck=bounds', '-c', str(source), '-o', str(obj)]}]))
    (build / 'package-build-success').touch()
    monkeypatch.setattr(guard, 'build_id', lambda path: 'matching-id')
    monkeypatch.setattr(guard.subprocess, 'check_output', lambda *a, **k: 'compiler-or-checkout')
    return build, source, installed, obj


def test_package_proof_records_exact_artifacts(package_proof):
    build, source, installed, obj = package_proof
    record = guard.verify(build, source, installed)
    assert record['source_sha256'] == guard.sha256(source)
    assert record['object_sha256'] == guard.sha256(obj)
    assert record['installed_sha256'] == guard.sha256(installed)
    assert record['elf_build_id'] == 'matching-id'


@pytest.mark.parametrize('damage', ['no-success', 'no-log', 'no-object', 'stale-object',
                                   'early-success', 'wrong-build-id'])
def test_incomplete_package_proof_fails(package_proof, monkeypatch, damage):
    build, source, installed, obj = package_proof
    if damage == 'no-success':
        (build / 'package-build-success').unlink()
    elif damage == 'no-log':
        (build / 'package-build.log').unlink()
    elif damage == 'no-object':
        obj.unlink()
    elif damage == 'stale-object':
        (build / 'package-build-start-ns').write_text(str(obj.stat().st_mtime_ns + 1))
    elif damage == 'early-success':
        os.utime(build / 'package-build-success', ns=(0, 0))
    else:
        monkeypatch.setattr(guard, 'build_id', lambda path: 'other' if path == installed else 'id')
    with pytest.raises(ValueError):
        guard.verify(build, source, installed)
