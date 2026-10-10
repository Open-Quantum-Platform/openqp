"""Run the same isolated backend checks on all three GitLab runner platforms."""
import os
from pathlib import Path
import platform
import subprocess
import sys
import venv

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = os.environ.get('STUDIO_CI_PLATFORM')
ACTUAL = {'Linux': 'linux', 'Darwin': 'macos', 'Windows': 'windows'}[platform.system()]
if EXPECTED != ACTUAL:
    raise SystemExit(f'Refusing mismatched runner: expected {EXPECTED}, found {ACTUAL}')
print(f'Backend CI: {ACTUAL}, {sys.version}', flush=True)
# Each GitLab job owns a separate checkout and venv. No user/site installation.
env_path = ROOT / '.cache' / 'ci-backend-venv'
if sys.version_info[:2] == (3, 12):
    venv.EnvBuilder(with_pip=True, clear=True).create(env_path)
else:
    # Shell runners may use Python 3.11 for OpenQP itself. Install Studio's
    # required interpreter inside this checkout without changing that runtime.
    bootstrap = ROOT / '.cache' / 'ci-bootstrap'
    venv.EnvBuilder(with_pip=True, clear=True).create(bootstrap)
    bootstrap_python = bootstrap / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
    subprocess.run([str(bootstrap_python), '-m', 'pip', 'install', 'uv>=0.8,<1'], check=True)
    uv = bootstrap / ('Scripts/uv.exe' if os.name == 'nt' else 'bin/uv')
    env = {**os.environ, 'UV_PYTHON_INSTALL_DIR': str(ROOT / '.cache' / 'ci-python'),
           'UV_CACHE_DIR': str(ROOT / '.cache' / 'uv'), 'UV_NATIVE_TLS': 'true'}
    subprocess.run([str(uv), 'venv', '--clear', '--seed', '--python', '3.12',
                    '--python-preference', 'only-managed', str(env_path)], env=env, check=True)
python = env_path / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
subprocess.run([str(python), '-c', 'import sys; print(sys.version); assert sys.version_info[:2] == (3, 12)'], check=True)
report = ROOT / 'reports'
report.mkdir(exist_ok=True)
commands = [
    ['-m', 'pip', 'install', '-e', '.[dev]'],
    ['-m', 'ruff', 'check', '.'],
    ['-m', 'pytest', '-v', f'--junitxml={report / "backend.xml"}'],
]
for args in commands:
    subprocess.run([str(python), *args], cwd=ROOT / 'backend', check=True)
