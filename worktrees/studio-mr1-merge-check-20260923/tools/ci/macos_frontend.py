"""Install verified Node 24 into this job's checkout and run frontend CI."""
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
if platform.system() != 'Darwin':
    raise SystemExit('This bootstrap is for a native macOS runner')
ARCH = {'arm64': 'arm64', 'x86_64': 'x64'}[platform.machine()]
VERSION = 'v24.19.0'
NAME = f'node-{VERSION}-darwin-{ARCH}'
ARCHIVE = NAME + '.tar.gz'
BASE = f'https://nodejs.org/dist/{VERSION}/'
CACHE = ROOT / '.cache' / 'ci-node'
CACHE.mkdir(parents=True, exist_ok=True)
with urllib.request.urlopen(BASE + 'SHASUMS256.txt', timeout=60) as response:
    checksums = dict(line.split()[::-1] for line in response.read().decode().splitlines())
archive = CACHE / ARCHIVE
urllib.request.urlretrieve(BASE + ARCHIVE, archive)
if hashlib.sha256(archive.read_bytes()).hexdigest() != checksums[ARCHIVE]:
    raise SystemExit('Node archive checksum mismatch')
with tarfile.open(archive) as stream:
    stream.extractall(CACHE, filter='data')
bin_dir = CACHE / NAME / 'bin'
env = {**os.environ, 'PATH': str(bin_dir) + os.pathsep + os.environ.get('PATH', '')}
subprocess.run([str(bin_dir / 'node'), str(ROOT / 'tools' / 'ci' / 'frontend.mjs')],
               cwd=ROOT, env=env, check=True)
