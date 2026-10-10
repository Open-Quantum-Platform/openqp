import asyncio
import io
import json
import subprocess
from pathlib import Path

from oqp_studio import engine
from oqp_studio.runners.wsl import WslRunner
from oqp_studio.server_main import _serve_requests


def test_slow_request_does_not_block_health_reply():
    async def scenario():
        release = asyncio.Event()
        replies = []

        async def dispatch(request):
            if request['path'] == '/slow':
                await asyncio.wait_for(release.wait(), timeout=2)
            else:
                release.set()
            return {'path': request['path']}

        data = ''.join(json.dumps({'id': i, 'request': {'path': p}}) + '\n'
                       for i, p in [(1, '/slow'), (2, '/api/health')])
        await _serve_requests(io.StringIO(data).readline,
                              lambda line: replies.append(json.loads(line)), dispatch)
        assert [reply['id'] for reply in replies] == [2, 1]
    asyncio.run(scenario())


def test_requests_are_bounded_and_malformed_id_is_not_reused():
    async def scenario():
        active = peak = 0
        replies = []

        async def dispatch(_request):
            nonlocal active, peak
            active += 1
            peak = max(peak, active)
            await asyncio.sleep(0.01)
            active -= 1
            return {}

        data = ''.join(json.dumps({'id': i, 'request': {}}) + '\n' for i in range(12))
        data += 'invalid-json\n[]\n'
        await _serve_requests(io.StringIO(data).readline,
                              lambda line: replies.append(json.loads(line)), dispatch, limit=2)
        assert peak == 2
        assert {r['id'] for r in replies if 'result' in r} == set(range(12))
        assert [r['id'] for r in replies if 'error' in r] == [None, None]
    asyncio.run(scenario())


def test_version_does_not_read_an_entire_binary(tmp_path, monkeypatch):
    executable = tmp_path / 'openqp'
    executable.write_bytes(b'\x7fELF' + b'\x00' * 8192)
    original_open = Path.open

    class LimitedReader(io.BytesIO):
        def read(self, size=-1):
            assert 0 < size <= 4096
            return super().read(size)

    def guarded_open(path, *args, **kwargs):
        if path == executable:
            return LimitedReader(b'\x7fELF' + b'\x00' * 8192)
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', guarded_open)
    assert engine.version(str(executable)) is None


def test_wsl_timeout_is_unavailable_not_discovery_failure(monkeypatch):
    from oqp_studio.runners import wsl
    monkeypatch.setattr(wsl.sys, 'platform', 'win32')
    monkeypatch.setattr(wsl.shutil, 'which', lambda _name: 'wsl.exe')

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired('wsl', 15)

    monkeypatch.setattr(wsl.subprocess, 'run', timeout)
    assert WslRunner().is_available() is False
