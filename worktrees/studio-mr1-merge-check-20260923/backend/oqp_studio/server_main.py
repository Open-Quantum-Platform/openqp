"""Standalone server entry point — used by the PyInstaller-frozen backend.

The Tauri shell spawns this binary as a sidecar:
    oqp-studio-backend --port 8814
"""

from __future__ import annotations

import argparse
import asyncio
import base64
import json
import os
import sys
from pathlib import Path


def log_path() -> Path:
    """Where a startup failure is recorded.

    Worked out here rather than by importing the app's own settings module,
    because the failure being recorded may be that module failing to import.
    """
    override = os.environ.get("OQP_STUDIO_CONFIG")
    if override:
        return Path(override).parent / "backend.log"
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support" / "OQP Studio"
    elif os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home())) / "OQP Studio"
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME",
                                   Path.home() / ".config")) / "oqp-studio"
    return base / "backend.log"


async def _dispatch(request: dict) -> dict:
    """Run one API request in-process without opening a TCP listener."""
    import httpx

    from oqp_studio.main import app

    body = request.get("body")
    content = base64.b64decode(body) if body else None
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://oqp-studio.internal") as client:
        response = await client.request(
            request["method"],
            request["path"],
            headers=request.get("headers") or {},
            content=content,
        )
    return {
        "status": response.status_code,
        "headers": dict(response.headers),
        "body": base64.b64encode(response.content).decode("ascii"),
    }


async def _serve_requests(readline, write, dispatch=_dispatch, limit: int = 8) -> None:
    """Bound concurrent requests so a slow probe cannot block health/cancel.

    All replies are written on this event loop, never interleaved by threads.
    Request IDs allow replies to arrive in a different order from requests.
    """
    pending: set[asyncio.Task] = set()

    async def respond(line: str) -> None:
        message = {}
        try:
            message = json.loads(line)
            result = await dispatch(message["request"])
            payload = {"id": message["id"], "result": result}
        except Exception as exc:  # noqa: BLE001 - isolate per-request failures
            payload = {"id": message.get("id") if isinstance(message, dict) else None,
                       "error": str(exc)}
        write(json.dumps(payload, separators=(",", ":")) + "\n")

    while True:
        if len(pending) >= limit:
            done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
            for task in done:
                task.result()
        line = await asyncio.to_thread(readline)
        if not line:
            break
        # Reap completed tasks, preserving exceptions from the output stream.
        done = {task for task in pending if task.done()}
        for task in done:
            task.result()
        pending -= done
        pending.add(asyncio.create_task(respond(line)))
    if pending:
        await asyncio.gather(*pending)


async def _serve_stdio() -> None:
    from oqp_studio.main import app

    def write(payload: str) -> None:
        sys.stdout.write(payload)
        sys.stdout.flush()

    # httpx.ASGITransport does not start lifespan itself. Use one lifecycle
    # and event loop for the process, including deferred runner warm-up.
    async with app.router.lifespan_context(app):
        await _serve_requests(sys.stdin.readline, write)


def serve_stdio() -> None:
    """Serve newline-delimited JSON RPC over the sidecar's stdio streams."""
    asyncio.run(_serve_stdio())


def main() -> None:
    parser = argparse.ArgumentParser(description="OQP Studio backend server")
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("OQP_STUDIO_PORT", "8814")),
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--stdio", action="store_true",
                        help="serve API requests through stdin/stdout instead of TCP")
    parser.add_argument("--postprocess", nargs=2, metavar=("CODE_FILE", "JOB_DIR"),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    try:
        if args.postprocess:
            from oqp_studio.postprocess_worker import run

            run(Path(args.postprocess[0]).resolve(), Path(args.postprocess[1]).resolve())
        elif args.stdio:
            serve_stdio()
        else:
            import uvicorn

            from oqp_studio.main import app

            uvicorn.run(app, host=args.host, port=args.port, log_level="info")
    except BaseException:
        # The shell can only see that the port never opened. Whatever went
        # wrong is written where it can be read afterwards -- guessing at a
        # sidecar's stack trace from "failed to start" is not debugging.
        import traceback

        report = traceback.format_exc()
        try:
            path = log_path()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(report)
        except OSError:
            pass
        print(report, file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
