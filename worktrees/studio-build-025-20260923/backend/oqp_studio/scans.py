"""Generate serial relaxed or rigid bond-distance scans."""

from __future__ import annotations

import math
import re
import uuid

from pydantic import BaseModel, Field

from .jobs import JobRequest
from .structure_io import parse_oqp


class BondScanRequest(BaseModel):
    input_text: str
    input_name: str | None = None
    name: str = "bond-scan"
    runner: str = "local"
    threads: int = Field(default=1, ge=1)
    atom_a: int = Field(ge=1)
    atom_b: int = Field(ge=1)
    start: float = Field(gt=0)
    end: float = Field(gt=0)
    points: int = Field(default=7, ge=2, le=101)
    relaxed: bool = True
    recipe_id: str | None = None


def _geometry_text(atoms: list[tuple[str, float, float, float]]) -> str:
    rows = [
        f"{symbol:<2} {x:11.6f} {y:11.6f} {z:11.6f}"
        for symbol, x, y, z in atoms
    ]
    return "\n".join(rows)


def _replace_geometry(text: str, atoms: list[tuple[str, float, float, float]]) -> str:
    pattern = re.compile(
        r'(geom(?:etry)?\s*=\s*""")(.*?)(""")', re.IGNORECASE | re.DOTALL
    )
    if not pattern.search(text):
        raise ValueError("bond scans require an inline triple-quoted geom block")
    return pattern.sub(lambda match: f'{match.group(1)}\n{_geometry_text(atoms)}\n{match.group(3)}',
                       text, count=1)


def _ensure_distance_constraint(text: str, atom_a: int, atom_b: int) -> str:
    lines = text.splitlines()
    driver_index = next((index for index, line in enumerate(lines[1:], start=1)
                         if line.strip() and not line.lstrip().startswith("#")), None)
    if driver_index is None or not re.match(
        r"\s*opt\b", lines[driver_index], re.IGNORECASE
    ):
        raise ValueError("a relaxed bond scan requires an opt driver")
    line = lines[driver_index].strip()
    freeze = f"freeze=distance({atom_a},{atom_b})"
    if re.search(r"\bfreeze\s*=", line, re.IGNORECASE):
        if not re.search(r"\bfreeze\s*=\s*distance\([^)]*\)", line, re.IGNORECASE):
            raise ValueError(
                "a relaxed bond scan cannot replace an existing non-distance freeze"
            )
        line = re.sub(r"freeze\s*=\s*distance\([^)]*\)", freeze, line,
                      count=1, flags=re.IGNORECASE)
    elif line.endswith(")"):
        line = f"{line[:-1]},{freeze})"
    else:
        line = f"{line}({freeze})"
    lines[driver_index] = line
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def _ensure_rigid_driver(text: str) -> str:
    lines = text.splitlines()
    driver_index = next((index for index, line in enumerate(lines[1:], start=1)
                         if line.strip() and not line.lstrip().startswith("#")), None)
    if driver_index is None:
        raise ValueError("a rigid bond scan requires an energy driver")
    line = lines[driver_index].strip()
    if re.match(r"energy\b", line, re.IGNORECASE):
        return text
    if not re.match(r"opt\b", line, re.IGNORECASE):
        raise ValueError("a rigid bond scan requires an energy driver")
    state = re.search(r"\(\s*(S\d+)\b", line, re.IGNORECASE)
    lines[driver_index] = f"energy({state.group(1)})" if state else "energy"
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def target_state(text: str) -> int | None:
    match = re.search(r"\b(?:energy|opt)\s*\(\s*S(\d+)\b", text, re.IGNORECASE)
    return int(match.group(1)) if match else None


def build(request: BondScanRequest) -> tuple[str, list[JobRequest], list[float]]:
    if request.atom_a == request.atom_b:
        raise ValueError("scan atoms must differ")
    frames = parse_oqp(request.input_text)
    if not frames or not frames[0].atoms:
        raise ValueError("no geometry was found in the scan input")
    source = list(frames[0].atoms)
    if max(request.atom_a, request.atom_b) > len(source):
        raise ValueError(f"scan atom index exceeds the {len(source)}-atom geometry")
    a = request.atom_a - 1
    b = request.atom_b - 1
    ax, ay, az = source[a][1:]
    bx, by, bz = source[b][1:]
    dx, dy, dz = bx - ax, by - ay, bz - az
    length = math.sqrt(dx * dx + dy * dy + dz * dz)
    if length < 1.0e-10:
        raise ValueError("scan atoms occupy the same position")
    direction = (dx / length, dy / length, dz / length)
    values = [
        request.start + (request.end - request.start) * index / (request.points - 1)
        for index in range(request.points)
    ]
    group_id = uuid.uuid4().hex[:12]
    jobs: list[JobRequest] = []
    for index, value in enumerate(values, start=1):
        atoms = list(source)
        symbol = atoms[b][0]
        atoms[b] = (
            symbol,
            ax + direction[0] * value,
            ay + direction[1] * value,
            az + direction[2] * value,
        )
        input_text = _replace_geometry(request.input_text, atoms)
        if request.relaxed:
            input_text = _ensure_distance_constraint(input_text, request.atom_a, request.atom_b)
        else:
            input_text = _ensure_rigid_driver(input_text)
        jobs.append(JobRequest(
            input_text=input_text,
            input_name=request.input_name,
            name=f"{request.name} {index}/{request.points} ({value:.4f} A)",
            runner=request.runner,
            threads=request.threads,
            recipe_id=request.recipe_id,
        ))
    return group_id, jobs, values
