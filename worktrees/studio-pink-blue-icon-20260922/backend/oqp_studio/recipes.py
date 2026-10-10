from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator

from .network import settings_path

SCHEMA = "oqp-studio-recipe/v1"
MAX_RECIPE_BYTES = 256_000
MAX_CODE_BYTES = 64_000


class PostprocessSpec(BaseModel):
    enabled: bool = False
    code: str = ""
    timeout_seconds: int = Field(default=60, ge=1, le=300)

    @field_validator("code")
    @classmethod
    def code_size(cls, value: str) -> str:
        if len(value.encode("utf-8")) > MAX_CODE_BYTES:
            raise ValueError("post-calculation Python exceeds 64 KiB")
        return value


class Recipe(BaseModel):
    schema_: Literal["oqp-studio-recipe/v1"] = Field(default=SCHEMA, alias="schema")
    id: str = Field(pattern=r"^[a-z0-9][a-z0-9._-]{1,63}$")
    name: str = Field(min_length=1, max_length=80)
    description: str = Field(default="", max_length=600)
    tags: list[str] = Field(default_factory=list, max_length=20)
    available: bool = True
    unavailable_reason: str = Field(default="", max_length=400)
    workflow: str
    method: dict[str, str | int | float | bool] = Field(default_factory=dict)
    execution: dict[str, str | int | float | bool] = Field(default_factory=dict)
    analysis: dict[str, str | int | float | bool] = Field(default_factory=dict)
    art: dict[str, str | int | float | bool] = Field(default_factory=dict)
    postprocess: PostprocessSpec = Field(default_factory=PostprocessSpec)

    model_config = {"populate_by_name": True, "extra": "forbid"}

    @field_validator("tags")
    @classmethod
    def normalize_tags(cls, values: list[str]) -> list[str]:
        return [value.strip()[:40] for value in values if value.strip()]


class RecipeRecord(BaseModel):
    recipe: Recipe
    builtin: bool = False
    trusted: bool = False


class SaveRecipeRequest(BaseModel):
    recipe: Recipe
    trust_python: bool = False
    replace: bool = False


def _builtin(
    recipe_id: str, name: str, description: str, tags: list[str], workflow: str,
    *, method: dict[str, Any] | None = None, analysis: dict[str, Any] | None = None,
    art: dict[str, Any] | None = None, available: bool = True, reason: str = "",
) -> Recipe:
    return Recipe(
        id=recipe_id, name=name, description=description, tags=tags, workflow=workflow,
        method=method or {}, analysis=analysis or {}, art=art or {}, available=available,
        unavailable_reason=reason,
    )


BUILTINS = {
    recipe.id: recipe for recipe in (
        _builtin(
            "vertical-absorption", "Vertical absorption and NTO",
            "Calculate vertical excited states, draw the absorption spectrum, and prepare NTO content for Art.",
            ["photochemistry", "absorption", "excited states", "NTO"], "abs",
            method={"theory": "mrsf", "nstate": 6},
            analysis={"spectrum_kind": "absorption", "spectrum_shape": "lorentzian", "spectrum_width": 20},
            art={"content": "excited:nto_hole", "quality": "0.75,5", "background": "studio"},
        ),
        _builtin(
            "excited-state-emission", "Excited-state relaxation and spectra",
            "Optimize S1, then prepare emission and excited-state absorption views from the relaxed structure.",
            ["photochemistry", "emission", "ESA", "excited-state optimization"], "exopt",
            method={"theory": "mrsf", "targetState": 1, "nstate": 6},
            analysis={"spectrum_kind": "emission", "spectrum_shape": "lorentzian", "spectrum_width": 20, "spectrum_state": 1},
            art={"content": "molecule", "quality": "0.75,5", "background": "neutral"},
        ),
        _builtin(
            "ip-ea-dyson", "IP/EA states and Dyson orbitals",
            "Calculate EKT ionization/electron-affinity states and make state-specific Dyson orbitals available for analysis and Art.",
            ["photoelectron", "ionization", "electron affinity", "Dyson"], "ekt",
            method={"theory": "mrsf", "nstate": 6},
            analysis={"spectrum_kind": "photoelectron", "spectrum_shape": "gaussian", "spectrum_width": 20},
            art={"content": "molecule", "quality": "0.75,5", "background": "studio"},
        ),
        _builtin(
            "vibrational-characterization", "Vibrational characterization",
            "Calculate an analytic Hessian where supported, then show frequencies, IR intensities, normal modes, and displacement arrows.",
            ["frequency", "Hessian", "IR", "normal modes"], "hess",
            method={"theory": "dft", "hessType": "analytical"},
            analysis={"spectrum_kind": "ir", "spectrum_shape": "lorentzian", "spectrum_width": 12},
            art={"content": "molecule", "quality": "0.75,5", "background": "neutral"},
        ),
        _builtin(
            "nmr-shielding-map", "NMR shielding and 3D map",
            "Calculate nuclear magnetic shielding tensors and select the total coupled shielding map for three-dimensional inspection.",
            ["NMR", "shielding", "magnetic response", "3D map"], "nmr",
            method={"theory": "dft", "nmrGauge": "giao"},
            analysis={"map_kind": "nmr"},
            art={"content": "molecule", "quality": "1,8", "background": "white"},
        ),
        _builtin(
            "nics-aromaticity", "NICS aromaticity",
            "Place magnetic shielding probes at ring centers and report NICS = -isotropic shielding.",
            ["NICS", "aromaticity", "ring current", "shielding probe"], "nmr",
            available=False,
            reason="OpenQP currently calculates shielding only at real nuclei and does not accept ghost/probe centers.",
        ),
    )
}


def directory() -> Path:
    override = os.environ.get("OQP_STUDIO_RECIPE_DIR")
    return Path(override) if override else settings_path().parent / "recipes"


def _trust_path() -> Path:
    return directory() / ".trusted-python.json"


def _code_digest(recipe: Recipe) -> str:
    content = f"{recipe.id}\0{recipe.postprocess.code}".encode()
    return hashlib.sha256(content).hexdigest()


def _trusted_digests() -> set[str]:
    try:
        values = json.loads(_trust_path().read_text())
    except (OSError, ValueError):
        return set()
    return {str(value) for value in values} if isinstance(values, list) else set()


def is_trusted(recipe: Recipe) -> bool:
    return bool(recipe.postprocess.code) and _code_digest(recipe) in _trusted_digests()


def list_records() -> list[RecipeRecord]:
    records = [RecipeRecord(recipe=recipe, builtin=True) for recipe in BUILTINS.values()]
    root = directory()
    if root.is_dir():
        for path in sorted(root.glob("*.json")):
            try:
                recipe = Recipe.model_validate_json(path.read_text())
            except (OSError, ValueError):
                continue
            records.append(RecipeRecord(recipe=recipe, trusted=is_trusted(recipe)))
    return records


def get(recipe_id: str) -> RecipeRecord | None:
    if recipe_id in BUILTINS:
        return RecipeRecord(recipe=BUILTINS[recipe_id], builtin=True)
    if not re.fullmatch(r"[a-z0-9][a-z0-9._-]{1,63}", recipe_id):
        return None
    path = directory() / f"{recipe_id}.json"
    try:
        recipe = Recipe.model_validate_json(path.read_text())
    except (OSError, ValueError):
        return None
    return RecipeRecord(recipe=recipe, trusted=is_trusted(recipe))


def save(request: SaveRecipeRequest) -> RecipeRecord:
    recipe = request.recipe
    if recipe.id in BUILTINS:
        raise ValueError("built-in recipes cannot be replaced")
    encoded = recipe.model_dump_json(by_alias=True, indent=2)
    if len(encoded.encode("utf-8")) > MAX_RECIPE_BYTES:
        raise ValueError("recipe exceeds 256 KiB")
    root = directory()
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{recipe.id}.json"
    if path.exists() and not request.replace:
        raise ValueError("a recipe with this id already exists; choose another name")
    path.write_text(encoded + "\n")
    trusted = _trusted_digests()
    digest = _code_digest(recipe)
    if request.trust_python and recipe.postprocess.enabled and recipe.postprocess.code:
        trusted.add(digest)
    else:
        trusted.discard(digest)
    _trust_path().write_text(json.dumps(sorted(trusted), indent=2) + "\n")
    return RecipeRecord(recipe=recipe, trusted=is_trusted(recipe))


def delete(recipe_id: str) -> None:
    if recipe_id in BUILTINS:
        raise ValueError("built-in recipes cannot be deleted")
    record = get(recipe_id)
    if record is None:
        raise FileNotFoundError(recipe_id)
    (directory() / f"{recipe_id}.json").unlink()


def snapshot(recipe_id: str, target: Path) -> RecipeRecord:
    record = get(recipe_id)
    if record is None:
        raise ValueError(f"unknown recipe '{recipe_id}'")
    if not record.recipe.available:
        raise ValueError(record.recipe.unavailable_reason or "recipe is unavailable")
    target.write_text(record.recipe.model_dump_json(by_alias=True, indent=2) + "\n")
    return record


def run_snapshot(snapshot_path: Path, log_path: Path) -> tuple[str | None, str | None]:
    recipe = Recipe.model_validate_json(snapshot_path.read_text())
    spec = recipe.postprocess
    if not spec.enabled or not spec.code.strip():
        return None, None
    if not is_trusted(recipe):
        return "not_trusted", "Python was not run because this recipe has not been trusted locally."

    code_path = snapshot_path.parent / ".oqp-studio-postprocess.py"
    code_path.write_text(spec.code)
    if getattr(sys, "frozen", False):
        command = [sys.executable, "--postprocess", str(code_path), str(snapshot_path.parent)]
    else:
        command = [sys.executable, "-I", str(Path(__file__).with_name("postprocess_worker.py")),
                   str(code_path), str(snapshot_path.parent)]
    try:
        completed = subprocess.run(
            command, cwd=snapshot_path.parent, text=True, capture_output=True,
            timeout=spec.timeout_seconds, env={**os.environ, "OQP_STUDIO_JOB_DIR": str(snapshot_path.parent)},
            check=False,
        )
        log_path.write_text(completed.stdout + completed.stderr)
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + (exc.stderr or "")
        log_path.write_text(output + f"\nTimed out after {spec.timeout_seconds} seconds.\n")
        return "failed", f"post-calculation Python timed out after {spec.timeout_seconds} seconds"
    finally:
        code_path.unlink(missing_ok=True)
    if completed.returncode != 0:
        return "failed", f"post-calculation Python exited with code {completed.returncode}"
    return "done", None
