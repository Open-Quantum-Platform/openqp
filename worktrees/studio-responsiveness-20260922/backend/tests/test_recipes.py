import json

from fastapi.testclient import TestClient

from oqp_studio import recipes
from oqp_studio.main import app

client = TestClient(app)


def custom_recipe(**overrides):
    values = {
        "schema": recipes.SCHEMA,
        "id": "custom-spectrum",
        "name": "Custom spectrum",
        "description": "A portable calculation and analysis recipe.",
        "tags": ["absorption", "NTO"],
        "workflow": "abs",
        "method": {"theory": "mrsf", "nstate": 6},
        "analysis": {"spectrum_kind": "absorption"},
        "art": {"background": "white"},
        "postprocess": {
            "enabled": True,
            "code": "(JOB_DIR / 'derived.txt').write_text('complete')\nprint('derived')\n",
            "timeout_seconds": 10,
        },
    }
    values.update(overrides)
    return values


def test_catalog_marks_nics_unavailable_until_openqp_accepts_probe_centers(tmp_path, monkeypatch):
    monkeypatch.setenv("OQP_STUDIO_RECIPE_DIR", str(tmp_path))
    records = client.get("/api/recipes").json()
    nics = next(record for record in records if record["recipe"]["id"] == "nics-aromaticity")

    assert nics["recipe"]["available"] is False
    assert "ghost/probe" in nics["recipe"]["unavailable_reason"]


def test_imported_python_is_not_trusted_by_recipe_contents(tmp_path, monkeypatch):
    monkeypatch.setenv("OQP_STUDIO_RECIPE_DIR", str(tmp_path))
    payload = custom_recipe()
    response = client.post("/api/recipes", json={"recipe": payload, "trust_python": False})

    assert response.status_code == 200
    assert response.json()["trusted"] is False
    saved = json.loads((tmp_path / "custom-spectrum.json").read_text())
    assert "trusted" not in saved

    duplicate = client.post(
        "/api/recipes", json={"recipe": payload, "trust_python": True, "replace": False}
    )
    assert duplicate.status_code == 400
    assert "already exists" in duplicate.json()["detail"]

    trusted = client.post(
        "/api/recipes", json={"recipe": payload, "trust_python": True, "replace": True}
    )
    assert trusted.status_code == 200
    assert trusted.json()["trusted"] is True


def test_trusted_python_runs_from_a_snapshot_and_records_output(tmp_path, monkeypatch):
    recipe_dir = tmp_path / "recipes"
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    monkeypatch.setenv("OQP_STUDIO_RECIPE_DIR", str(recipe_dir))
    record = recipes.save(recipes.SaveRecipeRequest(
        recipe=recipes.Recipe.model_validate(custom_recipe()), trust_python=True,
    ))
    assert record.trusted is True
    snapshot = job_dir / ".oqp-studio-recipe.json"
    recipes.snapshot(record.recipe.id, snapshot)

    status, error = recipes.run_snapshot(snapshot, job_dir / "postprocess.log")

    assert status == "done"
    assert error is None
    assert (job_dir / "derived.txt").read_text() == "complete"
    assert (job_dir / "postprocess.log").read_text() == "derived\n"


def test_untrusted_snapshot_never_executes_python(tmp_path, monkeypatch):
    recipe_dir = tmp_path / "recipes"
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    monkeypatch.setenv("OQP_STUDIO_RECIPE_DIR", str(recipe_dir))
    record = recipes.save(recipes.SaveRecipeRequest(
        recipe=recipes.Recipe.model_validate(custom_recipe()), trust_python=False,
    ))
    snapshot = job_dir / ".oqp-studio-recipe.json"
    recipes.snapshot(record.recipe.id, snapshot)

    status, error = recipes.run_snapshot(snapshot, job_dir / "postprocess.log")

    assert status == "not_trusted"
    assert "not been trusted" in error
    assert not (job_dir / "derived.txt").exists()


def test_job_keeps_calculation_success_when_trusted_postprocessing_finishes(
    tmp_path, monkeypatch,
):
    from oqp_studio import jobs

    recipe_dir = tmp_path / "recipes"
    job_root = tmp_path / "jobs"
    monkeypatch.setenv("OQP_STUDIO_RECIPE_DIR", str(recipe_dir))
    monkeypatch.setattr(jobs, "JOBS_ROOT", job_root)
    recipes.save(recipes.SaveRecipeRequest(
        recipe=recipes.Recipe.model_validate(custom_recipe()), trust_python=True,
    ))

    class Runner:
        def run(self, *_args, **_kwargs):
            return 0

    monkeypatch.setattr(jobs, "get_runner", lambda _name: Runner())
    manager = jobs.JobManager()
    manager._ready = True
    info = manager._prepare(jobs.JobRequest(
        input_text="hf/sto-3g\nenergy\n", recipe_id="custom-spectrum",
    ))

    manager._run(info.id)

    assert info.status == jobs.JobStatus.done
    assert info.postprocess_status == "done"
    assert info.postprocess_error is None
    assert (job_root / info.id / "derived.txt").read_text() == "complete"
