import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import ClassVar

import pytest
from fastapi.testclient import TestClient

from oqp_studio.main import app

client = TestClient(app)


def test_health():
    res = client.get("/api/health")
    assert res.status_code == 200
    assert res.json()["status"] == "ok"


def test_live_log_prefers_the_openqp_calculation_log(tmp_path, monkeypatch):
    """The Execution panel must show OpenQP's record, not launcher stdout."""
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    job_dir = tmp_path / "calculation"
    job_dir.mkdir()
    (job_dir / "input.oqp").write_text("hf/sto-3g\nenergy\n")
    (job_dir / "job.log").write_text("launcher progress\n")
    (job_dir / "input.log").write_text("TOTAL energy = -76.01074651\n")

    assert jobs.manager.log_tail("calculation") == "TOTAL energy = -76.01074651\n"


def test_custom_input_name_controls_the_calculation_log(tmp_path, monkeypatch):
    """The record follows the selected input name, not a fixed oqp.log name."""
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    job_dir = tmp_path / "calculation"
    job_dir.mkdir()
    (job_dir / "water.oqp").write_text("hf/sto-3g\nenergy\n")
    (job_dir / "water.log").write_text("TOTAL energy = -76.01074651\n")
    (job_dir / "job.log").write_text("launcher progress\n")

    assert jobs.manager.log_tail("calculation") == "TOTAL energy = -76.01074651\n"


def test_custom_input_name_is_saved(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    res = client.post(
        "/api/jobs",
        json={
            "input_text": "hf/sto-3g\nenergy\n",
            "input_name": "water.oqp",
            "name": "water",
            "runner": "local",
        },
    )
    assert res.status_code == 200
    assert res.json()["name"] == "water"
    names = {f["name"] for f in client.get(f"/api/jobs/{res.json()['id']}/files").json()}
    assert "water.oqp" in names


def test_recovered_jobs_keep_their_project_name(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    job_dir = tmp_path / "opaque-id"
    job_dir.mkdir()
    (job_dir / "water.oqp").write_text("hf/sto-3g\nenergy\n")
    (job_dir / ".oqp-studio.json").write_text(
        '{"name":"water TS candidate","runner":"bundled","threads":4,"created_at":"2026-08-22T00:00:00+00:00"}'
    )

    manager = jobs.JobManager()
    manager._ready = True
    manager._recover()

    info = manager.get("opaque-id")
    assert info is not None
    assert info.name == "water TS candidate"
    assert info.runner == "bundled"
    assert info.threads == 4


def test_recovered_jobs_fall_back_to_the_input_stem(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    job_dir = tmp_path / "opaque-id"
    job_dir.mkdir()
    (job_dir / "ammonia.oqp").write_text("hf/sto-3g\nenergy\n")

    manager = jobs.JobManager()
    manager._ready = True
    manager._recover()

    assert manager.get("opaque-id").name == "ammonia"


def test_native_optimizer_nonconvergence_is_not_reported_as_done(tmp_path):
    from oqp_studio.jobs import JobManager

    job_dir = tmp_path / "water"
    job_dir.mkdir()
    (job_dir / "water.oqp").write_text("hf/sto-3g\nopt\n")
    (job_dir / "water.log").write_text(
        "PyOQP: Native optimization after internal recovery did not converge "
        "(RMS gradient 4.916e-02). The best geometry was retained.\n"
    )

    diagnostic = JobManager._optimization_diagnostic(job_dir)

    assert diagnostic is not None
    assert "did not converge" in diagnostic.lower()


def test_restart_input_replaces_inline_geometry_with_retained_optimum(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    job_dir = tmp_path / "water"
    job_dir.mkdir()
    (job_dir / "water.oqp").write_text(
        "hf/sto-3g\nopt\ngeom=\"\"\"\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n\"\"\"\n"
    )
    (job_dir / "opt.xyz").write_text(
        "2\nretained optimum\nO 0.100000 0.200000 0.300000\nH 0.400000 0.500000 0.600000\n"
    )
    manager._jobs["water"] = jobs.JobInfo(
        id="water", name="water", status=jobs.JobStatus.not_converged,
        runner="local", threads=3, created_at="2026-08-22T00:00:00+00:00",
    )
    restart = manager.restart_input("water")

    assert restart["runner"] == "local"
    assert restart["threads"] == 3
    assert 'geom="""\nO     0.100000    0.200000    0.300000' in restart["input_text"]
    assert "H     0.400000    0.500000    0.600000" in restart["input_text"]


def test_cancel_terminates_the_running_process_group(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    job_dir = tmp_path / "slow"
    job_dir.mkdir()
    (job_dir / "input.oqp").write_text("hf/sto-3g\nenergy\n")
    manager._jobs["slow"] = jobs.JobInfo(
        id="slow", name="slow", status=jobs.JobStatus.queued,
        runner="local", created_at="2026-08-22T00:00:00+00:00",
    )

    class SlowRunner:
        def run(self, _job_dir, _threads, on_start):
            process = subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                start_new_session=True,
            )
            on_start(process)
            return process.wait()

    monkeypatch.setattr(jobs, "get_runner", lambda _name: SlowRunner())
    worker = threading.Thread(target=manager._run, args=("slow",))
    worker.start()
    for _ in range(100):
        if "slow" in manager._processes:
            break
        time.sleep(0.01)
    manager.cancel("slow")
    worker.join(timeout=3)

    assert not worker.is_alive()
    assert manager.get("slow").status == jobs.JobStatus.cancelled


def test_cancelling_one_scan_point_cancels_the_queued_group(tmp_path, monkeypatch):
    import json
    from datetime import datetime, timezone

    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    for index in range(3):
        job_id = f"point-{index}"
        (tmp_path / job_id).mkdir()
        manager._jobs[job_id] = jobs.JobInfo(
            id=job_id, name=job_id, status=jobs.JobStatus.queued, runner="bundled",
            created_at=datetime.now(timezone.utc).isoformat(), group_id="scan",
            scan_value=float(index), scan_unit="A",
        )

    manager.cancel("point-0")

    assert {info.status for info in manager.list()} == {jobs.JobStatus.cancelled}
    assert all(json.loads((tmp_path / info.id / ".oqp-studio.json").read_text())["status"]
               == "cancelled" for info in manager.list())


def test_cancel_keeps_a_starting_process_in_cancelling_state(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    job_dir = tmp_path / "starting"
    job_dir.mkdir()
    (job_dir / "input.oqp").write_text("hf/sto-3g\nenergy\n")
    manager._jobs["starting"] = jobs.JobInfo(
        id="starting", name="starting", status=jobs.JobStatus.queued,
        runner="local", created_at="2026-08-22T00:00:00+00:00",
    )
    runner_entered = threading.Event()
    allow_start = threading.Event()

    class StartingRunner:
        def run(self, _job_dir, _threads, on_start):
            runner_entered.set()
            allow_start.wait(timeout=3)
            process = subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                start_new_session=True,
            )
            on_start(process)
            return process.wait()

    monkeypatch.setattr(jobs, "get_runner", lambda _name: StartingRunner())
    worker = threading.Thread(target=manager._run, args=("starting",))
    worker.start()
    assert runner_entered.wait(timeout=3)

    manager.cancel("starting")
    assert manager.get("starting").status == jobs.JobStatus.cancelling
    allow_start.set()
    worker.join(timeout=3)

    assert not worker.is_alive()
    assert manager.get("starting").status == jobs.JobStatus.cancelled


def test_batch_preparation_failure_removes_earlier_points(tmp_path, monkeypatch):
    import pytest

    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    monkeypatch.setattr(manager, "_validate_request", lambda _request: None)
    original_prepare = manager._prepare
    calls = 0

    def prepare(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("disk full")
        return original_prepare(*args, **kwargs)

    monkeypatch.setattr(manager, "_prepare", prepare)
    requests = [jobs.JobRequest(input_text="hf/sto-3g\nenergy\n") for _ in range(3)]

    with pytest.raises(OSError, match="disk full"):
        manager.submit_batch(requests, group_id="scan", values=[1.0, 1.1, 1.2], unit="A")

    assert manager.list() == []
    assert list(tmp_path.iterdir()) == []


def test_prepare_rolls_back_the_current_point_when_metadata_write_fails(tmp_path, monkeypatch):
    import pytest

    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    def fail_write(_path, _info):
        raise OSError("disk full")

    monkeypatch.setattr(manager, "_write_metadata", fail_write)

    with pytest.raises(OSError, match="disk full"):
        manager._prepare(jobs.JobRequest(input_text="hf/sto-3g\nenergy\n"))

    assert manager.list() == []
    assert list(tmp_path.iterdir()) == []


def test_running_metadata_failure_leaves_a_terminal_failed_job(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    job_dir = tmp_path / "metadata-failure"
    job_dir.mkdir()
    manager._jobs["metadata-failure"] = jobs.JobInfo(
        id="metadata-failure", name="metadata failure", status=jobs.JobStatus.queued,
        runner="local", created_at="2026-08-22T00:00:00+00:00",
    )
    runner_called = False

    class Runner:
        def run(self, *_args, **_kwargs):
            nonlocal runner_called
            runner_called = True
            return 0

    monkeypatch.setattr(jobs, "get_runner", lambda _name: Runner())
    def fail_write(_path, _info):
        raise OSError("read only")

    monkeypatch.setattr(manager, "_write_metadata", fail_write)

    manager._run("metadata-failure")

    assert not runner_called
    assert manager.get("metadata-failure").status == jobs.JobStatus.failed
    assert manager.get("metadata-failure").error == "read only"


def test_recovery_does_not_mark_an_unstarted_scan_point_done(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    point = tmp_path / "queued-point"
    point.mkdir()
    (point / "point.oqp").write_text("hf/sto-3g\nenergy\n")
    (point / ".oqp-studio.json").write_text(
        '{"name":"queued point","runner":"bundled","threads":1,'
        '"created_at":"2026-08-22T00:00:00+00:00","status":"queued",'
        '"group_id":"scan","scan_value":0.9,"scan_unit":"A"}'
    )

    manager = jobs.JobManager()
    manager._ready = True
    manager._recover()

    recovered = manager.get("queued-point")
    assert recovered.status == jobs.JobStatus.cancelled
    assert "Interrupted" in recovered.error


def test_completed_project_can_be_deleted_with_its_result_files(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    adopted = manager.adopt("water", [("water.log", b"TOTAL energy = -76.0\n")])

    manager.delete(adopted.id)

    assert manager.get(adopted.id) is None
    assert not (tmp_path / adopted.id).exists()


def test_qmmm_pdb_asset_is_saved_with_its_input(tmp_path, monkeypatch):
    """The PDB referenced by a QM/MM route must reach the runner's directory."""
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    res = client.post(
        "/api/jobs",
        json={
            "input_text": 'hf/sto-3g qmmm_flag=true\ngeom="system.pdb 0"\n',
            "input_name": "system.oqp",
            "pdb_name": "system.pdb",
            "pdb_text": "ATOM      1  O   HOH A   1       0.0 0.0 0.0\nEND\n",
        },
    )
    assert res.status_code == 200
    job_id = res.json()["id"]
    assert (tmp_path / job_id / "system.pdb").read_text().endswith("END\n")


def test_runners_listed():
    res = client.get("/api/runners")
    assert res.status_code == 200
    assert set(res.json()) >= {"local", "bundled", "wsl"}


def test_memory_admission_refuses_a_calculation_above_available_ram(monkeypatch):
    """The execution preflight must stop before a memory-starved launch."""
    from oqp_studio import host

    monkeypatch.setattr(host, "snapshot", lambda: {
        "platform": "darwin",
        "physical_cores": 4,
        "logical_cores": 8,
        "memory_total_bytes": 8 * 1024**3,
        "memory_available_bytes": 512 * 1024**2,
    })
    check = host.admission("ccsd/cc-pvtz\nenergy\ngeom=\"\"\"\nO 0 0 0\n\"\"\"\n", 4)
    assert not check["permitted"]
    assert "currently available RAM" in check["reason"]


def test_bundled_engine_version_is_read_from_its_readme(tmp_path):
    from oqp_studio import engine

    executable = tmp_path / engine.EXECUTABLE
    executable.write_text("#!/bin/sh\n")
    (tmp_path / "README.txt").write_text("OpenQP version : 1.3.1\n")

    assert engine.version(str(executable)) == "1.3.1"


def test_submit_without_openqp_fails_gracefully(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    res = client.post(
        "/api/jobs",
        json={"input_text": "[input]\nsystem=O 0 0 0\n", "runner": "local"},
    )
    assert res.status_code == 200
    job_id = res.json()["id"]

    import time

    for _ in range(50):
        info = client.get(f"/api/jobs/{job_id}").json()
        if info["status"] in ("done", "failed"):
            break
        time.sleep(0.1)
    # Without an OpenQP install the job must fail with a clear error,
    # never hang or crash the server.
    assert info["status"] in ("done", "failed")

    files = client.get(f"/api/jobs/{job_id}/files").json()
    names = {f["name"] for f in files}
    assert "input.inp" in names

    res = client.get(f"/api/jobs/{job_id}/files/input.inp")
    assert res.status_code == 200
    assert "system=O 0 0 0" in res.text

    # Path escapes must 404, not leak files outside the job directory.
    res = client.get(f"/api/jobs/{job_id}/files/..%2F..%2Fpyproject.toml")
    assert res.status_code == 404


def test_oqp_style_input_saved_as_oqp(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    res = client.post(
        "/api/jobs",
        json={"input_text": 'hf/sto-3g\nenergy\ngeom="""\nO 0 0 0\n"""\n', "runner": "local"},
    )
    job_id = res.json()["id"]
    names = {f["name"] for f in client.get(f"/api/jobs/{job_id}/files").json()}
    assert "input.oqp" in names


def test_structure3d():
    res = client.post("/api/structure3d", json={"smiles": "O"})
    try:
        import rdkit  # noqa: F401
    except ImportError:
        assert res.status_code == 501
        return
    assert res.status_code == 200
    atoms = res.json()["atoms"]
    assert len(atoms) == 3  # water with added hydrogens
    assert sorted(a[0] for a in atoms) == ["H", "H", "O"]


def test_molden_parse_and_cube():
    from pathlib import Path

    from oqp_studio import molden

    sample = (
        Path(__file__).resolve().parents[2]
        / "frontend" / "public" / "viewer" / "samples" / "water-hessian-mo.molden"
    )
    data = molden.parse_molden(sample.read_text())
    assert data.supported
    assert len(data.atoms) == 3
    assert data.orbitals, "no MOs parsed"
    assert len(data.basis) == len(data.orbitals[0].coefficients)

    cube = molden.orbital_cube(data, mo_index=1, max_points=40_000)
    header = cube.splitlines()
    assert int(header[2].split()[0]) == 3  # atom count
    assert "E=" in header[1]


def test_excited_state_analysis_uses_the_physical_root_transition_density():
    """NTOs come from S0->S1, not the auxiliary-reference response vector."""
    import numpy as np

    from oqp_studio import excited_state, molden

    densities = np.zeros((2, 2, 2, 2))
    densities[1, 0, 0, 1] = 1.0  # particle MO 2, hole MO 1
    densities[:, :, 1, 1] = np.diag([-0.7, 0.7])
    data = excited_state.ExcitedStateData(
        energies=np.array([-0.2, 0.1]),
        transition_densities=densities,
        coefficients=np.eye(2),
        reference_occupations=np.array([2.0, 0.0]),
        molden_data=molden.MoldenData(atoms=[]),
        json_name="mrsf.json",
        molden_name="mrsf.molden",
    )

    report = data.summary(0, 1)
    nto = data.nto(0, 1)

    assert report["nto_pairs"][0]["fraction"] == 1.0
    assert report["nto_participation_ratio"] == 1.0
    assert abs(report["n_promoted"] - 0.7) < 1e-12
    np.testing.assert_allclose(np.abs(nto["holes_ao"][:, 0]), [1.0, 0.0])
    np.testing.assert_allclose(np.abs(nto["particles_ao"][:, 0]), [0.0, 1.0])


def test_excited_state_json_restores_openqp_fortran_array_order(tmp_path):
    import json

    import numpy as np

    from oqp_studio import excited_state, molden

    sample = (
        Path(__file__).resolve().parents[2]
        / "frontend" / "public" / "viewer" / "samples" / "water-hessian-mo.molden"
    )
    nbf = len(molden.parse_molden(sample.read_text()).basis)
    densities = np.zeros((nbf, nbf, 2, 2))
    densities[2, 1, 0, 1] = 0.625
    densities[0, 0, 1, 1] = -0.4
    densities[1, 1, 1, 1] = 0.4
    result = tmp_path / "mrsf.json"
    result.write_text(json.dumps({
        "OQP::td_energies": [-10.0, -9.7],
        "OQP::VEC_MO_A": np.eye(nbf).ravel().tolist(),
        "OQP::td_trans_density_mo": densities.ravel(order="F").tolist(),
    }))

    loaded = excited_state.ExcitedStateData.load([result, sample])

    assert loaded.tdm_mo(0, 1)[2, 1] == 0.625
    np.testing.assert_allclose(loaded.difference_mo(0, 1)[:2, :2], [[-0.4, 0], [0, 0.4]])


def test_arbitrary_ao_density_matrix_can_be_exported_as_cube():
    import numpy as np

    from oqp_studio import molden

    sample = (
        Path(__file__).resolve().parents[2]
        / "frontend" / "public" / "viewer" / "samples" / "water-hessian-mo.molden"
    )
    data = molden.parse_molden(sample.read_text())
    cube = molden.matrix_density_cube(
        data, np.eye(len(data.basis)) * 0.01,
        "test excited density", "density", max_points=40_000,
    )

    assert cube.startswith("test excited density\ndensity\n")
    assert int(cube.splitlines()[2].split()[0]) == 3


def test_relaxed_bond_scan_generates_target_geometries_and_native_constraints():
    from oqp_studio import scans
    from oqp_studio.structure_io import parse_oqp

    request = scans.BondScanRequest(
        input_text=(
            "dft/b3lyp/6-31g*\nopt(maxit=20)\ngeom=\"\"\"\n"
            "H 0 0 0\nH 0 0 0.7\n\"\"\"\n"
        ),
        atom_a=1, atom_b=2, start=0.8, end=1.2, points=3, relaxed=True,
    )

    group_id, jobs, values = scans.build(request)

    assert len(group_id) == 12
    assert values == [0.8, 1.0, 1.2]
    assert len(jobs) == 3
    for job, target in zip(jobs, values):
        atoms = parse_oqp(job.input_text)[0].atoms
        assert abs(atoms[1][3] - target) < 1.0e-12
        assert job.input_text.count("freeze=distance(1,2)") == 1


def test_relaxed_bond_scan_rejects_an_existing_non_distance_freeze():
    import pytest

    from oqp_studio import scans

    request = scans.BondScanRequest(
        input_text=(
            "dft/b3lyp/6-31g*\nopt(freeze=angle(1,2,3))\ngeom=\"\"\"\n"
            "H 0 0 0\nO 0 0 1\nH 1 0 1\n\"\"\"\n"
        ),
        atom_a=1, atom_b=2, start=0.8, end=1.2, points=3, relaxed=True,
    )

    with pytest.raises(ValueError, match="existing non-distance freeze"):
        scans.build(request)


def test_bond_scan_rejects_an_atom_outside_the_geometry():
    import pytest

    from oqp_studio import scans

    request = scans.BondScanRequest(
        input_text="hf/sto-3g\nenergy\ngeom=\"\"\"\nH 0 0 0\n\"\"\"\n",
        atom_a=1, atom_b=2, start=0.8, end=1.0,
    )
    with pytest.raises(ValueError, match="exceeds the 1-atom geometry"):
        scans.build(request)


def test_bond_scan_api_submits_one_persisted_job_per_point(monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    class ScanManager:
        def __init__(self):
            self.infos = []

        def submit_batch(self, requests, *, group_id, values, unit, state=None):
            self.infos = [
                jobs.JobInfo(
                    id=f"point-{index}", name=request.name, status=jobs.JobStatus.queued,
                    runner=request.runner, threads=request.threads,
                    created_at=datetime.now(timezone.utc).isoformat(),
                    group_id=group_id, scan_value=value, scan_unit=unit, scan_state=state,
                )
                for index, (request, value) in enumerate(zip(requests, values), start=1)
            ]
            return self.infos

        def list(self):
            return self.infos

    scan_manager = ScanManager()
    monkeypatch.setattr(main, "manager", scan_manager)
    response = client.post("/api/scans", json={
        "input_text": "hf/sto-3g\nenergy\ngeom=\"\"\"\nH 0 0 0\nH 0 0 0.7\n\"\"\"\n",
        "name": "H2 stretch", "atom_a": 1, "atom_b": 2,
        "start": 0.7, "end": 1.1, "points": 3, "relaxed": False,
    })

    assert response.status_code == 200, response.text
    payload = response.json()
    assert len(payload["jobs"]) == 3
    assert {job["group_id"] for job in payload["jobs"]} == {payload["group_id"]}
    status = client.get(f"/api/scans/{payload['group_id']}")
    assert [point["value"] for point in status.json()["points"]] == [0.7, 0.9, 1.1]


def test_bond_scan_api_reads_total_energy_from_an_openqp_log(tmp_path, monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    job_dir = tmp_path / "scan-point"
    job_dir.mkdir()
    (job_dir / "point.oqp").write_text("hf/sto-3g\nenergy\n")
    (job_dir / "point.log").write_text("TOTAL energy = -1.12000000\n")
    manager = jobs.JobManager()
    manager._ready = True
    manager._jobs["scan-point"] = jobs.JobInfo(
        id="scan-point", name="H2 0.9 A", status=jobs.JobStatus.done,
        runner="bundled", created_at=datetime.now(timezone.utc).isoformat(),
        group_id="energy-scan", scan_value=0.9, scan_unit="A",
    )
    monkeypatch.setattr(main, "manager", manager)

    response = client.get("/api/scans/energy-scan")

    assert response.status_code == 200
    assert response.json()["points"][0]["energy"] == -1.12


def test_rigid_scan_replaces_an_optimization_driver_with_energy():
    from oqp_studio import scans

    request = scans.BondScanRequest(
        input_text=("mrsf(nstate=3)/bhhlyp/6-31g*\nopt(S1,maxit=20)\n"
                    "geom=\"\"\"\nH 0 0 0\nH 0 0 0.7\n\"\"\"\n"),
        atom_a=1, atom_b=2, start=0.8, end=1.0, points=2, relaxed=False,
    )

    _, jobs, _ = scans.build(request)

    assert all("\nenergy(S1)\n" in job.input_text for job in jobs)
    assert all("\nopt" not in job.input_text for job in jobs)
    assert scans.target_state(request.input_text) == 1


def test_scan_endpoint_uses_the_requested_excited_state_energy(monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    class ExcitedScanManager:
        def list(self):
            return [jobs.JobInfo(
                id="excited", name="S1 scan", status=jobs.JobStatus.done,
                runner="bundled", created_at=datetime.now(timezone.utc).isoformat(),
                group_id="s1-scan", scan_value=1.0, scan_unit="A", scan_state=1,
            )]

        def get(self, job_id):
            return self.list()[0] if job_id == "excited" else None

    monkeypatch.setattr(main, "manager", ExcitedScanManager())
    monkeypatch.setattr(main, "_job_summary", lambda _job_id: {
        "energy": {"total": -76.3},
        "states": [{"index": 0, "total": -76.3}, {"index": 1, "total": -76.0}],
    })
    monkeypatch.setattr(main, "_job_paths", lambda _job_id: [])
    main._scan_energy_cache.clear()

    response = client.get("/api/scans/s1-scan")

    assert response.status_code == 200
    assert response.json()["points"][0]["energy"] == -76.0


def test_relaxed_scan_uses_the_final_optimized_state_energy(monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import analysis, jobs, main

    info = jobs.JobInfo(
        id="relaxed", name="S1 relaxed", status=jobs.JobStatus.done,
        runner="bundled", created_at=datetime.now(timezone.utc).isoformat(),
        group_id="relaxed-scan", scan_value=1.0, scan_unit="A", scan_state=1,
    )

    class RelaxedManager:
        def list(self):
            return [info]

        def get(self, job_id):
            return info if job_id == info.id else None

    monkeypatch.setattr(main, "manager", RelaxedManager())
    monkeypatch.setattr(main, "_job_summary", lambda _job_id: {
        "energy": {"total": -76.3},
        "states": [{"index": 1, "total": -76.0}],
    })
    monkeypatch.setattr(main, "_job_paths", lambda _job_id: [])
    monkeypatch.setattr(analysis, "optimization_history", lambda _paths: {"steps": [
        {"index": 1, "states": [{"index": 1, "total": -76.0}]},
        {"index": 2, "states": [{"index": 1, "total": -76.2}]},
    ]})
    main._scan_energy_cache.clear()

    response = client.get("/api/scans/relaxed-scan")

    assert response.status_code == 200
    assert response.json()["points"][0]["energy"] == -76.2


def test_scan_polling_parses_each_completed_point_only_once(monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    infos = [jobs.JobInfo(
        id=f"point-{index}", name=f"point {index}", status=jobs.JobStatus.done,
        runner="bundled", created_at=datetime.now(timezone.utc).isoformat(),
        group_id="long-scan", scan_value=float(index), scan_unit="A",
    ) for index in range(12)]

    class LongScanManager:
        def list(self):
            return infos

    calls = 0

    def summary(_job_id):
        nonlocal calls
        calls += 1
        return {"energy": {"total": -float(calls)}}

    monkeypatch.setattr(main, "manager", LongScanManager())
    monkeypatch.setattr(main, "_job_summary", summary)
    main._scan_energy_cache.clear()

    assert client.get("/api/scans/long-scan").status_code == 200
    assert client.get("/api/scans/long-scan").status_code == 200
    assert calls == len(infos)


def test_project_comparison_aligns_structures_and_keeps_energy_differences(tmp_path, monkeypatch):
    import json
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    for job_id, name, energy, xyz, states in (
        ("left", "reference", -1.1, "H 0 0 -0.4\nH 0 0 0.4\n", [-1.1, -0.9]),
        ("right", "candidate", -1.0, "H 4 -0.4 2\nH 4 0.4 2\n", [-1.0, -0.75]),
    ):
        directory = tmp_path / job_id
        directory.mkdir()
        (directory / "result.xyz").write_text(f"2\n{name}\n{xyz}")
        (directory / "result.json").write_text(json.dumps({
            "energy": energy, "td_energies": states, "dipole": [0.0, 0.0, 1.0],
        }))
        manager._jobs[job_id] = jobs.JobInfo(
            id=job_id, name=name, status=jobs.JobStatus.done, runner="bundled",
            created_at=datetime.now(timezone.utc).isoformat(),
        )
    monkeypatch.setattr(main, "manager", manager)
    main._summary_cache.clear()
    main._summary_cache.update({
        "left": {"energy": {"total": 10.0}},
        "right": {"energy": {"total": 20.0}},
    })

    response = client.get("/api/comparison?left=left&right=right")

    assert response.status_code == 200, response.text
    data = response.json()
    assert abs(data["energy_delta_hartree"] - 0.1) < 1.0e-12
    assert abs(data["energy_delta_kcal_mol"] - 62.7509474) < 1.0e-9
    assert data["geometry"]["rmsd_angstrom"] < 1.0e-12
    assert data["states"][1]["delta_ev"] > 1.0


def test_project_comparison_uses_the_optimized_excited_state_energy():
    from oqp_studio import main

    summary = {
        "excited_state_optimized": 1,
        "energy": {
            "components": {"total": -76.3},
            "final_states": {1: -76.05},
        },
        "scf": {"energy": -76.3},
    }

    assert main._summary_energy(summary) == -76.05


def test_comparison_geometry_skips_cube_files(tmp_path, monkeypatch):
    from oqp_studio import main, structure_io

    cube = tmp_path / "density.cube"
    cube.write_text("cube data that resembles coordinates")
    monkeypatch.setattr(structure_io, "parse", lambda *_args, **_kwargs: None)

    assert main._comparison_frame([cube]) is None


def test_comparison_geometry_includes_supported_import_formats(tmp_path, monkeypatch):
    from oqp_studio import main, structure_io

    structure = tmp_path / "reference.sdf"
    structure.write_text("structure")
    monkeypatch.setattr(
        structure_io,
        "parse",
        lambda *_args, **_kwargs: structure_io.Structure(
            "sdf", [structure_io.Frame([("C", 0.0, 0.0, 0.0)])],
        ),
    )

    assert main._comparison_frame([structure]).atoms[0][0] == "C"


def test_comparison_geometry_prefers_openqp_optimization_trajectory(tmp_path, monkeypatch):
    from oqp_studio import main, structure_io

    trajectory = tmp_path / "opt_geom.xyz"
    log = tmp_path / "calculation.log"
    trajectory.write_text("trajectory")
    log.write_text("log")

    def parse(name, *_args, **_kwargs):
        x = 2.0 if name == trajectory.name else 1.0
        return structure_io.Structure("test", [structure_io.Frame([("H", x, 0.0, 0.0)])])

    monkeypatch.setattr(structure_io, "parse", parse)

    assert main._comparison_frame([log, trajectory]).atoms[0][1] == 2.0


def test_comparison_geometry_does_not_eagerly_read_packed_trajectories(tmp_path, monkeypatch):
    from oqp_studio import main, structure_io

    trajectory = tmp_path / "dynamics.namd.trj"
    trajectory.write_bytes(b"packed trajectory")

    def parse(name, data, path=None):
        assert name == trajectory.name
        assert data == b""
        assert path == str(trajectory)
        return structure_io.Structure("namd.trj", [
            structure_io.Frame([("H", 0.0, 0.0, 0.0)]),
        ])

    monkeypatch.setattr(structure_io, "parse", parse)

    assert main._comparison_frame([trajectory]).atoms[0][0] == "H"


def test_comparison_geometry_prefers_an_explicit_result_xyz(tmp_path):
    import json

    from oqp_studio import main

    input_json = tmp_path / "input.json"
    input_json.write_text(json.dumps({
        "elements": ["H"], "coordinates": [0.0, 0.0, 0.0], "unit": "angstrom",
    }))
    result_xyz = tmp_path / "result.xyz"
    result_xyz.write_text("1\nfinal\nH 1.0 2.0 3.0\n")

    assert main._comparison_frame([input_json, result_xyz]).atoms[0][1:] == (1.0, 2.0, 3.0)


def test_comparison_geometry_ranks_generic_text_after_molden(tmp_path, monkeypatch):
    from oqp_studio import main, structure_io

    table = tmp_path / "charges.txt"
    molden = tmp_path / "result.molden"
    table.write_text("C 99.0 0.0 0.0\n")
    molden.write_text("molden")

    def parse(name, *_args, **_kwargs):
        x = 1.0 if name.endswith(".molden") else 99.0
        return structure_io.Structure("test", [structure_io.Frame([("C", x, 0.0, 0.0)])])

    monkeypatch.setattr(structure_io, "parse", parse)

    assert main._comparison_frame([table, molden]).atoms[0][1] == 1.0


def test_comparison_geometry_promotes_a_verified_openqp_text_log(tmp_path):
    from oqp_studio import main

    input_file = tmp_path / "input.oqp"
    input_file.write_text('hf/sto-3g\nenergy\ngeom="""\nH 0.0 0.0 0.0\n"""\n')
    output_file = tmp_path / "results.txt"
    output_file.write_text(
        "Cartesian Coordinate in Angstrom\n"
        "--------------------------------\n"
        "1 1.0 2.0 0.0 0.0\n"
    )

    assert main._comparison_frame([input_file, output_file]).atoms[0][1] == 2.0


def test_comparison_geometry_does_not_read_text_when_result_xyz_exists(tmp_path, monkeypatch):
    from oqp_studio import main

    result_xyz = tmp_path / "result.xyz"
    result_xyz.write_text("1\nfinal\nH 1.0 2.0 3.0\n")
    large_text = tmp_path / "results.txt"
    large_text.write_text("must not be read")
    original_read_bytes = type(large_text).read_bytes

    def guarded_read(path):
        if path == large_text:
            raise AssertionError("lower-priority text was read eagerly")
        return original_read_bytes(path)

    monkeypatch.setattr(type(large_text), "read_bytes", guarded_read)
    assert main._comparison_frame([large_text, result_xyz]).atoms[0][1] == 1.0


@pytest.mark.parametrize("suffix", [".log", ".out", ".txt"])
def test_comparison_geometry_streams_openqp_logs(tmp_path, monkeypatch, suffix):
    from oqp_studio import main

    output = tmp_path / f"result{suffix}"
    output.write_text(
        "Cartesian Coordinate in Angstrom\n"
        "--------------------------------\n"
        "1 1.0 4.0 0.0 0.0\n"
    )
    original_read_bytes = type(output).read_bytes

    def guarded_read(path):
        if path == output:
            raise AssertionError("OpenQP log was materialized instead of streamed")
        return original_read_bytes(path)

    monkeypatch.setattr(type(output), "read_bytes", guarded_read)
    assert main._comparison_frame([output]).atoms[0][1] == 4.0


def test_comparison_geometry_uses_xyz_content_from_a_text_file(tmp_path):
    from oqp_studio import main

    structure = tmp_path / "geometry.txt"
    structure.write_text("1\ngeometry\nH 3.0 2.0 1.0\n")

    assert main._comparison_frame([structure]).atoms[0][1:] == (3.0, 2.0, 1.0)


def test_comparison_geometry_discards_oversized_log_blocks(tmp_path, monkeypatch):
    from oqp_studio import main

    monkeypatch.setattr(main, "MAX_COMPARISON_ATOMS", 2)
    output = tmp_path / "result.log"
    output.write_text(
        "Cartesian Coordinate in Angstrom\n-----\n"
        "1 1.0 0.0 0.0 0.0\n2 1.0 1.0 0.0 0.0\n3 1.0 2.0 0.0 0.0\n\n"
        "Cartesian Coordinate in Angstrom\n-----\n1 1.0 7.0 0.0 0.0\n"
    )

    frame = main._comparison_frame([output])
    assert len(frame.atoms) == 1
    assert frame.atoms[0][1] == 7.0


def test_comparison_geometry_discards_overlong_log_lines(tmp_path, monkeypatch):
    from oqp_studio import main

    monkeypatch.setattr(main, "MAX_COMPARISON_LINE", 64)
    output = tmp_path / "result.log"
    output.write_text(
        "x" * 1000 + "\n"
        "Cartesian Coordinate in Angstrom\n-----\n1 1.0 5.0 0.0 0.0\n"
    )

    assert main._comparison_frame([output]).atoms[0][1] == 5.0


def test_comparison_geometry_uses_unknown_result_suffix_as_a_last_resort(tmp_path):
    from oqp_studio import main

    structure = tmp_path / "geometry.dat"
    structure.write_text("1\ngeometry\nH 0.0 0.0 0.0\n")

    assert main._comparison_frame([structure]).atoms[0][0] == "H"


def test_comparison_geometry_detects_an_openqp_log_with_an_unknown_suffix(tmp_path):
    from oqp_studio import main

    output = tmp_path / "result.dat"
    output.write_text(
        "Cartesian Coordinate in Angstrom\n"
        "--------------------------------\n"
        "1 8.0 4.0 0.0 0.0\n"
    )

    frame = main._comparison_frame([output])
    assert frame.atoms == [("O", 4.0, 0.0, 0.0)]


def test_project_comparison_rejects_an_active_calculation(monkeypatch):
    from datetime import datetime, timezone

    from oqp_studio import jobs, main

    class ActiveManager:
        def get(self, job_id):
            return jobs.JobInfo(
                id=job_id, name=job_id,
                status=jobs.JobStatus.running if job_id == "active" else jobs.JobStatus.done,
                runner="bundled", created_at=datetime.now(timezone.utc).isoformat(),
            )

    monkeypatch.setattr(main, "manager", ActiveManager())

    response = client.get("/api/comparison?left=done&right=active")

    assert response.status_code == 409
    assert "completed projects" in response.json()["detail"]


def test_dyson_strength_is_reported_as_twice_its_occupation(monkeypatch):
    """Molden's Occup field is an OpenQP Dyson strength, not an occupancy."""
    import numpy as np

    from oqp_studio import main, molden

    data = molden.MoldenData(
        atoms=[],
        orbitals=[
            molden.Orbital(1, -0.5, "Alpha", 1.0004,
                           "Dyson-IP-state-2", np.zeros(0)),
        ],
    )
    monkeypatch.setattr(main, "_load_molden", lambda _job, _name: data)

    orbital = main.molden_orbitals("job", "dyson.molden")["orbitals"][0]
    assert orbital["kind"] == "dyson"
    assert orbital["dyson_kind"] == "IP"
    assert orbital["state_index"] == 2
    assert orbital["strength"] == 1.0004
    assert orbital["occupation"] == 2.0008


def test_normal_modes_and_trajectory():
    from pathlib import Path

    from oqp_studio import molden

    sample = (
        Path(__file__).resolve().parents[2]
        / "frontend" / "public" / "viewer" / "samples" / "water-hessian-mo.freq.molden"
    )
    vib = molden.parse_vibrations(sample.read_text())
    assert len(vib.atoms) == 3
    assert len(vib.modes) == 3  # water: 3N-6
    assert vib.modes[0].frequency > 0
    assert vib.modes[0].displacements.shape == (3, 3)

    traj = molden.mode_trajectory(vib, 1, frames=8)
    lines = traj.splitlines()
    assert len(lines) == 8 * (3 + 2)  # frames x (count + comment + atoms)
    assert lines[0].strip() == "3"
    assert "cm-1" in lines[1]


def test_runner_hands_the_engine_an_absolute_input_path(tmp_path):
    """The engine runs with cwd=job_dir, so a relative path would not resolve.

    It also decides where results land: an engine writes its .log and .json
    next to the input it was given.
    """
    from oqp_studio.runners.base import Runner

    class Recording(Runner):
        name = "recording"
        seen: ClassVar[list[str]] = []

        def is_available(self) -> bool:
            return True

        def build_command(self, input_file):
            Recording.seen.append(str(input_file))
            return ["true"]

    job_dir = tmp_path / "jobs_data" / "abc123"
    job_dir.mkdir(parents=True)
    (job_dir / "input.oqp").write_text("system=water\n")

    assert Recording().run(job_dir) == 0
    passed = Path(Recording.seen[-1])
    assert passed.is_absolute()
    assert passed == (job_dir / "input.oqp").resolve()


def test_imported_results_analyse_like_a_run_of_our_own(tmp_path, monkeypatch):
    """A run made elsewhere must be analysable, not just downloadable.

    Analysis is keyed on a job directory, so importing means giving foreign
    output files one; the summary and spectrum endpoints then work unchanged.
    """
    import json

    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    export = json.dumps({
        "energy": -76.0107465151,
        "frequencies_cm-1": [1600.0, 3700.0, 3800.0],
        "infrared_intensities": [70.0, 10.0, 60.0],
    })
    res = client.post(
        "/api/jobs/import",
        files=[
            ("files", ("hf.log", "TOTAL energy =  -76.0107465151\n", "text/plain")),
            ("files", ("hf.json", export, "application/json")),
        ],
        data={"name": "cluster run"},
    )
    assert res.status_code == 200, res.text
    job_id = res.json()["id"]
    assert res.json()["runner"] == "imported"

    summary = client.get(f"/api/jobs/{job_id}/summary").json()
    assert summary["energy"]["total"] == -76.0107465151
    assert summary["has_frequencies"]
    assert client.get(f"/api/jobs/{job_id}/spectrum?kind=ir").json()["x"]


def test_ekt_dyson_roots_produce_a_photoelectron_spectrum(tmp_path):
    """IP roots are electron-removal energies, not optical excited states."""
    from oqp_studio import analysis

    log = tmp_path / "water.log"
    log.write_text("""
MRSF-EKT ionization potentials (Dyson roots)
Rows index EKT Dyson roots/orbitals, not TDDFT excited states.
dyson      eig(ha)        eBE(ha)        eBE(eV)       metric       strength
    1       -0.591801        0.591801       16.1037      1.000000      1.006008
    2       -0.453861        0.453861       12.3502      1.000000      1.003115
""")

    summary = analysis.summarize([log])
    assert not summary["has_states"]
    assert summary["has_ekt_ip"]
    assert summary["ekt"]["ip"][0]["binding_ev"] == 16.1037
    data = analysis.spectrum(summary, "photoelectron", fwhm=0.2)
    assert data["available"]
    assert data["x_label"] == "Electron binding energy (eV)"
    assert max(data["y"]) == 1.0


def test_mrsf_state_table_produces_an_optical_absorption_spectrum(tmp_path):
    """MRSF IP/EA runs print an optical state table before their EKT roots."""
    from oqp_studio import analysis

    log = tmp_path / "water.log"
    log.write_text("""
State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
           Hartree           eV        rel. S0                  X          Y          Z        Abs.      strength
    S0    -76.3600325966    -7.428327     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
   REF    -76.0870465998     0.000000     7.428327        (triplet ROHF/UHF internal working reference)
    S1    -76.0412293756     1.246750     8.675078      0.000     0.0000    -0.0000     0.2401     0.2401      0.0123
    S2    -75.9853724356     2.766695    10.195022      0.000    -0.0000     0.0000     0.0001     0.0001      0.0000
""")

    summary = analysis.summarize([log])
    assert summary["has_states"]
    assert summary["states"][1]["excitation_ev"] == 8.675078
    assert summary["states"][1]["oscillator"] == 0.0123
    data = analysis.spectrum(summary, "absorption")
    assert data["available"]
    assert data["title"] == "Absorption from S0"


def test_mrsf_transition_table_uses_the_state_specific_esa_strength(tmp_path):
    from oqp_studio import analysis

    log = tmp_path / "water.log"
    log.write_text("""
State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
    S0    -76.3600325966    -7.428327     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
    S1    -76.0412293756     1.246750     8.675078      0.000     0.0000    -0.0000     0.2401     0.2401      0.0123
    S2    -75.9853724356     2.766695    10.195022      0.000    -0.0000     0.0000     0.0001     0.0001      0.0000

Transition   Excitation         Transition dipole, a.u.                   Oscillator
                 eV              x          y          z         Abs.       strength
   S0 -> S1     8.675078        0.0000    -0.0000     0.2401      0.2401       0.0123
   S0 -> S2    10.195022       -0.0000     0.0000     0.0001      0.0001       0.0000
   S1 -> S2     1.519945        1.7982    -0.0046    -0.0000      1.7982       0.1204
""")

    summary = analysis.summarize([log])
    assert summary["transitions"][-1]["oscillator"] == 0.1204
    data = analysis.spectrum(summary, "esa", state=1)
    assert data["available"]
    assert data["sticks"][0]["position"] == analysis.NM_EV / 1.519945
    assert max(data["x"]) <= 2 * analysis.NM_EV / 1.519945 + 1e-6


def test_optimization_history_keeps_each_geometrys_electronic_spectrum(tmp_path):
    """An optimization trajectory must not borrow the final spectrum for every step."""
    from oqp_studio import analysis

    trajectory = tmp_path / "opt_geom.xyz"
    trajectory.write_text("""3
Geom 1 -76.0412293829
O 0.000 0.000 0.000
H 0.000 0.757 -0.469
H 0.000 -0.757 -0.469
3
Geom 2 -76.0685554095
O 0.010 0.000 0.000
H 0.000 0.757 -0.469
H 0.000 -0.757 -0.469
""")
    status = tmp_path / "opt_status.txt"
    status.write_text("""Step Energy Shift RMSD Step Max Step RMSD Grad Max Grad
1 -76.04122938 -76.04122938 0.000000 0.000000 0.061977 0.119177
2 -76.06855541 -0.02732603 0.100467 0.200173 0.013577 0.022500
""")
    log = tmp_path / "water.log"
    log.write_text("""Geometry Optimization Step 1
State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
    S0    -76.3600325966    -7.428327     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
    S1    -76.0412293756     1.246750     8.675078      0.000     0.0000    -0.0000     0.2401     0.2401      0.0123
    S2    -75.9853724356     2.766695    10.195022      0.000    -0.0000     0.0000     0.0001     0.0001      0.0000
Transition   Excitation         Transition dipole, a.u.                   Oscillator
   S1 -> S2     1.519945        1.7982    -0.0046    -0.0000     1.7982       0.1204
Geometry Optimization Step 2
State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
    S0    -76.3285260128    -6.828327     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
    S1    -76.0685554095     0.246750     7.074160      0.000     0.0000    -0.0000     0.2401     0.2401      0.0048
    S2    -76.0227501977     1.766695     8.320584      0.000    -0.0000     0.0000     0.0001     0.0001      0.0000
Transition   Excitation         Transition dipole, a.u.                   Oscillator
   S1 -> S2     1.246423        1.7982    -0.0046    -0.0000     1.7982       0.0915
""")

    history = analysis.optimization_history([trajectory, status, log])
    first, second = history["steps"]
    assert len(history["steps"]) == 2
    assert first["energy"] == -76.04122938
    assert first["rmsd_grad"] == 0.061977
    assert first["states"][1]["excitation_ev"] == 8.675078
    assert second["states"][1]["excitation_ev"] == 7.07416
    assert second["transitions"][0]["oscillator"] == 0.0915


def test_emission_is_available_only_after_excited_state_optimization(tmp_path):
    from oqp_studio import analysis

    log = tmp_path / "water.log"
    log.write_text("""
State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
    S0    -76.3600325966    -7.428327     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
    S1    -76.0412293756     1.246750     8.675078      0.000     0.0000    -0.0000     0.2401     0.2401      0.0123
""")
    input_file = tmp_path / "water.oqp"
    input_file.write_text("mrsf/bhhlyp/6-31g\nopt(S1)\n")

    vertical = analysis.summarize([log])
    assert not analysis.spectrum(vertical, "emission")["available"]
    optimized = analysis.summarize([input_file, log])
    assert optimized["excited_state_optimized"] == 1
    assert analysis.spectrum(optimized, "emission")["available"]


def test_importing_a_folder_takes_the_results_and_leaves_the_rest(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path / "jobs")
    run = tmp_path / "run"
    run.mkdir()
    (run / "water.log").write_text("TOTAL energy =  -76.0\n")
    (run / "water.json").write_text("{}")
    (run / "restart.bin").write_bytes(b"\0" * 16)

    res = client.post("/api/jobs/import-path", json={"path": str(run)})
    assert res.status_code == 200, res.text
    names = [f["name"] for f in client.get(f"/api/jobs/{res.json()['id']}/files").json()]
    assert names == ["water.json", "water.log"]

    assert client.post("/api/jobs/import-path", json={"path": str(tmp_path / "no")}).status_code == 404


def test_an_import_cannot_write_outside_its_job_directory(tmp_path, monkeypatch):
    from oqp_studio import jobs

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path / "jobs")
    res = client.post(
        "/api/jobs/import",
        files=[("files", ("../../escaped.log", "x", "text/plain"))],
    )
    assert res.status_code == 200, res.text
    assert not (tmp_path / "escaped.log").exists()
    names = [f["name"] for f in client.get(f"/api/jobs/{res.json()['id']}/files").json()]
    assert names == ["escaped.log"]


def test_an_installer_that_carries_the_engine_needs_no_download(tmp_path, monkeypatch):
    """An all-in-one installer must just work: no download, no configuration.

    The shell tells the backend where its resources went, because that
    directory is in a different place on every platform.
    """
    import os
    import stat

    from oqp_studio import engine, environment

    resources = tmp_path / "resources"
    (resources / "engine").mkdir(parents=True)
    executable = resources / "engine" / engine.EXECUTABLE
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(executable.stat().st_mode | stat.S_IEXEC)

    monkeypatch.setenv("OQP_STUDIO_RESOURCES", str(resources))
    monkeypatch.setattr(engine, "_bundled", None)
    monkeypatch.setattr(engine, "_bundled_resolved", False)
    monkeypatch.setattr(engine, "data_dir", lambda: tmp_path / "unused")
    # No engine of the user's own on PATH, so the bundled one is what is found.
    monkeypatch.setattr(environment, "locate", lambda command: None)

    assert engine.bundled_dir() == resources / "engine"
    assert engine.locate() == str(executable)
    assert engine.status()["source"] == "included with the installer"
    assert os.access(executable, os.X_OK)


def test_updating_an_all_in_one_install_keeps_the_engine(monkeypatch):
    """The slim installer would replace the directory the engine lives in.

    Someone who installed the all-in-one must be offered the all-in-one
    again, or the update would silently leave them unable to compute.
    """
    from oqp_studio import engine, update

    monkeypatch.setattr(update, "_platform_asset", lambda: ("windows-x64", "-setup.exe"))
    assets = [
        {"name": "OQP-Studio-0.2.0-windows-x64-setup.exe"},
        {"name": "OQP-Studio-0.2.0-windows-x64-with-engine-setup.exe"},
    ]

    monkeypatch.setattr(engine, "bundled_dir", lambda: None)
    assert update.pick_asset(assets)["name"].endswith("windows-x64-setup.exe")

    monkeypatch.setattr(engine, "bundled_dir", lambda: object())
    assert update.pick_asset(assets)["name"].endswith("with-engine-setup.exe")

    # A release that published no all-in-one still updates the app.
    assert update.pick_asset(assets[:1])["name"].endswith("windows-x64-setup.exe")


def test_results_land_in_documents_where_the_user_can_find_them(tmp_path, monkeypatch):
    """Two things at once, both of which broke a real install.

    The working directory is never it: an app started from Finder has "/",
    where the old relative default resolved to /jobs_data and macOS refused
    the mkdir, so every run died at submission with a bare 500. And of the
    directories the user does own, results belong in Documents -- Finder
    hides ~/Library, so anything written there is invisible.
    """
    from pathlib import Path

    from oqp_studio import workspace

    home = tmp_path / "home"
    (home / "Documents").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("OQP_STUDIO_CONFIG", str(tmp_path / "cfg" / "network.json"))
    monkeypatch.delenv("OQP_STUDIO_JOBS", raising=False)
    monkeypatch.chdir("/")

    root = workspace.resolve()
    assert root.is_absolute()
    assert root == home / "Documents" / "OQP Studio" / "jobs"
    assert root.is_dir()


def test_a_refused_documents_folder_does_not_stop_the_app(tmp_path, monkeypatch):
    """macOS asks before an app may write to Documents, and can be refused."""
    from pathlib import Path

    from oqp_studio import workspace

    home = tmp_path / "home"
    documents = home / "Documents"
    documents.mkdir(parents=True)
    # Something in the way that no uid can mkdir through -- a refusal that
    # reproduces whether or not the test runs as root.
    (documents / "OQP Studio").write_text("not a directory")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("OQP_STUDIO_CONFIG", str(tmp_path / "cfg" / "network.json"))
    monkeypatch.delenv("OQP_STUDIO_JOBS", raising=False)

    root = workspace.resolve()
    assert root.is_dir()
    assert documents not in root.parents


def test_the_user_chooses_where_results_are_written(tmp_path, monkeypatch):
    """Results are the user's data, so the folder is theirs to pick."""
    from oqp_studio import jobs, main, workspace

    monkeypatch.setenv("OQP_STUDIO_CONFIG", str(tmp_path / "cfg" / "network.json"))
    monkeypatch.delenv("OQP_STUDIO_JOBS", raising=False)
    # This API-level test must not inherit a real job still running in the
    # developer's workspace from an earlier test or local server session.
    isolated_manager = jobs.JobManager()
    monkeypatch.setattr(jobs, "manager", isolated_manager)
    monkeypatch.setattr(main, "manager", isolated_manager)
    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path / "initial-jobs")
    chosen = tmp_path / "My Calculations"

    res = client.post("/api/workspace", json={"jobs_dir": str(chosen)})
    assert res.status_code == 200, res.text
    assert res.json()["active"] == str(chosen.resolve())
    assert jobs.JOBS_ROOT == chosen.resolve()
    assert chosen.is_dir()
    # It has to survive a restart, so it is on disk, not just in memory.
    assert workspace.configured() == str(chosen)

    # And the job that follows lands there.
    job = client.post("/api/jobs", json={"input_text": "energy\n", "runner": "local"})
    assert job.status_code == 200, job.text
    assert (chosen / job.json()["id"]).is_dir()
    # Submission is asynchronous; this test only needs to prove placement.
    isolated_manager.get(job.json()["id"]).status = jobs.JobStatus.done

    # A directory that cannot be written is refused with the reason, not a 500.
    blocked = tmp_path / "blocked"
    blocked.write_text("a file, not a directory")
    bad = client.post("/api/workspace", json={"jobs_dir": str(blocked)})
    assert bad.status_code == 400
    assert str(blocked) in bad.json()["detail"]

    # Clearing it hands the choice back to the app.
    back = client.post("/api/workspace", json={"jobs_dir": ""})
    assert back.status_code == 200
    assert workspace.configured() == ""


def test_starting_the_server_touches_no_directory_it_may_be_denied(tmp_path, monkeypatch):
    """Importing must not create anything, or the app never finishes starting.

    macOS asks the user before an app may write to Documents and blocks until
    they answer. Doing that while the module is imported happens before the
    server binds its port, and the shell -- which waits thirty seconds --
    reports the backend as failing to start.
    """
    from pathlib import Path

    from oqp_studio import jobs, workspace

    home = tmp_path / "home"
    (home / "Documents").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("OQP_STUDIO_CONFIG", str(tmp_path / "cfg" / "network.json"))
    monkeypatch.delenv("OQP_STUDIO_JOBS", raising=False)

    wanted = workspace.preferred()
    assert wanted == home / "Documents" / "OQP Studio" / "jobs"
    assert not wanted.exists(), "preferred() created a directory"

    # A fresh manager is what import builds; it must not reach the disk either.
    manager = jobs.JobManager()
    assert not wanted.exists(), "constructing the job manager created a directory"

    # The first job is what creates it.
    monkeypatch.setattr(jobs, "JOBS_ROOT", wanted)
    assert manager._ensure() == wanted.resolve()
    assert wanted.is_dir()


def _cube(values: list[float], *, origin: float = 0.0) -> str:
    return (
        "test cube\nvalues\n"
        f"0 {origin:.6f} 0.0 0.0\n"
        "2 1.0 0.0 0.0\n"
        "1 0.0 1.0 0.0\n"
        "1 0.0 0.0 1.0\n"
        + " ".join(str(value) for value in values) + "\n"
    )


def test_cube_arithmetic_preserves_the_grid_and_combines_values():
    from oqp_studio import cube

    result = cube.parse(cube.combine(_cube([1.5, -2.0]), _cube([0.5, 3.0]), "difference"))

    assert result.shape == (2, 1, 1)
    assert result.values == [1.0, -5.0]


def test_cube_geometry_is_extracted_in_angstrom():
    from io import StringIO

    import pytest

    from oqp_studio import cube

    text = (
        "test cube\nvalues\n"
        "1 0.0 0.0 0.0\n"
        "1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n"
        "8 0.0 1.0 0.0 0.0\n0.0\n"
    )
    xyz = cube.geometry_xyz(cube.parse(text))

    assert xyz.startswith("1\ncube geometry\nO ")
    assert "0.5291772109 0.0000000000 0.0000000000" in xyz
    assert cube.geometry_xyz(cube.parse_header(StringIO(text[:-4] + "not-grid-data\n"))) == xyz
    with pytest.raises(ValueError, match="header lines are limited"):
        cube.parse_header(StringIO("x" * (cube.MAX_HEADER_LINE + 1) + "\n" + text))
    orbital_header = text.replace("1 0.0 0.0 0.0", "-1 0.0 0.0 0.0 2000000", 1)
    assert cube.geometry_xyz(cube.parse_header(StringIO(orbital_header))) == xyz
    fractional = text.replace("8 0.0 1.0 0.0 0.0", "1.6 0.0 1.0 0.0 0.0")
    with pytest.raises(ValueError, match="nonintegral atomic number"):
        cube.geometry_xyz(cube.parse(fractional))


def test_cube_header_has_an_aggregate_size_limit():
    from io import StringIO

    import pytest

    from oqp_studio import cube

    atom_line = "1" + " " * (cube.MAX_HEADER_LINE - 16) + "0 0 0 0\n"
    text = (
        "cube\nvalues\n"
        "300 0 0 0\n1 1 0 0\n1 0 1 0\n1 0 0 1\n"
        + atom_line * 300
    )

    with pytest.raises(ValueError, match="header is limited to 4 MiB"):
        cube.parse_header(StringIO(text))
    with pytest.raises(ValueError, match="header is limited to 4 MiB"):
        cube.parse(text)

    padded_identifiers = (
        "cube\nvalues\n"
        "-1 0 0 0\n1 1 0 0\n1 0 1 0\n1 0 0 1\n"
        "1 0 0 0 0\n"
        "1" + " " * cube.MAX_CUBE_HEADER + "7\n0.0\n"
    )
    with pytest.raises(ValueError, match="header is limited to 4 MiB"):
        cube.parse(padded_identifiers)


def test_cube_arithmetic_rejects_different_grids():
    import pytest

    from oqp_studio import cube

    with pytest.raises(ValueError, match="different origins"):
        cube.combine(_cube([1.0, 2.0]), _cube([1.0, 2.0], origin=0.1), "sum")


def test_cube_arithmetic_supports_multiline_datasets_and_fortran_headers():
    from oqp_studio import cube

    def multi(values):
        return (
            "orbital cube\nvalues\n"
            "-1 0.0D+00 0.0 0.0 2\n"
            "1 1.0D+00 0.0 0.0\n"
            "1 0.0 1.0D+00 0.0\n"
            "1 0.0 0.0 1.0D+00\n"
            "1 0.0D+00 0.0 0.0 0.0\n"
            "2 10\n20\n"
            + " ".join(values) + "\n"
        )

    result = cube.parse(cube.combine(
        multi(["1.0D+00", "2.0D+00"]),
        multi(["0.5D+00", "4.0D+00"]),
        "difference",
    ))

    assert result.datasets == 2
    assert result.dataset_ids == (10, 20)
    assert result.values == [0.5, -2.0]

    wrapped_differently = multi(["0.5D+00", "4.0D+00"]).replace("2 10\n20\n", "2 10 20\n")
    result = cube.parse(cube.combine(
        multi(["1.0D+00", "2.0D+00"]), wrapped_differently, "difference",
    ))
    assert result.values == [0.5, -2.0]


def test_cube_parser_rejects_non_finite_header_and_grid_values():
    import pytest

    from oqp_studio import cube

    with pytest.raises(ValueError, match="non-finite"):
        cube.parse(_cube([float("nan"), 1.0]))
    with pytest.raises(ValueError, match="non-finite"):
        cube.combine(
            _cube([1.0, 2.0]).replace("0.000000 0.0 0.0", "Inf 0.0 0.0", 1),
            _cube([1.0, 2.0]),
            "sum",
        )
    with pytest.raises(ValueError, match="arithmetic produced a non-finite"):
        cube.combine(_cube([1.0e308, 1.0]), _cube([1.0e308, 1.0]), "sum")


def test_cube_parser_stops_at_the_first_grid_value_beyond_the_declared_shape(monkeypatch):
    import pytest

    from oqp_studio import cube

    converted = 0
    original = cube._number

    def count_number(token):
        nonlocal converted
        converted += 1
        return original(token)

    monkeypatch.setattr(cube, "_number", count_number)
    with pytest.raises(ValueError, match="more than the expected 2"):
        cube.parse(_cube([1.0, 2.0, *([0.0] * 10_000)]))

    assert converted == 2


def test_cube_parser_does_not_materialize_every_grid_line():
    import pytest

    from oqp_studio import cube

    class NoSplit(str):
        def splitlines(self, *_args, **_kwargs):
            raise AssertionError("the complete cube was split into a line list")

    with pytest.raises(ValueError, match="more than the expected 2"):
        cube.parse(NoSplit(_cube([1.0, 2.0]) + "0\n" * 10_000))


def test_symmetry_handles_many_leading_blank_rows():
    from oqp_studio import symmetry

    assert symmetry.analyze("\n" * 50_000 + "He 0 0 0\n")["point_group"] == "Kh"


def test_cube_parser_rejects_a_negative_dataset_count():
    import pytest

    from oqp_studio import cube

    malformed = _cube([1.0, 2.0]).replace(
        "0 0.000000 0.0 0.0", "0 0.000000 0.0 0.0 -2", 1,
    )

    with pytest.raises(ValueError, match="dataset count must be positive"):
        cube.parse(malformed)

    fractional = _cube([1.0, 2.0]).replace(
        "0 0.000000 0.0 0.0", "1 0.000000 0.0 0.0", 1,
    ).replace("1 0.0 0.0 1.0\n", "1 0.0 0.0 1.0\n1.6 0.0 0.0 0.0 0.0\n", 1)
    with pytest.raises(ValueError, match="nonintegral atomic number"):
        cube.combine(fractional, fractional, "sum")


def test_cube_parser_bounds_dataset_identifiers_before_accumulating_them():
    import pytest

    from oqp_studio import cube

    malformed = (
        "orbital cube\nvalues\n"
        "-1 0.0 0.0 0.0\n"
        "1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n"
        "1 0.0 0.0 0.0 0.0\n"
        "999999999 " + "1 " * 10_000 + "\n0.0\n"
    )

    with pytest.raises(ValueError, match="identifier count exceeds"):
        cube.parse(malformed)


def test_cube_parser_stops_fixed_header_records_after_the_first_extra_field(monkeypatch):
    import pytest

    from oqp_studio import cube

    original = cube._tokens
    seen = 0

    def count_tokens(line):
        nonlocal seen
        for token in original(line):
            seen += 1
            if seen > 6:
                raise AssertionError("fixed cube record was consumed without a bound")
            yield token

    monkeypatch.setattr(cube, "_tokens", count_tokens)
    malformed = _cube([1.0, 2.0]).replace(
        "0 0.000000 0.0 0.0", " ".join(["0"] * 10_000), 1,
    )

    with pytest.raises(ValueError, match="cube header is invalid"):
        cube.parse(malformed)
    assert seen == 6


def test_cube_parser_rejects_malformed_geometry_records_and_empty_axes():
    import pytest

    from oqp_studio import cube

    valid = _cube([1.0, 2.0])
    malformed_atom = valid.replace(
        "0 0.000000 0.0 0.0", "1 0.000000 0.0 0.0", 1,
    ).replace(
        "1 0.0 0.0 1.0\n", "1 0.0 0.0 1.0\n8 0.0 0.0 0.0\n", 1,
    )
    malformed = [
        valid.replace("0 0.000000 0.0 0.0", "0 0.000000 0.0"),
        valid.replace("2 1.0 0.0 0.0", "2 1.0 0.0"),
        valid.replace("2 1.0 0.0 0.0", "0 1.0 0.0 0.0"),
        valid.replace("0 0.000000 0.0 0.0", "0 0.000000 0.0 0.0 1 2"),
        malformed_atom,
    ]

    for text in malformed:
        with pytest.raises(ValueError, match="header is invalid"):
            cube.parse(text)


def test_cube_arithmetic_rejects_excessive_grids_and_mixed_axis_units():
    import pytest

    from oqp_studio import cube

    oversized = _cube([1.0, 2.0]).replace("2 1.0 0.0 0.0", "2000001 1.0 0.0 0.0")
    with pytest.raises(ValueError, match="limited to 2,000,000"):
        cube.parse(oversized)
    angstrom_cube = (
        _cube([1.0, 2.0])
        .replace("2 1.0 0.0 0.0", "-2 1.0 0.0 0.0")
        .replace("1 0.0 1.0 0.0", "-1 0.0 1.0 0.0")
        .replace("1 0.0 0.0 1.0", "-1 0.0 0.0 1.0")
    )
    with pytest.raises(ValueError, match="different coordinate units"):
        cube.combine(_cube([1.0, 2.0]), angstrom_cube, "difference")
    with pytest.raises(ValueError, match="header is invalid"):
        cube.parse(_cube([1.0, 2.0]).replace("2 1.0 0.0 0.0", "-2 1.0 0.0 0.0"))


def test_cube_arithmetic_endpoint_uses_only_job_files(tmp_path, monkeypatch):
    from oqp_studio import cube, jobs, main

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    job = manager.adopt("cube pair", [
        ("left.cube", _cube([1.0, 2.0]).encode()),
        ("right.cube", _cube([0.25, 0.5]).encode()),
    ])
    monkeypatch.setattr(main, "manager", manager)

    response = client.get(
        f"/api/jobs/{job.id}/cube-combine",
        params={"left": "left.cube", "right": "right.cube", "operation": "sum"},
    )
    escaped = client.get(
        f"/api/jobs/{job.id}/cube-combine",
        params={"left": "../left.cube", "right": "right.cube"},
    )

    assert response.status_code == 200
    assert cube.parse(response.text).values == [1.25, 2.5]
    assert escaped.status_code == 404


def test_cube_arithmetic_endpoint_caps_combined_input_size(tmp_path, monkeypatch):
    from oqp_studio import jobs, main

    monkeypatch.setattr(jobs, "JOBS_ROOT", tmp_path)
    manager = jobs.JobManager()
    manager._ready = True
    job = manager.adopt("large cubes", [
        ("left.cube", _cube([1.0, 2.0]).encode()),
        ("right.cube", _cube([1.0, 2.0]).encode()),
    ])
    for name in ("left.cube", "right.cube"):
        (tmp_path / job.id / name).open("ab").truncate(33 * 1024 * 1024)
    monkeypatch.setattr(main, "manager", manager)

    response = client.get(
        f"/api/jobs/{job.id}/cube-combine",
        params={"left": "left.cube", "right": "right.cube"},
    )

    assert response.status_code == 413
    assert "64 MiB combined" in response.json()["detail"]


def test_symmetry_identifies_water_and_equivalent_hydrogens():
    from oqp_studio import symmetry

    result = symmetry.analyze(
        "O 0.000000 0.000000 0.117300\n"
        "H 0.000000 0.757200 -0.469200\n"
        "H 0.000000 -0.757200 -0.469200\n",
        tolerance=0.01,
    )

    assert result["point_group"] == "C2v"
    assert [2, 3] in result["equivalent_atoms"]
    assert result["max_deviation_angstrom"] < 1.0e-8


def test_symmetry_identifies_linear_centrosymmetric_co2():
    from oqp_studio import symmetry

    result = symmetry.analyze(
        "O 0.0 0.0 -1.16\nC 0.0 0.0 0.0\nO 0.0 0.0 1.16\n",
        tolerance=0.01,
    )

    assert result["point_group"] == "Dinfh"
    assert [1, 3] in result["equivalent_atoms"]


def test_symmetry_uses_maximum_per_atom_deviation_for_linearity():
    from oqp_studio import symmetry

    result = symmetry.analyze(
        "C -2.0 0.04 0.0\nC -1.0 -0.04 0.0\n"
        "C 1.0 -0.04 0.0\nC 2.0 0.04 0.0\n",
        tolerance=0.05,
    )

    assert result["point_group"] in {"Cinfv", "Dinfh"}


def test_symmetry_classifies_an_isolated_atom_as_spherical():
    from oqp_studio import symmetry

    assert symmetry.analyze("He 2.0 -3.0 4.0\n")["point_group"] == "Kh"


def test_symmetry_uses_heavy_element_masses_for_principal_axes():
    import pytest

    from oqp_studio import symmetry

    result = symmetry.analyze("H 0.0 0.0 0.0\nCs 1.0 0.0 0.0\n")

    assert result["center_angstrom"][0] == pytest.approx(132.905 / (132.905 + 1.008))


def test_symmetry_rejects_incomplete_or_malformed_coordinate_rows():
    import pytest

    from oqp_studio import symmetry

    with pytest.raises(ValueError, match="every coordinate row"):
        symmetry.analyze("3\nwater\nO 0 0 0\nH 0 1 0\nnot-an-atom\n")
    with pytest.raises(ValueError, match="every coordinate row"):
        symmetry.analyze("O 0 0 0\nH invalid 1 0\n")
    with pytest.raises(ValueError, match="at most 300"):
        symmetry.analyze("301\ntoo many\n")
    with pytest.raises(ValueError, match="exact element symbol"):
        symmetry.analyze("Hx 0.0 0.0 0.0\n")


def test_symmetry_rejects_a_structure_containing_only_dummy_sites():
    import pytest

    from oqp_studio import symmetry

    with pytest.raises(ValueError, match="positive mass"):
        symmetry.analyze("X 0.0 0.0 0.0\nX 1.0 0.0 0.0\n")


def test_symmetry_rejects_rotations_indistinguishable_at_the_tolerance():
    import numpy as np

    from oqp_studio import symmetry

    coordinates = np.asarray([[1.0, 0.0, 0.0]])
    almost_identity = symmetry._rotation(np.asarray([0.0, 0.0, 1.0]), 1.0e-6)

    assert not symmetry._moves_coordinates(coordinates, almost_identity, 0.01)
    high_order = symmetry._rotation(np.asarray([0.0, 0.0, 1.0]), 2 * np.pi / 100)
    assert not symmetry._moves_coordinates(np.asarray([[7.0, 0.0, 0.0]]), high_order, 0.5)


def test_symmetry_retains_atoms_near_an_axis_when_rotation_moves_them():
    import numpy as np

    from oqp_studio import symmetry

    coordinates = np.asarray([[0.0, -0.3, 0.0], [0.0, 0.3, 0.0]])
    orders = symmetry._rotation_orders(
        ["H", "H"], coordinates, np.asarray([1.0, 0.0, 0.0]), 0.5,
    )

    assert 2 in orders


def test_symmetry_clusters_indistinguishable_polyhedral_axes():
    import numpy as np

    from oqp_studio import symmetry

    operation = symmetry.Operation("C5", np.eye(3), [0], 0.0)
    axes = [
        np.asarray([0.0, 0.0, 1.0]),
        np.asarray([0.001, 0.0, 0.9999995]),
        np.asarray([-0.001, 0.0, 0.9999995]),
    ]
    distinct = symmetry._distinct_rotation_axes(
        [(axis, operation) for axis in axes], np.asarray([[1.0, 0.0, 0.0]]), 0.01,
    )

    assert len(distinct) == 1


def test_symmetry_preserves_exact_polyhedral_axes_when_tolerance_exceeds_radius():
    import numpy as np

    from oqp_studio import symmetry

    operation = symmetry.Operation("C3", np.eye(3), [0], 0.0)
    axes = [
        np.asarray(values, dtype=float) / np.sqrt(3.0)
        for values in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    ]
    distinct = symmetry._distinct_rotation_axes(
        [(axis, operation) for axis in axes], np.asarray([[0.4, 0.0, 0.0]]), 0.5,
    )

    assert len(distinct) == 4


def test_symmetry_rejects_oversized_json_before_endpoint_binding():
    from oqp_studio import main

    response = client.post(
        "/api/symmetry",
        content=b"x" * (main.MAX_SYMMETRY_REQUEST_BYTES + 1),
        headers={"Content-Type": "application/json"},
    )

    assert response.status_code == 413
    assert response.json()["detail"] == "symmetry request body is too large"
    assert main.MAX_SYMMETRY_REQUEST_BYTES < len(response.request.content)


def test_symmetry_assignment_handles_dense_hall_deficiency_without_backtracking():
    import numpy as np
    import pytest

    from oqp_studio import symmetry

    distances = np.zeros((20, 20))
    distances[:, -1] = 2.0

    assert symmetry._assignment(distances, 1.0) is None
    with pytest.raises(ValueError, match="work limit"):
        symmetry._assignment(np.zeros((101, 101)), 1.0)


def test_symmetry_stops_before_assignment_work_can_block_the_sidecar(monkeypatch):
    import pytest

    from oqp_studio import symmetry

    monkeypatch.setattr(symmetry, "MAX_MATCH_WORK", 0)
    monkeypatch.setattr(symmetry, "_match", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(symmetry, "_rotation_orders", lambda *_args, **_kwargs: list(range(2, 20)))
    xyz = "\n".join(
        f"C {index * 0.37:.6f} {(index ** 2) * 0.013:.6f} {(index ** 3) * 0.001:.6f}"
        for index in range(10)
    )

    with pytest.raises(ValueError, match="work limit"):
        symmetry.analyze(xyz, tolerance=0.01)


def test_symmetry_budgets_rotation_order_screening(monkeypatch):
    import pytest

    from oqp_studio import symmetry

    monkeypatch.setattr(symmetry, "MAX_ROTATION_SCREEN_WORK", 0)

    with pytest.raises(ValueError, match="rotation screening exceeded"):
        symmetry.analyze(
            "C 1.0 0.0 0.0\nC -1.0 0.0 0.0\nC 0.0 1.0 0.0\n",
            tolerance=0.001,
        )


def test_symmetry_finds_rotated_ammonia_mirror_planes():
    from math import cos, pi, sin

    from oqp_studio import symmetry

    angle = 17 * pi / 180
    rows = ["N 0.0 0.0 0.15"]
    for index in range(3):
        theta = angle + index * 2 * pi / 3
        rows.append(f"H {cos(theta):.12f} {sin(theta):.12f} -0.35")

    result = symmetry.analyze("\n".join(rows), tolerance=0.001)

    assert result["point_group"] == "C3v"
    assert result["operation_count"] == 6


def test_symmetry_uses_the_improper_axis_to_classify_d2d():
    from oqp_studio import symmetry

    allene = (
        "C 0.0 0.0 -1.3\nC 0.0 0.0 0.0\nC 0.0 0.0 1.3\n"
        "H 0.9 0.0 -1.8\nH -0.9 0.0 -1.8\n"
        "H 0.0 0.9 1.8\nH 0.0 -0.9 1.8\n"
    )

    assert symmetry.analyze(allene, tolerance=0.001)["point_group"] == "D2d"


def test_symmetry_closes_odd_improper_rotations_to_twice_the_order():
    from math import cos, pi, sin

    from oqp_studio import symmetry

    rows = ["B 0.0 0.0 0.0"] + [
        f"F {cos(index * 2 * pi / 3):.12f} {sin(index * 2 * pi / 3):.12f} 0.0"
        for index in range(3)
    ]
    result = symmetry.analyze("\n".join(rows), tolerance=0.001)

    assert result["point_group"] == "D3h"
    assert result["operation_count"] == 12


def test_symmetry_rejects_a_generator_when_a_required_power_fails():
    from math import cos, pi, sin

    import numpy as np

    from oqp_studio import symmetry

    coordinates = np.asarray([
        [cos(angle * pi / 180), sin(angle * pi / 180), 0.0]
        for angle in (0, 85, 180, 275)
    ])
    matrix = symmetry._rotation(np.asarray([0.0, 0.0, 1.0]), pi / 2)
    matched = symmetry._match(["C"] * 4, coordinates, matrix, 0.1)

    assert matched is not None
    assert symmetry._operation_powers("C4", matrix, matched, 4, coordinates, 0.1) is None


def test_symmetry_searches_rotation_orders_above_eight():
    from math import cos, pi, sin

    from oqp_studio import symmetry

    ring = "\n".join(
        f"C {cos(index * pi / 5):.12f} {sin(index * pi / 5):.12f} 0.0"
        for index in range(10)
    )

    assert symmetry.analyze(ring, tolerance=0.001)["point_group"] == "D10h"


def test_symmetry_detects_a_pure_improper_rotation_group():
    from math import cos, pi, sin

    from oqp_studio import symmetry

    rows = []
    for symbol, radius, height, offset in (("C", 1.0, 0.4, 0.0), ("H", 1.4, 0.7, 0.37)):
        for index in range(4):
            angle = offset + index * pi / 2
            z = height if index % 2 == 0 else -height
            rows.append(f"{symbol} {radius * cos(angle):.12f} {radius * sin(angle):.12f} {z:.12f}")

    assert symmetry.analyze("\n".join(rows), tolerance=0.001)["point_group"] == "S4"


def test_symmetry_screens_compact_improper_operations_by_their_full_displacement():
    from math import cos, pi, sin

    import numpy as np

    from oqp_studio import symmetry

    coordinates = []
    symbols = []
    for symbol, radius, height, offset in (("C", 0.30, 0.40, 0.0), ("H", 0.34, 0.70, 0.37)):
        for index in range(4):
            angle = offset + index * pi / 2
            z = height if index % 2 == 0 else -height
            symbols.append(symbol)
            coordinates.append([radius * cos(angle), radius * sin(angle), z])

    xyz = np.asarray(coordinates)
    axis = np.asarray([0.0, 0.0, 1.0])

    assert 4 not in symmetry._rotation_orders(symbols, xyz, axis, 0.5)
    assert 4 in symmetry._improper_orders(symbols, xyz, axis, 0.5)


def test_symmetry_improper_screening_allows_mixed_cycle_lengths():
    import numpy as np

    from oqp_studio import symmetry

    coordinates = np.asarray([
        [1.0, 0.0, 0.5], [0.0, 1.0, -0.5], [-1.0, 0.0, 0.5], [0.0, -1.0, -0.5],
        [0.0, 0.0, 0.8], [0.0, 0.0, -0.8],
    ])

    assert 4 in symmetry._improper_orders(
        ["C"] * 6, coordinates, np.asarray([0.0, 0.0, 1.0]), 0.01,
    )


def test_symmetry_distinguishes_an_icosahedral_structure():
    from oqp_studio import symmetry

    phi = (1 + 5 ** 0.5) / 2
    vertices = []
    for first in (-1.0, 1.0):
        for second in (-phi, phi):
            vertices.extend([(0.0, first, second), (first, second, 0.0), (second, 0.0, first)])
    xyz = "\n".join(f"C {x:.12f} {y:.12f} {z:.12f}" for x, y, z in vertices)

    assert symmetry.analyze(xyz, tolerance=0.001)["point_group"] == "Ih"


def test_symmetry_finds_face_centered_fivefold_axes_in_c60():
    from itertools import product

    from oqp_studio import symmetry

    phi = (1 + 5 ** 0.5) / 2
    seeds = [
        (0.0, 1.0, 3 * phi),
        (1.0, 2 + phi, 2 * phi),
        (phi, 2.0, 2 * phi + 1),
    ]
    vertices = set()
    for seed_index, seed in enumerate(seeds):
        for signs in product((-1.0, 1.0), repeat=3):
            signed = tuple(value * sign for value, sign in zip(seed, signs))
            for shift in range(3):
                vertex = signed[shift:] + signed[:shift]
                vertices.add(tuple(round(value, 12) for value in vertex))
    assert len(vertices) == 60
    xyz = "\n".join(f"C {x} {y} {z}" for x, y, z in sorted(vertices))

    assert symmetry.analyze(xyz, tolerance=0.001)["point_group"] == "Ih"
    reversed_xyz = "\n".join(
        f"C {x} {y} {z}" for x, y, z in sorted(vertices, reverse=True)
    )
    assert symmetry.analyze(reversed_xyz, tolerance=0.001)["point_group"] == "Ih"


def test_symmetry_falls_back_to_c1_and_centers_aligned_coordinates():
    import numpy as np

    from oqp_studio import symmetry

    result = symmetry.analyze(
        "C 0.1 0.2 0.3\nN 1.0 0.1 -0.2\nO -0.2 1.3 0.4\nH 0.5 -0.4 1.7\n",
        tolerance=0.001,
    )
    aligned = np.asarray([atom[1:] for atom in result["aligned_atoms"]])
    weights = np.asarray([12.011, 14.007, 15.999, 1.008], dtype=float)

    assert result["point_group"] == "C1"
    assert np.linalg.norm(np.average(aligned, axis=0, weights=weights)) < 1.0e-10
