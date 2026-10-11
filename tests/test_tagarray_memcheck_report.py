"""Reject false lifetime qualification without importing/running native OpenQP."""
import importlib.util
from pathlib import Path

import pytest


def _reporter():
    path = Path(__file__).parent / "native" / "run_tagarray_memcheck.py"
    spec = importlib.util.spec_from_file_location("tagarray_reporter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _xml(path, *, finished=True, owner=False, kind="Leak_StillReachable"):
    frames = ""
    if owner:
        frames = """<frame><ip>0x1</ip><obj>/installed/liboqp.so</obj>
        <fn>TA_Container_new</fn></frame><frame><ip>0x2</ip>
        <obj>/installed/liboqp.so</obj><fn>occupied_cofactors_persistent_selftest</fn></frame>"""
    else:
        frames = "<frame><ip>0x3</ip><obj>/dependency/lib.so</obj><fn>dependency_init</fn></frame>"
    status = "FINISHED" if finished else "RUNNING"
    path.write_text(f"""<valgrindoutput><protocoltool>memcheck</protocoltool>
      <status><state>{status}</state></status><error><kind>{kind}</kind>
      <xwhat><text>102923 bytes retained</text><leakedbytes>102923</leakedbytes>
      <leakedblocks>56</leakedblocks></xwhat><stack>{frames}</stack></error></valgrindoutput>""")
    return path


def _arm(reporter, tmp_path, control="omit-cleanup", *, owner=True, raw=3, vg=97):
    status, checks = reporter.CONTROLS[control]
    line = f"Cofactor status={status} checks={','.join(map(str, checks))}\n"
    line += "TagArray check failed at line 56: status == 0\n"
    return {"control": control, "raw_returncode": raw, "valgrind_returncode": vg,
            "raw_output": line, "valgrind_output": line,
            "memory": reporter.memory_report(_xml(tmp_path / "report.xml", owner=owner))}


@pytest.mark.parametrize("control", ["positive", "missing-owner", "recreated-owner", "missing-record"])
def test_clean_owner_arms_reject_possibly_lost_memory(tmp_path, control):
    reporter = _reporter()
    if control == "positive":
        arm = {"control": control, "raw_returncode": 0, "valgrind_returncode": 0}
    else:
        arm = _arm(reporter, tmp_path, control=control, owner=False, vg=3)
    arm["memory"] = reporter.memory_report(
        _xml(tmp_path / "possibly-lost.xml", kind="Leak_PossiblyLost"))
    # Simulate even an erroneous successful child exit: parsed losses must
    # independently reject the clean-owner contract.
    assert arm["memory"]["non_leak_errors"] == 0
    assert arm["memory"]["leaks"]["Leak_DefinitelyLost"]["bytes"] == 0
    assert arm["memory"]["leaks"]["Leak_IndirectlyLost"]["bytes"] == 0
    assert arm["memory"]["leaks"]["Leak_PossiblyLost"]["bytes"] == 102923
    assert reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})["pass"] is False


def test_reachable_counts_without_owner_stacks_are_inconclusive(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path, owner=False)
    result = reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})
    assert arm["memory"]["leaks"]["Leak_StillReachable"] == {"bytes": 102923, "blocks": 56}
    assert result["guard"] == "CONFIRMED"
    assert result["memory_sensitivity"] == "INCONCLUSIVE"
    assert result["pass"] is False


def test_symbol_attributed_delta_confirms_memory_sensitivity(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path)
    result = reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})
    assert result["pass"] is True
    assert result["memory_sensitivity"] == "CONFIRMED"
    frames = result["owner_retention"]["records"][0]["stacks"][0]
    assert frames[0]["obj"] == "/installed/liboqp.so"


def test_same_positive_retention_does_not_prove_missing_cleanup(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path)
    result = reporter.verdict(arm, reporter.owner_retention(arm["memory"]))
    assert result["memory_sensitivity"] == "INCONCLUSIVE"
    assert result["pass"] is False


@pytest.mark.parametrize("raw,vg", [(0, 97), (3, 0), (-11, 97), (2, 2)])
def test_nonmatching_child_failures_cannot_qualify_guard(tmp_path, raw, vg):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path, raw=raw, vg=vg)
    assert reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})["guard"] == "FAILED"


def test_wrong_native_status_cannot_qualify_guard(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path)
    arm["raw_output"] = arm["raw_output"].replace("status=7", "status=4")
    assert reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})["pass"] is False


def test_raw_guard_rejection_alone_cannot_qualify_valgrind(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path)
    arm["valgrind_output"] = "dlopen failure\n"
    assert reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})["pass"] is False


def test_native_memory_error_is_not_expected_negative_control(tmp_path):
    reporter = _reporter()
    arm = _arm(reporter, tmp_path)
    arm["memory"]["non_leak_errors"] = 1
    assert reporter.verdict(arm, {"bytes": 0, "blocks": 0, "records": []})["pass"] is False


def test_incomplete_or_missing_diagnostics_fail(tmp_path):
    reporter = _reporter()
    with pytest.raises(ValueError, match="did not finish"):
        reporter.memory_report(_xml(tmp_path / "unfinished.xml", finished=False))
    with pytest.raises(FileNotFoundError):
        reporter.memory_report(tmp_path / "missing.xml")


def test_child_timeout_preserves_partial_output_and_fails(monkeypatch, tmp_path):
    import json
    import subprocess
    reporter = _reporter()
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 120, output=b"partial native output\n")
    monkeypatch.setattr(reporter.subprocess, "run", timeout)
    outcome = reporter.run_child(["not-executed"], {}, tmp_path, "probe")
    assert outcome["timed_out"] is True
    assert outcome["returncode"] is None
    assert (tmp_path / "probe.log").read_text() == "partial native output\n"
    assert json.loads((tmp_path / "probe-outcome.json").read_text()) == outcome
    assert reporter.finish(tmp_path, {}, [{"verdict": {"pass": False}}]) == 1


@pytest.mark.parametrize("broken_xml", [False, True])
def test_main_isolates_controls_and_preserves_outcomes(monkeypatch, tmp_path, broken_xml):
    import json
    import subprocess
    reporter = _reporter()
    library = tmp_path / "not-a-native-library"
    library.write_bytes(b"synthetic input: no native execution")
    output = tmp_path / "diagnostics"
    monkeypatch.setattr(reporter.sys, "argv", ["reporter", "--library", str(library), "--output", str(output)])
    monkeypatch.setenv(reporter.CONTROL_ENV, "missing-record")
    monkeypatch.setattr(reporter.subprocess, "check_output", lambda argv, **kw:
                        "valgrind-3.22.0" if argv[0] == "valgrind" else "mock compiler")
    calls = []
    def run(argv, **kwargs):
        if argv[0] == "cc":
            return subprocess.CompletedProcess(argv, 0, "")
        calls.append((argv, kwargs["env"].copy()))
        control = kwargs["env"].get(reporter.CONTROL_ENV)
        if argv[0] == "valgrind":
            assert "--errors-for-leak-kinds=definite,indirect,possible" in argv
            xml = Path(next(arg.split("=", 1)[1] for arg in argv if arg.startswith("--xml-file=")))
            if broken_xml:
                xml.write_text("<truncated>")
            else:
                _xml(xml, owner=control == "omit-cleanup")
        if control:
            status, checks = reporter.CONTROLS[control]
            line = f"Cofactor status={status} checks={','.join(map(str, checks))}\n"
            line += "TagArray check failed at line 56: status == 0\n"
            code = 97 if argv[0] == "valgrind" and control == "omit-cleanup" else 3
        else:
            line, code = "retained child output\n", 0
        return subprocess.CompletedProcess(argv, code, line)
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.setattr(reporter.subprocess, "run", run)
    assert reporter.main() == (1 if broken_xml else 0)
    arm = json.loads((output / "omp1" / "positive" / "result.json").read_text())
    assert arm["raw_returncode"] == 0 and arm["valgrind_returncode"] == 0
    result = json.loads((output / "result.json").read_text())
    if broken_xml:
        assert all(reporter.CONTROL_ENV not in env for _, env in calls)
        assert arm["diagnostics_error"].startswith("ParseError:")
        assert arm["verdict"]["pass"] is False and result["pass"] is False
    else:
        assert result["pass"] is True and len(result["arms"]) == 10
        assert len(calls) == 20  # independent raw and Valgrind child for every arm
        for threads in (1, 4):
            selected = [env for _, env in calls if env["OMP_NUM_THREADS"] == str(threads)]
            assert len(selected) == 10
            assert [env.get(reporter.CONTROL_ENV) for env in selected] == [
                value for value in (None, *reporter.CONTROLS) for _ in range(2)]
            assert all(env["MKL_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == "1" for env in selected)
    assert reporter.os.environ[reporter.CONTROL_ENV] == "missing-record"
