#!/usr/bin/env python3
"""Check installed TagArray lifetimes and isolated cofactor fault controls.

Linux CI requires Valgrind >= 3.22. Preserve full allocation stacks, raw child
exits and Valgrind exits independently; a rejected fault is not a leak proof.
No copied dependency, allocator replacement, or suppression is used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET

CONTROL_ENV = "OQP_COFACTOR_LIFETIME_CONTROL"
CONTROLS = {
    "missing-owner": (2, [0, 0, 0, 0, 0, 2]),
    "recreated-owner": (3, [1, 0, 0, 0, 0, 2]),
    "missing-record": (4, [1, 1, 0, 0, 0, 2]),
    "omit-cleanup": (7, [1, 4, 16, 4, 4, 2]),
}
LEAK_KINDS = {"Leak_DefinitelyLost", "Leak_IndirectlyLost",
              "Leak_PossiblyLost", "Leak_StillReachable"}


def memory_report(path):
    """Keep every record and frame; never infer owner identity from heap counts."""
    root = ET.parse(path).getroot()  # absent/truncated diagnostics fail, never skip
    if root.tag != "valgrindoutput" or root.findtext("protocoltool") != "memcheck":
        raise ValueError("Not a Valgrind memcheck XML report")
    if not any(s.findtext("state") == "FINISHED" for s in root.findall("status")):
        raise ValueError("Valgrind report did not finish")
    records = []
    for error in root.findall("error"):
        stacks = [[{key: frame.findtext(key) for key in
                    ("ip", "obj", "fn", "dir", "file", "line")}
                   for frame in stack.findall("frame")]
                  for stack in error.findall("stack")]
        records.append({"kind": error.findtext("kind"),
                        "bytes": int(error.findtext("xwhat/leakedbytes", "0")),
                        "blocks": int(error.findtext("xwhat/leakedblocks", "0")),
                        "text": error.findtext("xwhat/text") or error.findtext("what"),
                        "stacks": stacks})
    return {"records": records,
            "non_leak_errors": sum(r["kind"] not in LEAK_KINDS for r in records),
            "leaks": {kind: {"bytes": sum(r["bytes"] for r in records if r["kind"] == kind),
                             "blocks": sum(r["blocks"] for r in records if r["kind"] == kind)}
                      for kind in sorted(LEAK_KINDS)},
            "reachable_owners": "UNCLASSIFIED: inspect retained allocation stacks and object identities"}


def owner_retention(report):
    """Identify only stacks naming both the shim and a TagArray allocation site."""
    records = []
    for record in report["records"]:
        if record["kind"] not in LEAK_KINDS:
            continue
        for stack in record["stacks"]:
            names = [frame["fn"] or "" for frame in stack]
            shim = any("occupied_cofactors_persistent_selftest" in name for name in names)
            allocation = any("TA_Container_new" in name or "tagarray::hidden::Record::Record" in name
                             or "tagarray::Container::create" in name for name in names)
            if shim and allocation:
                records.append(record)
                break
    return {"bytes": sum(r["bytes"] for r in records),
            "blocks": sum(r["blocks"] for r in records), "records": records}


def first_status(output):
    match = re.search(r"^Cofactor status=(-?\d+) checks=((?:-?\d+,){5}-?\d+)$",
                      output, re.MULTILINE)
    if match is None:
        return None
    return int(match[1]), list(map(int, match[2].split(",")))


def verdict(arm, positive_owner):
    report = arm["memory"]
    if arm["control"] == "positive":
        losses = report["leaks"]
        return {"guard": "NOT_APPLICABLE", "memory_sensitivity": "NOT_APPLICABLE",
                "pass": arm["raw_returncode"] == 0 and arm["valgrind_returncode"] == 0
                and report["non_leak_errors"] == 0
                and losses["Leak_DefinitelyLost"]["bytes"] == 0
                and losses["Leak_IndirectlyLost"]["bytes"] == 0
                and losses["Leak_PossiblyLost"]["bytes"] == 0}
    expected = CONTROLS[arm["control"]]
    # A raw child rejection and a Valgrind override must remain distinguishable.
    guard = (arm["raw_returncode"] == 3 and arm["valgrind_returncode"] in (3, 97)
             and first_status(arm["raw_output"]) == expected
             and first_status(arm["valgrind_output"]) == expected
             and "status == 0" in arm["raw_output"]
             and "status == 0" in arm["valgrind_output"])
    result = {"guard": "CONFIRMED" if guard else "FAILED",
              "memory_sensitivity": "NOT_APPLICABLE", "pass": guard}
    if arm["control"] == "omit-cleanup":
        retained = owner_retention(report)
        confirmed = (retained["bytes"] > positive_owner["bytes"]
                     and retained["blocks"] > positive_owner["blocks"])
        result.update(memory_sensitivity="CONFIRMED" if confirmed else "INCONCLUSIVE",
                      owner_retention=retained, positive_owner_retention=positive_owner,
                      pass_=guard and confirmed and report["non_leak_errors"] == 0)
        result["pass"] = result.pop("pass_")
    else:
        # These faults clean their owner; incidental memory errors are never accepted.
        result["pass"] = guard and arm["valgrind_returncode"] == 3 and report["non_leak_errors"] == 0
        result["pass"] &= all(report["leaks"][kind]["bytes"] == 0
                              for kind in ("Leak_DefinitelyLost", "Leak_IndirectlyLost", "Leak_PossiblyLost"))
    return result


def run_child(argv, env, folder, name):
    """Persist terminal/timeout/launch evidence before interpreting diagnostics."""
    outcome = {"argv": argv, "returncode": None, "timed_out": False, "launch_error": None}
    try:
        child = subprocess.run(argv, env=env, timeout=120, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outcome.update(returncode=child.returncode, output=child.stdout)
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode("utf-8", errors="replace")
        outcome.update(timed_out=True, output=partial)
    except OSError as error:
        outcome.update(launch_error=str(error), output="")
    (folder / f"{name}.log").write_text(outcome["output"])
    (folder / f"{name}-outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")
    return outcome


def finish(output, provenance, arms):
    result = {"provenance": provenance, "arms": arms,
              "pass": len(arms) == 2 * (1 + len(CONTROLS))
              and all(a["verdict"]["pass"] for a in arms),
              "scope": "isolated test-shim defects; reachable ownership requires stacks, no whole-engine leak claim"}
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return 0 if result["pass"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--output", type=Path, default=Path("tagarray-memory"))
    args = parser.parse_args()
    version = subprocess.check_output(["valgrind", "--version"], text=True).strip()
    numbers = re.search(r"(\d+)\.(\d+)", version)
    if not numbers or tuple(map(int, numbers.groups())) < (3, 22):
        raise RuntimeError(f"Valgrind >= 3.22 required for alignment checks: {version}")
    library = args.library
    if library is None:
        from oqp.runtime import library_path, resolve_oqp_root
        library = library_path(*resolve_oqp_root())
    if library is None or not library.is_file():
        raise FileNotFoundError(f"Installed liboqp not found: {library}")
    library = library.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("tagarray_lifetime.c")
    executable = output / "tagarray_lifetime"
    compiler = shlex.split(os.environ.get("CC", "cc"))
    subprocess.run(compiler + ["-std=c11", "-Wall", "-Wextra", "-Werror", "-g",
                              str(source), "-ldl", "-o", str(executable)], check=True)
    provenance = {"valgrind": version,
                  "compiler": subprocess.check_output(compiler + ["--version"], text=True).splitlines()[0],
                  "library": str(library), "library_sha256": hashlib.sha256(library.read_bytes()).hexdigest(),
                  "driver_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance), flush=True)
    arms = []
    for threads in (1, 4):
        positive_owner = None
        for control in ("positive", *CONTROLS):
            folder = output / f"omp{threads}" / control
            folder.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            env.update(OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads),
                       OMP_MAX_ACTIVE_LEVELS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       BLIS_NUM_THREADS="1")
            env.pop(CONTROL_ENV, None)  # inherited faults cannot contaminate another child
            if control != "positive":
                env[CONTROL_ENV] = control
            raw_argv = [str(executable), str(library)]
            vg_argv = ["valgrind", "--tool=memcheck", "--track-origins=yes", "--keep-debuginfo=yes",
                       "--leak-check=full", "--show-leak-kinds=all", "--num-callers=40",
                       "--leak-resolution=high", "--xml=yes", f"--xml-file={folder / 'valgrind.xml'}",
                       "--errors-for-leak-kinds=definite,indirect,possible", "--error-exitcode=97",
                       f"--log-file={folder / 'valgrind.log'}", *raw_argv]
            (folder / "commands.json").write_text(json.dumps({"raw": raw_argv, "valgrind": vg_argv,
                                                              "control_env": env.get(CONTROL_ENV), "threads": threads}, indent=2) + "\n")
            arm = {"control": control, "threads": threads, "memory": None,
                   "verdict": {"pass": False, "diagnostics": "INVALID"}}
            raw = run_child(raw_argv, env, folder, "probe-raw")
            arm.update(raw_returncode=raw["returncode"], raw_output=raw["output"], raw_child=raw)
            (folder / "result.json").write_text(json.dumps(arm, indent=2) + "\n")
            if raw["timed_out"] or raw["launch_error"]:
                arms.append(arm)
                return finish(output, provenance, arms)
            checked = run_child(vg_argv, env, folder, "probe")
            arm.update(valgrind_returncode=checked["returncode"], valgrind_output=checked["output"],
                       valgrind_child=checked)
            (folder / "result.json").write_text(json.dumps(arm, indent=2) + "\n")
            if checked["timed_out"] or checked["launch_error"]:
                arms.append(arm)
                return finish(output, provenance, arms)
            try:
                arm["memory"] = memory_report(folder / "valgrind.xml")
            except (OSError, ET.ParseError, ValueError) as error:
                arm["diagnostics_error"] = f"{type(error).__name__}: {error}"
                (folder / "result.json").write_text(json.dumps(arm, indent=2) + "\n")
                arms.append(arm)
                return finish(output, provenance, arms)
            if control == "positive":
                positive_owner = owner_retention(arm["memory"])
            arm["verdict"] = verdict(arm, positive_owner)
            (folder / "result.json").write_text(json.dumps(arm, indent=2) + "\n")
            arms.append(arm)
            print(json.dumps({"control": control, "threads": threads, "raw_returncode": raw["returncode"],
                              "valgrind_returncode": checked["returncode"], "verdict": arm["verdict"]}), flush=True)
    return finish(output, provenance, arms)


if __name__ == "__main__":
    sys.exit(main())
