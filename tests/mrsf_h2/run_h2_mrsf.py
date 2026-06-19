#!/usr/bin/env python3
"""Run the H2 ROHF-MRSF regression inputs and print energies.
Verifies the zero-closed-shell fix: each must complete without the DGEMM crash.
Run: .venv/bin/python tests/mrsf_h2/run_h2_mrsf.py"""
import os
import sys

from oqp.pyoqp import Runner

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = ["h2_mrsfcis_rohf", "h2_mrsftddft_rohf"]


def main():
    failures = 0
    for case in CASES:
        inp = os.path.join(HERE, case + ".inp")
        log = os.path.join(HERE, case + ".log")
        print(f"\n===== {case} =====", flush=True)
        try:
            runner = Runner(project=case, input_file=inp, log=log, silent=1)
            runner.run()
            energy = runner.results().get("energy")
            print(f"  OK  energy={energy}", flush=True)
        except Exception as exc:
            failures += 1
            print(f"  FAIL  {type(exc).__name__}: {exc}", flush=True)

    if failures:
        print(f"\n{failures} H2 MRSF regression case(s) failed", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
