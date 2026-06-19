#!/usr/bin/env python3
"""Run the four H2 MRSF regression inputs (CIS/TDDFT x ROHF/UHF) and print energies.
Verifies the zero-closed-shell fix: each must complete without the DGEMM crash.
Run: .venv/bin/python tests/mrsf_h2/run_h2_mrsf.py"""
import os, sys, glob
HERE=os.path.dirname(os.path.abspath(__file__))
from oqp.pyoqp import Runner
cases=['h2_mrsfcis_rohf','h2_mrsftddft_rohf']
for c in cases:
    inp=os.path.join(HERE,c+'.inp'); log=os.path.join(HERE,c+'.log')
    print(f"\n===== {c} =====", flush=True)
    try:
        r=Runner(project=c, input_file=inp, log=log, silent=1)
        r.run(); res=r.results()
        e=res.get('energy')
        print(f"  OK  energy={e}", flush=True)
    except Exception as ex:
        print(f"  FAIL  {type(ex).__name__}: {ex}", flush=True)
