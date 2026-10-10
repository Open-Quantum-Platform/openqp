#!/usr/bin/env python3
# Run a small DF-RHF/LDA OpenQP job on the ownxc lane with OQP_OWNXC_DUMP set,
# so the lane exporter writes grid.bin / basis.bin / cart.bin / dens.bin / ref.bin
# into the dump dir on the LAST SCF iteration (we read the final converged ones).
import os, sys, json, subprocess

WD   = "/bighome/cheolho.choi/clone/mixedprec/libintRot/sessions/20260617_ownxc/wd"
LANE = "/bighome/cheolho.choi/clone/mixedprec/libintRot/sessions/20260617_ownxc/oqp_root_lane"
PYOQP= "/bighome/cheolho.choi/openqp_routec/pyoqp"
os.makedirs(WD, exist_ok=True)
os.chdir(WD)

FUNC = os.environ.get("FUNC", "SLATER")
RAD  = os.environ.get("RAD", "96")
ANG  = os.environ.get("ANG", "302")
DUMP = os.environ.get("DUMPDIR", os.path.join(WD, "dump_"+FUNC))
os.makedirs(DUMP, exist_ok=True)

# small water molecule (Angstrom)
SYSTEM = """
 8   0.000000000  0.000000000  0.000000000
 1   0.957200000  0.000000000  0.000000000
 1  -0.239987656  0.926627480  0.000000000"""

cfg = {
  'input': {'system': SYSTEM, 'charge': '0', 'runtype': 'energy',
            'basis': 'cc-pvdz', 'method': 'hf',
            'functional': FUNC, 'd4': 'False'},
  'guess': {'type': 'huckel'},
  'dftgrid': {'rad_npts': RAD, 'ang_npts': ANG, 'pruned': 'none'},
  'scf': {'type': 'rhf', 'multiplicity': '1', 'maxit': '50',
          'conv': '1.0e-9', 'save_molden': 'False', 'incremental': 'False'},
}

RUNSCRIPT = r'''
import sys, json
from oqp.pyoqp import Runner
c = json.load(open(sys.argv[1]))
r = Runner(project='dump', input_dict=c, log=sys.argv[2], silent=1, usempi=False)
r.run()
print(f"EFINAL {r.mol.mol_energy.energy:.12f}")
'''
open("_run_dump.py","w").write(RUNSCRIPT)
json.dump(cfg, open("_cfg_dump.json","w"))

env = dict(os.environ)
env["OPENQP_ROOT"] = LANE
env["PYTHONPATH"]  = PYOQP + (":"+env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
env["OQP_OWNXC_DUMP"] = DUMP
env.pop("OQP_ROUTEC_XC_LIB", None)
env.pop("OQP_ROUTEC_SCF_LIB", None)
env["OMP_NUM_THREADS"]="16"; env["OPENBLAS_NUM_THREADS"]="16"; env["MKL_NUM_THREADS"]="16"

rr = subprocess.run([sys.executable,"_run_dump.py","_cfg_dump.json",f"dump_{FUNC}.log"],
                    capture_output=True, text=True, env=env)
print(rr.stdout[-2500:]);
if rr.returncode != 0: print("STDERR:", rr.stderr[-2500:])
print("DUMPDIR:", DUMP)
for f in sorted(os.listdir(DUMP)):
    print("  ", f, os.path.getsize(os.path.join(DUMP,f)))
