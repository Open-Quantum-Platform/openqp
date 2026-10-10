#!/usr/bin/env python3
# Dump grid+basis+density+native-ref for a w8 water cluster (matches the
# worker/g4p timing baseline scale) at the SCF-default grid (96x302), SLATER.
import os, sys, json, subprocess
import numpy as np

WD   = "/bighome/cheolho.choi/clone/mixedprec/libintRot/sessions/20260617_ownxc/wd"
LANE = "/bighome/cheolho.choi/clone/mixedprec/libintRot/sessions/20260617_ownxc/oqp_root_lane"
PYOQP= "/bighome/cheolho.choi/openqp_routec/pyoqp"
os.chdir(WD)
FUNC = os.environ.get("FUNC","SLATER")
DUMP = os.path.join(WD,"dump_w8_"+FUNC)
os.makedirs(DUMP, exist_ok=True)

def water(ox,oy,oz):
    r=0.9572; th=np.deg2rad(104.52)
    return [(8,(ox,oy,oz)),(1,(ox+r,oy,oz)),
            (1,(ox+r*np.cos(th),oy+r*np.sin(th),oz))]
atoms=[]; aa=2.8
for i in range(2):
  for j in range(2):
    for k in range(2):
      atoms+=water(i*aa,j*aa,k*aa)
SYSTEM="".join(f"\n {z}  {x:.9f} {y:.9f} {zz:.9f}" for z,(x,y,zz) in atoms)

cfg={'input':{'system':SYSTEM,'charge':'0','runtype':'energy','basis':'cc-pvdz',
              'method':'hf','functional':FUNC,'d4':'False'},
     'guess':{'type':'huckel'},
     'dftgrid':{'rad_npts':'96','ang_npts':'302','pruned':'none'},
     'scf':{'type':'rhf','multiplicity':'1','maxit':'50','conv':'1.0e-9',
            'save_molden':'False','incremental':'False'}}
open("_run_w8.py","w").write(
'''import sys,json
from oqp.pyoqp import Runner
c=json.load(open(sys.argv[1]))
r=Runner(project="w8",input_dict=c,log=sys.argv[2],silent=1,usempi=False)
r.run()
print(f"EFINAL {r.mol.mol_energy.energy:.12f}")
''')
json.dump(cfg,open("_cfg_w8.json","w"))
env=dict(os.environ)
env["OPENQP_ROOT"]=LANE
env["PYTHONPATH"]=PYOQP+(":"+env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
env["OQP_OWNXC_DUMP"]=DUMP
env.pop("OQP_ROUTEC_XC_LIB",None)
env["OMP_NUM_THREADS"]="16";env["OPENBLAS_NUM_THREADS"]="16";env["MKL_NUM_THREADS"]="16"
rr=subprocess.run([sys.executable,"_run_w8.py","_cfg_w8.json",f"dump_w8_{FUNC}.log"],
                  capture_output=True,text=True,env=env)
print(rr.stdout[-800:])
if rr.returncode!=0: print("ERR",rr.stderr[-1500:])
print("DUMP",DUMP)
for f in sorted(os.listdir(DUMP)): print("  ",f,os.path.getsize(os.path.join(DUMP,f)))
