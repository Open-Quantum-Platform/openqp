#!/usr/bin/env python3
"""HF/DFT/MRSF gradient through OpenQP, native vs GPU seam (routec_grad2 / _mrsf)."""
import argparse, json, os, subprocess, sys
import numpy as np
ap=argparse.ArgumentParser()
ap.add_argument("kind",choices=["hfgrad","dftgrad","mrsfgrad"])
ap.add_argument("leg",choices=["native","gpu"])
ap.add_argument("--grid",nargs=3,type=int,required=True)
ap.add_argument("--so",default=""); ap.add_argument("--b",default="")
ap.add_argument("--deck",default=""); ap.add_argument("--tables",default="")
ap.add_argument("--out",required=True)
a=ap.parse_args()
def water(ox,oy,oz):
    r=0.9572; th=np.deg2rad(104.52)
    return [(8,(ox,oy,oz)),(1,(ox+r,oy,oz)),(1,(ox+r*np.cos(th),oy+r*np.sin(th),oz))]
atoms=[]; aa=2.8; NX,NY,NZ=a.grid
for i in range(NX):
    for j in range(NY):
        for k in range(NZ): atoms+=water(i*aa,j*aa,k*aa)
SYSTEM="".join(f"\n   {z}   {x:.9f} {y:.9f} {zz:.9f}" for z,(x,y,zz) in atoms)
inp={'input':{'system':SYSTEM,'charge':'0','basis':'cc-pvdz','runtype':'grad','d4':'False'},
     'guess':{'type':'huckel'}}
if a.kind=="mrsfgrad":
    inp['input']['method']='tdhf'; inp['input']['functional']='bhhlyp'
    inp['scf']={'type':'rohf','multiplicity':'3','maxit':'200','conv':'1.0e-9','incremental':'False'}
    inp['tdhf']={'type':'mrsf','multiplicity':'1','nstate':'3','conv':'1.0e-8'}
    inp['properties']={'grad':'1'}
else:
    inp['input']['method']='hf'
    if a.kind=="dftgrad": inp['input']['functional']='bhhlyp'
    inp['scf']={'type':'rhf','multiplicity':'1','maxit':'100','conv':'1.0e-9','incremental':'False'}
RUN=r'''
import json,sys,numpy as np
from oqp.pyoqp import Runner
cfg=json.load(open(sys.argv[1]))
r=Runner(project='grad',input_dict=cfg,log=sys.argv[3],silent=1,usempi=False); r.run()
mol=r.mol; d={'E':float(mol.mol_energy.energy)}
try: d['grad']=np.asarray(mol.get_grad(),dtype=float).ravel()
except Exception as e: d['grad_err']=str(e)
np.savez(sys.argv[2],**d); print('RUN_OK')
'''
json.dump(inp,open("_g.json","w")); open("_rung.py","w").write(RUN)
env=dict(os.environ)
for k in ("OQP_ROUTEC_LIB","OQP_ROUTEC_XC_LIB","OQP_ROUTEC_GRAD_LIB","OQP_ROUTEC_B",
          "OQP_ROUTEC_GRAD_INP","OQP_ROUTEC_TABLES"): env.pop(k,None)
if a.leg=="gpu":
    env["OQP_ROUTEC_GRAD_LIB"]=a.so; env["OQP_ROUTEC_B"]=a.b
    env["OQP_ROUTEC_GRAD_INP"]=a.deck; env["OQP_ROUTEC_TABLES"]=a.tables
    env["OQP_ROUTEC_LIB"]=a.so; env["OQP_ROUTEC_XC_LIB"]=a.so   # so the SCF part also uses GPU
rr=subprocess.run([sys.executable,"_rung.py","_g.json",a.out,a.out+".log"],
                  capture_output=True,text=True,env=env)
for ln in rr.stderr.splitlines():
    if "routec" in ln.lower() and ("grad" in ln.lower() or "declin" in ln.lower() or "loaded" in ln.lower()): print("  [stderr]",ln)
if "RUN_OK" not in rr.stdout: print("FAILED:",rr.stdout[-1200:],rr.stderr[-1200:]); sys.exit(1)
z=dict(np.load(a.out)); print(f"  E={float(z['E']):.9f}", "grad" if "grad" in z else z.get("grad_err"))
