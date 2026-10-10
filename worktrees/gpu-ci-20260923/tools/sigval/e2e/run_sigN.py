#!/usr/bin/env python3
"""(H2O)N cc-pVDZ MRSF seam test: native vs GPU sigma-session, geometry
identical to stage2_gpu_inloop_e2e.py water grid NX x NY x NZ."""
import argparse, json, os, subprocess, sys
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument("leg",choices=["native","sig"])
ap.add_argument("--grid",nargs=3,type=int,required=True)
ap.add_argument("--so",default=""); ap.add_argument("--b",default="")
ap.add_argument("--nstate",type=int,default=3)
ap.add_argument("--out",required=True)
a=ap.parse_args()

def water(ox,oy,oz):
    r=0.9572; th=np.deg2rad(104.52)
    return [(8,(ox,oy,oz)),(1,(ox+r,oy,oz)),
            (1,(ox+r*np.cos(th),oy+r*np.sin(th),oz))]
atoms=[]; aa=2.8; NX,NY,NZ=a.grid
for i in range(NX):
    for j in range(NY):
        for k in range(NZ): atoms+=water(i*aa,j*aa,k*aa)
SYSTEM="".join(f"\n   {z}   {x:.9f} {y:.9f} {zz:.9f}" for z,(x,y,zz) in atoms)

CFG={'input':{'system':SYSTEM,'charge':'0','runtype':'energy',
              'basis':'cc-pvdz','method':'tdhf','functional':'bhhlyp','d4':'False'},
     'guess':{'type':'huckel'},
     'scf':{'type':'rohf','multiplicity':'3','maxit':'200','conv':'1.0e-9',
            'save_molden':'False','incremental':'False'},
     'tdhf':{'type':'mrsf','multiplicity':'1','nstate':str(a.nstate),
             'conv':'1.0e-8','maxit':'100'}}
RUNSCRIPT=r'''
import json,sys
import numpy as np
from oqp.pyoqp import Runner
cfg=json.load(open(sys.argv[1]))
r=Runner(project='sigN',input_dict=cfg,log=sys.argv[3],silent=1,usempi=False)
r.run(); mol=r.mol
d={'E_scf':mol.mol_energy.energy}
try: d['td']=np.asarray(mol.data['OQP::td_energies'],dtype=float)
except Exception as e: d['td_err']=str(e)
np.savez(sys.argv[2],**d); print('RUN_OK')
'''
json.dump(CFG,open("_sN.json","w")); open("_runN.py","w").write(RUNSCRIPT)
env=dict(os.environ)
for k in ("OQP_ROUTEC_SIG","OQP_ROUTEC_LIB","OQP_ROUTEC_B"): env.pop(k,None)
if a.leg=="sig": env["OQP_ROUTEC_SIG"]=a.so; env["OQP_ROUTEC_B"]=a.b
rr=subprocess.run([sys.executable,"_runN.py","_sN.json",a.out,a.out+".log"],
                  capture_output=True,text=True,env=env)
for ln in rr.stderr.splitlines():
    if "routec_sig" in ln or "sig-gpu" in ln: print("  [stderr]",ln)
if "RUN_OK" not in rr.stdout:
    print("FAILED\n",rr.stdout[-1500:],rr.stderr[-1500:]); sys.exit(1)
z=dict(np.load(a.out))
print(f"  E_scf={float(z['E_scf']):.9f}")
if 'td' in z: print("  exc:", " ".join(f"{v:.9f}" for v in z['td'][:a.nstate]))
else: print("  td MISSING:",z.get('td_err'))
