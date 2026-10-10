"""Validation: openqp-gpu RHF DF-SCF vs a pyscf DF-RHF reference on a small
water cluster.  Requires a built libopenqp_gpu.so (OPENQP_GPU_LIB) and pyscf.
Run on a CUDA host.  Passes if |E - E_ref| < 1e-8.
"""
import os, numpy as np
from openqp_gpu import solve_rhf, build_df_tensor

def test_water_cluster(n=4):
    from pyscf import gto, scf
    import itertools
    def wc(n):
        r=0.9572; th=np.deg2rad(104.52); a=2.8; at=[]; c=0
        for i,j,k in itertools.product(range(4),repeat=3):
            if c>=n: break
            o=(i*a,j*a,k*a)
            at+=[("O",o),("H",(o[0]+r,o[1],o[2])),
                 ("H",(o[0]+r*np.cos(th),o[1]+r*np.sin(th),o[2]))]; c+=1
        return at
    mol=gto.M(atom=wc(n),basis="cc-pvdz",unit="Angstrom",cart=False,verbose=0)
    os.environ["OPENQP_GPU_B"]=build_df_tensor(mol, "/tmp/B_test.bin")
    H=mol.intor("int1e_kin")+mol.intor("int1e_nuc"); S=mol.intor("int1e_ovlp")
    E,cyc=solve_rhf(H,S,mol.nelectron//2,mol.energy_nuc(),guess="gwh")
    Eref=scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit").kernel()
    assert abs(E-Eref)<1e-8, f"E={E} Eref={Eref} d={abs(E-Eref):.2e}"
    print(f"OK  E={E:.10f}  ref={Eref:.10f}  cycles={cyc}")

if __name__=="__main__":
    test_water_cluster()
