#!/usr/bin/env python3
"""Build matched-geometry inputs for the C GPU DF-RHF driver (routec_scf_solve),
for gpu4pyscf's water_cluster(n) molecule (the honest-comparison geometry).

Produces in <outdir>:
  B_w{n}.bin           dense B tensor, header {naux,nbf} int32 + float64 (naux,nbf,nbf)
                       (eigh-whitened V^-1/2, def2-universal-jkfit)
  H_w{n}.npy S_w{n}.npy  raw pyscf-frame 1e matrices (cc-pVDZ spherical)
  guess_minao_w{n}.bin packed lower-tri minao density (OpenQP t=i*(i+1)/2+j order)
  meta_w{n}.env        NOCC / ENUC / NAO / NAUX / EREF (shell-sourceable)

pyscf is used for INPUT-FILE GENERATION ONLY (same role as June's prep_w24_cdriver.py);
the C SCF itself is pyscf-free.  EREF = pyscf DF-RHF energy of the same molecule/aux —
the correctness bar for the driver's converged energy.
Usage: prep_matched.py <nwat> <outdir>
"""
import numpy as np, struct, os, sys, itertools


def water_cluster(n):
    r=0.9572; th=np.deg2rad(104.52); a=2.8; atoms=[]; c=0
    for (i,j,k) in itertools.product(range(4),repeat=3):
        if c>=n: break
        ox,oy,oz=i*a,j*a,k*a
        atoms+=[("O",(ox,oy,oz)),("H",(ox+r,oy,oz)),("H",(ox+r*np.cos(th),oy+r*np.sin(th),oz))]; c+=1
    return atoms

def build_df_tensor(nwat, outdir):
    """Generate the historical water-cluster DF reference files explicitly.

    Requires PySCF only when called. Returns the B tensor file path.
    """
    from pyscf import gto, scf, df
    n = int(nwat)
    if not 1 <= n <= 64:
        raise ValueError("nwat must be between 1 and 64")
    OUT = os.fspath(outdir)
    os.makedirs(OUT, exist_ok=True)
    mol = gto.M(atom=water_cluster(n), basis="cc-pvdz", unit="Angstrom", cart=False, verbose=0)
    nbf = mol.nao; nocc = mol.nelectron//2; enuc = mol.energy_nuc()
    S = mol.intor("int1e_ovlp"); H = mol.intor("int1e_kin") + mol.intor("int1e_nuc")
    np.save(f"{OUT}/H_w{n}.npy", H); np.save(f"{OUT}/S_w{n}.npy", S)
    print(f"[w{n}] nao={nbf} nocc={nocc} Enuc={enuc:.12f}")

    # B tensor: 3c integrals whitened by V^-1/2 (eigh, relative threshold), dense layout
    auxmol = df.addons.make_auxmol(mol, "def2-universal-jkfit")
    T3 = df.incore.aux_e2(mol, auxmol)              # (nao, nao, naux)
    V  = auxmol.intor("int2c2e")
    w, U = np.linalg.eigh(V)
    keep = w > 1e-13*w.max()
    Wm = U[:, keep] / np.sqrt(w[keep])              # (naux_full, naux_kept)
    Bt = np.tensordot(T3, Wm, axes=(2, 0))          # (nao, nao, naux_kept)
    B  = np.ascontiguousarray(np.moveaxis(Bt, 2, 0))
    naux = B.shape[0]
    with open(f"{OUT}/B_w{n}.bin", "wb") as f:
        f.write(struct.pack('ii', naux, nbf))
        B.astype(np.float64).tofile(f)
    print(f"[w{n}] B: naux_kept={naux}/{V.shape[0]}  {B.nbytes/1e9:.2f} GB")
    del T3, Bt, B

    # minao guess, packed lower-tri in OpenQP order
    D = scf.RHF(mol).get_init_guess(key="minao")
    print(f"[w{n}] guess tr(SD)={np.trace(S@D):.4f} (target {mol.nelectron})")
    ntri = nbf*(nbf+1)//2; out = np.zeros(ntri)
    for i in range(nbf):
        b = i*(i+1)//2; out[b:b+i+1] = D[i, :i+1]
    out.astype(np.float64).tofile(f"{OUT}/guess_minao_w{n}.bin")

    # correctness bar: pyscf CPU DF-RHF on the same molecule/aux
    mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit")
    mf.conv_tol = 1e-8; mf.verbose = 0
    eref = mf.kernel()
    print(f"[w{n}] pyscf DF-RHF reference E = {eref:.10f} ({mf.cycles if hasattr(mf,'cycles') else '?'} it)")

    with open(f"{OUT}/meta_w{n}.env", "w") as f:
        f.write(f"NAO={nbf}\nNOCC={nocc}\nENUC={enuc:.12f}\nNAUX={naux}\nEREF={eref:.10f}\n")
    print(f"[w{n}] DONE -> {OUT}")
    return os.path.join(OUT, f"B_w{n}.bin")


if __name__ == "__main__":
    build_df_tensor(int(sys.argv[1]), sys.argv[2])
