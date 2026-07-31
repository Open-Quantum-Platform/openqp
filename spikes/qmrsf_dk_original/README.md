# QMRSF-DK — production code

The QMRSF-DK engine used for every QMRSF-DK number in the paper. These five
files are byte-identical copies of the modules in
`~/Projects/QMRSF-DK-pyscf/` (md5 below); nothing was renamed or edited.

| file | role |
|---|---|
| `qmrsf_dk_pyscf.py` | quintet ROKS reference, active-space selection, active-space integral transform |
| `qmrsf_dk_paper_direct.py` | **the production builder** — `build_paper_matrix_f0_direct` |
| `qmrsf_dk_paper_f0.py` | F(0) / one-electron part of the dressed active-space operator |
| `qmrsf_dk_paper.py` | spin-adapted matrix definitions used by the builder |
| `qmrsf_dk_matrices_ref.py` | configuration labels, spin values, reference matrices |

Requirements: python 3 + numpy + pyscf (2.x). Put the folder on `sys.path`
(or run from inside it) — the modules import each other by plain module name.

## The protocol (this exact sequence produces the published numbers)

```python
import qmrsf_dk_pyscf as Q
import qmrsf_dk_paper_direct as D
import numpy as np

atom  = "\n".join(open("mol.xyz").read().splitlines()[2:])   # xyz coordinate block
basis = "cc-pvtz"        # cc-pVDZ was used for the polyene dimers
c_HF  = 0.5              # exact-exchange fraction of the dressed operator

mf, _ = Q.build_reference(atom, basis, xc="bhhlyp")   # S=2, Ms=+2 quintet ROKS
somo, ncore = Q.active_space(mf)                      # 4 energy-ordered SOMOs = CAS(4,4)
h, eri, ec, ci, mc = Q.active_integrals(mf, ncore, somo)

dk = D.build_paper_matrix_f0_direct(mf, ncore, somo, c_HF, h_act=h, eri_act=eri)

w, V = np.linalg.eigh(np.atleast_2d(dk["A_singlet"]))   # 20 singlets
wT   = np.linalg.eigvalsh(np.atleast_2d(dk["A_triplet"]))  # 15 triplets
Aq   = float(np.ravel(dk["A_quintet"])[0])              # 1 quintet

E = [mf.e_tot + (x - Aq) for x in w]                    # total energies of the singlets
VEE_eV = [(x - E[0]) * 27.211386245988 for x in E]      # vertical excitation energies
```

Two points that must not be changed, or the numbers will not match the paper:

1. **Use `build_paper_matrix_f0_direct`.** It assembles the hardcoded
   spin-adapted blocks directly; the 36-determinant Hamiltonian is never formed
   and no unitary transformation is applied. Other builders in the file history
   are development versions and do not reproduce the published values.
2. **Use the energy law `E_I = E_SCF(quintet) + (lambda_I - A_quintet)`**, not
   `e_core + lambda_I`.

State character comes from the eigenvector weights over `dk["labels"]` filtered
by `dk["Sval"] == 0` (singlets), squared and expressed in percent.

## Settings used in the paper

BHHLYP quintet reference; c_HF = 0.5; cc-pVTZ (cc-pVDZ for the polyene dimers);
SCF converged to 1e-8 Ha on the energy and 1e-6 on the RMS density, fine
exchange–correlation grid. The four SOMOs of the aufbau quintet define the
active space; systems whose aufbau quintet does not span the intended valence
space (cyclobutadiene is the only such case encountered) need a maximum-overlap
constrained reference — build it yourself and pass the converged `mf` into
`Q.active_space` / `Q.active_integrals`, the rest of the sequence is unchanged.

## Integrity

```
004d9f97f96f0b62e11fc7275eb5ef4f  qmrsf_dk_matrices_ref.py
58f9ddeec6c58eed77589c2c21c32049  qmrsf_dk_paper_direct.py
a4beaf9732b36fbc3e53d7216240660d  qmrsf_dk_paper_f0.py
97a732359b161e6996244a01e00ae264  qmrsf_dk_paper.py
a0b237d2d70618b8b9ca634feb5e3ad1  qmrsf_dk_pyscf.py
```
