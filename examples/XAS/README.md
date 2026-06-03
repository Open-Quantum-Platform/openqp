This is examples for MRSF and delta-CHP-MRSF described in the reference below.
doi : https://doi.org/10.1021/acs.jctc.2c00746

-------------------------------------------------------------------------------------

* MRSF
- Use "ixcore"
To get core excitation from normal MRSF calculation, which does not include core-hole relaxation, one need only ixcore keyword.  
The ixcore keywords shift the fock elements to first get core excitations. By using this, we don't need to request many number of state (nstate keyword).
One can either use huckel guess or read MO from previous calculation.

In the example, Carbon K-edge is compted from state 52. The 1 -> 52 transition has energy of 285.55 eV and oscillator strength of 0.0495.  The state 52 has 2  ->   7 spin flip excitation from the reference state which results in 2 -> 8 (pi*) transition. 
One can see that Carbon K-edge of HCN system is 1s --> pi* transition with 285.55 eV and 0.0495 oscillator strength. 


* Important note for MRSF
Because of MRSF implementation, one can only obtain 1s --> LUMO transition with two singly occupied orbital which is single excitation.
All the rest 1s --> LUMO+n transitions have four singly occupied orbitals (two alpha in HOMO and LUMO, and two beta in core and LUMO+n), which are double excitations.

-------------------------------------------------------------------------------------

* delta-CHP-MRSF
- Use "ixcore, swapmo, rstctmo"
The MRSF utilizes Ms=+1 triplet state as a reference state, which usually put two singly occupied electons in HOMO and LUMO. However, one can also choose any other MOs.
The delta-CHP-MRSF includes core-hole relaxation by restricting one electron in core orbital with the help of MOM method.
For example, one can put singly occupied electrons in one of core orbital and LUMO.

One first need to swap MO (core and HOMO) by "swapmo=2,7" and restrict the MO ordering during SCF by "rstctmo=True" keyword. Once the SCF is converged, one can get excitations from core (HOMO) by "ixcore=7" keyword.

In the example, Carbon K-edge is computed from state 2, which is a 1s --> pi* transition. The ground state obtained from CHP-MRSF is usually poorly described, therefore, it is better to use the ground state from KS DFT (or MRSF-TDDFT) calculation.

* Important note for delta-CHP-MRSF
In current implementation, we don't get core (HOMO after using swapmo and rstctmo) -> LUMO excitation.
Therefore, one must check if core->LUMO excitation has strong intensity by using core and other virtual orbital in the open shell, e.g. core and LUMO+1.

The robust regression example is `HCN_CHP-MRSF_v2.inp`, which uses `swapmo=2,7,8,21`. It obtains two degenerate 1s --> pi* transitions that were missing in the original calculation. Current reference values for `HCN_CHP-MRSF_v2.inp` should be used for automated `openqp --run_tests` validation.

* Numerical note: near-degenerate pi orbitals in linear HCN
The original delta-CHP-MRSF input with only `swapmo=2,7` is intentionally not kept as a `.inp` regression test because it is platform sensitive. Linear HCN has degenerate or nearly degenerate pi orbitals. During the core-hole ROHF SCF with `rstctmo=True`, the maximum-overlap reordering can choose different rotated pi-orbital combinations on different BLAS/LAPACK/compiler platforms. This changes the core-hole reference orbitals slightly before the TD calculation. For example, macOS gives the stored-reference-like core-hole energy near -82.765684506 Hartree, while a Linux/chc2 build can converge to about -82.76568725 Hartree and shifts the TD roots by about 0.0004 Hartree (roughly 0.01 eV). This is a numerical degeneracy/orbital-tracking issue, not a Davidson iteration convergence issue.

For documentation, the platform-sensitive input is preserved below, but the file is named `HCN_CHP-MRSF.inp.example` instead of `HCN_CHP-MRSF.inp` so that `openqp --run_tests` does not discover and run it automatically. Its old reference data is also kept only as `HCN_CHP-MRSF.json.example`, not as an active regression reference:

```ini
[input]
system=
   6     0.000000000      0.000000000     -0.505994000
   7     0.000000000      0.000000000      0.657328000
   1     0.000000000      0.000000000     -1.565330000
basis=6-31g*
functional=bhhlyp
method=tdhf
runtype=energy

[guess]
type=huckel
save_mol=False
continue_geom=False
swapmo=2,7

[scf]
type=rohf
multiplicity=3
save_molden=True
rstctmo=True
conv=1e-9

[dftgrid]
rad_npts=96
ang_npts=302

[tdhf]
type=mrsf
multiplicity=1
nstate=25
ixcore=7
conv=1e-8
zvconv=1e-8
maxit=100
```

