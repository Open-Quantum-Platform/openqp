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

One first need to swap MO (core and HOMO) by "swapmo=2,7" in the "HCN_CHP-MRSF.inp" example. Then one need to restrict the MO ordering during SCF by "rstctmo=True" keyword.  Once the SCF is converged, one can get excitations from core (HOMO) by "ixcore=7" keyword.

In the example, Carbon K-edge is computed from state 2, which is a 1s --> pi* transition. The ground state obtained from CHP-MRSF is usually poorly described, therefore, it is better to use the ground state from KS DFT (or MRSF-TDDFT) calculation.

Current reference values for `HCN_CHP-MRSF.inp` are: core-hole reference state 0 = -82.7656845061 Hartree; bright state 2 = -83.3333968722 Hartree; transition 1 -> 2 = 293.81 eV with oscillator strength 0.0530. The DFT ground-state reference from `HCN_DFT.inp` is -93.3641769363 Hartree. If an absolute delta-CHP-MRSF energy is formed as bright CHP state 2 minus KS DFT, the value is 10.0307800641 Hartree (272.95 eV). Do not use the core-hole reference state 0 energy as the bright XAS state energy.

* Important note for delta-CHP-MRSF
In current implementation, we don't get core (HOMO after using swapmo and rstctmo) -> LUMO excitation.  
Therefore, one must check if core->LUMO excitation has strong intensity by using core and other virtual orbital in the open shell, e.g. core and LUMO+1.

For example, in the HCN_CHP-MRSF_v2.inp, it obtained two degenerate 1s --> pi* transitions which were missing in the original calculation.

