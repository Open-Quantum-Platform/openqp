# Results and Molecule Data

After `runner.run()`, the `Runner` keeps the active calculation object at
`runner.mol`. This `Molecule` object stores the parsed configuration, native
OpenQP data handle, scalar results, and NumPy-accessible arrays.

```python
runner.run()

mol = runner.mol
print(mol.config["input"]["runtype"])
print(mol.get_scf_energy())
print(mol.get_system())
```

## Runner Summary

`runner.results()` returns a compact dictionary for Python-to-Python
communication:

| Key | Meaning |
| --- | --- |
| `atoms` | Atomic numbers. |
| `system` | Cartesian coordinates in Bohr. |
| `energy` | Energy array or scalar data accumulated by the workflow. |
| `grad` | Gradient data when available. |
| `dcm` | Derivative-coupling matrix data when available. |
| `nac` | Nonadiabatic coupling data when available. |
| `soc` | Spin-orbit coupling data when available. |
| `data` | Raw OpenQP data tags serialized into Python lists. |

For file-style result export, use `runner.mol.get_results()`. It returns a JSON
friendly dictionary with atoms, coordinates, total energy, symmetry metadata,
TDDFT energies, gradients, NAC, SOC, Hessian data, MP2 correlated totals, and
MRSF-EKT records when present.

Wavefunction calculations additionally populate `get_results()["energies"]`
with ordered total energies in Hartree. For FCI, CASCI, and CASSCF this is the
per-root total-energy list; a state-averaged CASSCF calculation still reports
the solved roots even though its scalar objective is a weighted average. For
CASPT2, NEVPT2, and QDPT variants the values are the corrected total energies,
not the perturbative corrections alone, ordered by the requested target roots
(or by ascending eigenvalue for a diagonalized multistate model space).
Corresponding native tags include `OQP::FCI_ENERGIES`, the CASCI/CASSCF energy
tags, and `OQP::CASPT2_ENERGIES`; `mol.energies` is the source used for this
portable result key.

For dedicated SA-CASSCF gradients, `gradient_state=averaged` stores the
weighted-objective derivative in gradient row zero, but the energy array still
contains the individual roots; there is no matching weighted-objective energy
row. An integer `gradient_state=J` stores the analytic Z-vector derivative in
physical CI-root row `J`, and the corresponding reported energy is that same
root. Input validation requires the direct-gradient or optimizer selector to
address the same row.

For `method=mp2`, the exported energy is the correlated MP2 total after the
SCF reference energy and MP2 correlation are combined. The run log also prints
the reference SCF energy plus same-spin and opposite-spin MP2 components.

For `method=ccsd` and `method=ccsd(t)`, the exported energy follows the same
rule: it is the correlated total, not the reference. `ccsd` exports
`E(SCF) + E(CCSD correlation)`, and `ccsd(t)` exports
`E(SCF) + E(CCSD correlation) + E((T) correction)` — so switching between them
changes the value under the same key, and the key alone does not say which
method produced it. The individual pieces are not exported separately; the run
log prints them as a labelled table:

| Log line | Meaning |
| --- | --- |
| `E(reference, SCF)` | The converged HF reference energy. |
| `E(CCSD, correlation)` | The CCSD correlation energy alone. |
| `E(CCSD, total)` | Reference plus CCSD correlation. |
| `E((T), correction)` | The perturbative triples correction (`ccsd(t)` only). |
| `E(CCSD(T), correlation)` | CCSD correlation plus triples (`ccsd(t)` only). |
| `E(CCSD(T), total)` | The exported total for `ccsd(t)`. |

Consumers that need the decomposition rather than the total should record
which method they requested and read these lines. The exported keys may gain a
dedicated coupled-cluster breakdown in a later release.

When the corresponding property is requested via `[properties] scf_prop`,
`get_results()` also includes the following (identically for file-based and
`input_dict`/scripting-API runs):

| Key | Requested by | Meaning |
| --- | --- | --- |
| `dipole` | `el_mom` | Electric dipole vector (a.u.). |
| `mulliken_charges` | `mulliken` | Mulliken atomic partial charges (e). |
| `lowdin_charges` | `lowdin` | Löwdin atomic partial charges (e). |
| `resp_charges` | `resp` | RESP/ESP-fitted atomic charges (e). |
| `nmr_shielding` | `nmr` | Isotropic NMR shielding (ppm), shape `(natom, 5)`: columns are `dia`, `para_uncoupled`, `para_coupled`, `total_uncoupled`, `total_coupled`. |

`nac` is populated with the NACME derivative-coupling matrix for `runtype=nacme`
and is empty otherwise.

## Common Molecule Methods

| Method | Returns |
| --- | --- |
| `get_atoms()` | Atomic numbers as a copied NumPy array. |
| `get_mass()` | Atomic masses as a copied NumPy array. |
| `get_system()` | Cartesian coordinates in Bohr as a copied NumPy array. |
| `get_scf_energy()` | Total SCF energy. |
| `get_scf_energy("all")` | Dictionary of available SCF energy components. |
| `get_grad()` | Gradient in Hartree/Bohr. |
| `get_hess()` | Cartesian Hessian matrix when available. |
| `get_soc()` | SOC eigenvalues in cm-1 when available. |
| `get_data()` | Raw registered OpenQP data tags converted to lists. |
| `get_results()` | JSON-friendly calculation summary. |
| `write_molden(filename)` | Writes Molden orbital data. |

## Dynamic Tag Accessors

OpenQP registers several native data tags and creates matching Python getter
and setter methods. The method names are lower-case versions of the tag names.

| Tag | Getter | Typical content |
| --- | --- | --- |
| `OQP::DM_A` | `get_dm_a()` | Alpha density matrix. |
| `OQP::DM_B` | `get_dm_b()` | Beta density matrix. |
| `OQP::FOCK_A` | `get_fock_a()` | Alpha Fock matrix. |
| `OQP::FOCK_B` | `get_fock_b()` | Beta Fock matrix. |
| `OQP::E_MO_A` | `get_e_mo_a()` | Alpha orbital energies. |
| `OQP::E_MO_B` | `get_e_mo_b()` | Beta orbital energies. |
| `OQP::VEC_MO_A` | `get_vec_mo_a()` | Alpha MO coefficients. |
| `OQP::VEC_MO_B` | `get_vec_mo_b()` | Beta MO coefficients. |
| `OQP::Hcore` | `get_hcore()` | Core Hamiltonian. |
| `OQP::SM` | `get_sm()` | AO overlap matrix. |
| `OQP::TM` | `get_tm()` | Kinetic-energy matrix. |
| `OQP::td_energies` | `get_td_energies()` | TDHF/TDDFT excitation energies. |
| `OQP::td_mrsf_density` | `get_td_mrsf_density()` | MRSF density data. |
| `OQP::td_states_overlap` | `get_td_states_overlap()` | State-overlap matrix. |
| `OQP::soc_eval` | `get_soc_eval()` | SOC eigenvalues. |

Each dynamic tag also has a setter named `set_<tag>()`, for example
`set_dm_a(array)`. Setters mutate the native data store and are intended for
advanced workflows such as custom SCF experiments or internal coupling between
OpenQP modules.

## Units and Shapes

- Coordinates returned by `get_system()` are in Bohr.
- Gradients returned by `get_grad()` are in Hartree/Bohr.
- `get_soc()` reports SOC eigenvalues in cm-1.
- Many raw tag arrays are stored in flattened native order. Convert or reshape
  them using the current basis size, atom count, or state count before analysis.

When possible, prefer `runner.results()` or `mol.get_results()` for automated
post-processing because those summaries are less tied to native storage layout.
