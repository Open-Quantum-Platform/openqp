# Python Wrapper for Open Quantum Platform

2024 Aug 5

## Features

- Native geometry optimization for minima, TS, MECI, MECP, TCI, NEB, IRC,
  and MEP paths. Concise `.oqp` workflows select it automatically; traditional
  sectioned `.inp` files may still select SciPy or the optional geomeTRIC backend
- Native Fortran initial guesses: `hcore`, `huckel`, `modhuckel`, `minao`, and `sap` (no external quantum-chemistry package required at runtime)
- Energy, gradient, state-overlap, and hop-driver interfaces for
  nonadiabatic molecular dynamics
- Native FSSH NAMD (`runtype=namd`), SOC-NAMD, and ESPF/OpenMM QM/MM NAMD

## Prerequisite

- python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.10.0
- dftd4 >= 3.5.0   Optional for source-level Python development

## Installation

PyOQP is installed together with the OpenQP native library by the top-level
pip install (see the main [README](../README.md) for details and build options):

      git clone https://github.com/Open-Quantum-Platform/openqp.git
      cd openqp
      pip install .

The default installation contains the complete native optimizer.  Install the
legacy geomeTRIC compatibility backend only for workflows that explicitly use
`[optimize] lib=geometric` (notably its constrained optimizer):

      pip install 'OpenQP[geometric]'

No environment variables are required afterwards: the installed package
locates its own native library and data files, so do not set `OPENQP_ROOT`.

The manual cmake/ninja development flow ("Detailed Compile" in the main
README) is also detected when PyOQP is imported directly from the source
checkout with the native library installed into that tree. Keep `OPENQP_ROOT`
only as a compatibility fallback for a custom root build; the retired
standalone `pyoqp/setup.py` entry point is not an installation path.

## Usage

- optional environment variable

      export OMP_NUM_THREADS=4

- test pyoqp

      openqp --run-tests other                 # examples/other, automatic input selection
      openqp --run-tests all                   # all examples, automatic input selection
      openqp --run-tests path_to_folder        # a specific folder
      openqp --run-tests all --input-format inp
      openqp --run-tests all --input-format oqp
      openqp --run-tests all --input-format both

The default `auto` selection prefers concise `.oqp`, retains legacy-only `.inp`
inputs, and keeps a representative legacy compatibility set. Use `inp`, `oqp`,
or `both` to choose the syntax within the requested test scope. The historical
`all` scope keeps its slow/non-self-contained exclusions; an explicit directory
selects every matching input below that directory. The historical `--run_tests`
spelling remains supported.

Results, including per-input output folders, log files, and `test_report.txt`,
will be stored in the current path in the `oqp_test_tmp_{date}_{time}` folder.

- run pyoqp with command

      openqp input

- import oqp_runner in Python

      from oqp.pyoqp import Runner
      oqp_runner = Runner(project=project_name,
                 input_file=input_name,
                 log=log_name)
      oqp_runner.run()

## Input manual

 Full input sections and keywords

<pre>
[input]
    charge=0
    basis=''
    functional=''
    method=hf
    runtype=energy
    system=''
    d4=False

[guess]
    type=huckel
    file=''
    save_mol=False
    continue_geom=False

[geometric]
    coordsys=tric
    trust=0.1
    tmax=0.3
    convergence_set=GAU
    prefix=geometric
    hessian=never
    irc_direction=forward
    constraints_file=''
    enforce=0.0
    conmethod=0

[scf]
    type=rhf
    maxit=30
    multiplicity=1
    conv=1.0e-6
    incremental=True
    init_scf=no
    init_it=0
    save_molden=True

[dftgrid]
    rad_type=mhl
    rad_npts=96
    ang_npts=302
    partfun=ssf
    pruned=''

[tdhf]:
    type='rpa'
    maxit=50
    multiplicity=1
    conv=1.0e-6
    nstate=1
    zvconv=1.0e-10
    nvdav=50

[properties]:
    scf_prop='el_mom,mulliken'
    td_prop=False
    grad=0
    nac=''
    soc=''
    export=False

[optimize]
    lib=oqp
    optimizer=bfgs
    step_size=0.1
    step_tol=1e-2
    maxit=30
    mep_maxit=10
    rmsd_grad=1e-4
    rmsd_step=1e-3
    max_grad=3e-4
    max_step=2e-3
    istate=1
    jstate=2
    kstate=3
    states=
    energy_shift=1e-6
    energy_gap=1e-5
    meci_search=auto
    mecp_search=auglag
    gap_sigma=10.0
    pen_sigma=1.0
    pen_alpha=0.0
    pen_incre=1.0
    pen_delta=0.025
    pen_jump=10,10,25,25,100,100,1000,1000,3000
    gap_weight=1.0
    init_scf=False

[hess]
    type=numerical
    state=0
    dx=0.01
    nproc=1
    read=False
    restart=False

[md]
    nstep=100
    dt=0.5
    active=1
    substep=200
    decoherence=edc
    thrshe=1.0e9
    tdc=fd
    trivial=False
    init_temp=300.0
    velocity=maxwell
    seed=1
    restart=False
    soc=False
    soc_basis=adiabatic
    soc_du_dt_corr=False
    soc_tdc_grad_corr=False
    grad_wthr=0.001
    init_state=''

</pre>

Details for keywords

### [input] 

input section handle the basic information of molecular system

- charge // set the total charge

      default is 0

- basis // set the basis set

      check options in opq/share/basis_sets
      no default value

- functional // set the functional 
      
      check options in libxc list 
      default is empty, which will do Hartree-Fork

- method // choose the type of HF/DFT calculation
      
      hf         time-independent calculations, HF, DFT (default)
      tdhf       timd-dependent calculation, TDDFT, MRSF-TDDFT

- runtype // choose the type of oqp calculation
       
      energy     single-point energy (default)
      grad       single-point energy and gradients
      hess       frequency calculation (numerical only)
      nac        nonadiabatic coupling
      soc        spin-orbit coupling
      namd       fewest-switches nonadiabatic molecular dynamics
      optimize   local minimum geometry optimization
      meci       minimum energy conical intersection optimization
      mep        minimum energy path calculation
      ts         transition state optimization
      irc        intrinsic reaction coordinate
      neb        nudged elastic band calculation
 
- system // specify molecular structure or xyz file

      system=filename.xyz  open a xyz file

      or input coordinate in the next line as the following
      system=
           O  -0.0000000000   0.0000000000  -0.0410615540
           H  -0.5331943294   0.5331943294  -0.6144692230
           H   0.5331943294  -0.5331943294  -0.6144692230

      note there should be at least one space in front of the element symbol

- d4 // apply DFTD4 dispersion correction

      False      do not compute DFTD4 corrections according to functional (default)
      True       compute DFTD4 corrections for energy and gradients. some functional might not be supported

### [guess]

guess section handle the guess orbitals 

- type // choose the type of guess orbital

      sap        superposition of atomic potentials (default, native Fortran)
      minao      projected atomic minimal-basis densities (native Fortran)
      huckel     extended Huckel guess (native Fortran)
      modhuckel  modified (weighted Wolfsberg-Helmholz) Huckel (native Fortran)
      hcore      bare core-Hamiltonian guess
      model      read orbital from molden
      json       load data from json
      auto       load json if the requested file exists; otherwise use huckel

- file // set the guess orbital or data file

      filename   name or absolute path to molden or json file

- save_mol // save complete data to a json file

For overlap-aligned excited-state calculations, the saved JSON also contains
a public `state_tracking` record.  `order` maps each saved current-root index
to its previous-step index (zero-based), while `lineage` carries the initial
physical-state identity through root exchanges.  `phase_step` is the sign
applied to each raw current response vector and `phase_initial` is that raw
vector's correction to the transported initial gauge.  External dynamics
drivers should consume this joint state gauge instead of independently fitting
the sign of every pairwise NAC vector.  `matched_overlap` and `margin` expose
weak or ambiguous correspondences.
`raw_order` preserves the pre-alignment solver-root map; `output_reordered` is
true only for internal numerical-NAC displacement workers, whose response
vectors are restored to the central structure's physical-root order before
the +dx/-dx finite difference.

      True       save complete calculation data to json file
      False      do not save data (default)
    
- continue_geom // choose structure for calculations

      True       use the structure saved in the json file
      False      use the input structure (default)

### [scf]

scf section handle the time-independent calculations

- type // choose the type of wavefunction

      rhf        restricted Hartree-Fork/Kohn-Sham (default)
      uhf        unrestricted Hartree-Fork/Kohn-Sham
      rohf       restricted Hartree-Fork/Kohn-Sham

- maxit // set the maximum number of SCF iterations

      30 (default)

- multiplicity // set the multiplicity of the reference state

      1 (default)
    
- conv // set energy convergence

      1.0e-6 (default)

- incremental // use incremental Fork method

      True       use the incremental method (default)
      False      do not use the incremental method

- init_scf // do initial SCF iteration to help convergence

      rhf        do initial RHF type HF calculation, regardless the functional
      uhf        do initial UHF type HF calculation, regardless the functional
      rohf       do initial ROHF type HF calculation, regardless the functional
      rks        do initial RKS type DFT calculation according to the functional
      uks        do initial UKS type DFT calculation according to the functional
      roks       do initial ROKS type DFT calculation according to the functional
      no         do not do initial SCF iteration (default)

- init_it // set the maximum number of initial SCF iteration

      0 (default)

- save_molden // save orbitals

      True       save orbitals to a molden file
      False      do not save orbitals

### [dftgrid]

dftgrid section handle the accuracy of the DFT calculations

- rad_type // choose the radial point sampling methods for electronic integrals

      mhl (default)

- rad_npts // set the number of radial point for electronic integrals

      96 (default)
    
- ang_npts // set the number of angular point for electronic integrals

      302 (default)
    
- partfun // choose the partition function for electronic integrals

      ssf (default)
    
- pruned // choose orbital prune method for electronic integrals

      empty      no prune (default)
      ao         prune atomic orbital with SG1 scheme

### [tdhf]

tdhf section handle the time-dependent calculations

- type // choose the type of time-dependent wavefunction

      rpa        use random phase approximation (default)
      tda        use Tam-Dancoff  approximation 
      sf         use spin-flip
      mrsf       use mixed-reference spin-flip
    
- maxit // set the maximum number CI iterations

      50 (default)
    
- multiplicity // set the multiplicity of the response state

      1 (default)
    
- conv // set the energy convergence for the response state

      1.0e-6 (default)
    
- nstate // set the number of response state

      1 (default)
    
- zvconv // set the convergence of Z-vector calculation

      1.0e-6 (default)
    
- nvdav // set the dimension of the Davidson subspace

      50 (default)

### [properties]

properties section handel the property calculation

- scf_prop // compute reference state properties

      el_mom     electronic momentum
      mulliken   mulliken charges
      
      scf_prop=el_mom,mulliken will compute both (default)
    
- td_prop // compute response state properties

      True       compute the electronic properties for response state (not available yet)
      False      do not compute the electronic properties for response state
    
- grad // compute the gradient for given states

      0           the reference state (default)
      1,2,3       the first, second and third response states
      
      currently, time-dependent calculations ([input]method=tdhf) do not compute
      the gradients for reference state (grad=0)

    
- nac // compute non-adiabatic coupling for two given states

      nac         compute derivative-coupling vectors
      nacme       compute nonadiabatic coupling matrix elements
    
- soc // compute the spin-orbit coupling for two given states

      use [input]runtype=soc for an SOC calculation
    
- export // save the computed data to text files

      True        save the energies, gradients
      False       do not save data

### [md]

md section controls fewest-switches nonadiabatic molecular dynamics
(`[input] runtype=namd`).

- soc // enable spin-orbit-coupled dynamics

      False       internal-conversion NAMD without SOC (default)
      True        SOC-NAMD / ISC dynamics

- soc_basis // choose the SOC propagation and force representation

      adiabatic   SHARC-like spin-adiabatic propagation with weighted MCH gradients (default)
      mch         MCH-basis SOC propagation with exact active-root MCH gradients; recommended for SOC-QM/MM production

- soc_du_dt_corr // optional spin-adiabatic diagnostic correction

      True        add finite-difference dU/dt force correction
      False       no dU/dt correction (default)

- soc_tdc_grad_corr // optional spin-adiabatic diagnostic correction

      True        add MCH TDC-projected derivative-coupling force correction
      False       no TDC-projected force correction (default)

For QM/MM NAMD, set `[input] qmmm_flag=True`; `runtype=namd` remains on the
normal OpenQP Runner path.  The legacy OpenMM ground-state MD path is selected
only by `runtype=md` with `qmmm_flag=True`.

### [optimize]

optimize section handle the geometry optimization

This is the traditional sectioned `.inp` interface.  New concise `.oqp` files
do not expose `lib`; all geometry drivers use the native OpenQP optimizer
automatically.

- lib // choose the optimization library

      oqp         use the native optimizer (default). Supports optimize, ts, meci, mecp, tci, neb, irc, and mep
      scipy       use scipy.optimize library for optimize, meci, mecp, and mep
      geometric   use geomeTRIC for optimize, meci, mecp, ts, irc, and neb

- optimizer // choose the scipy optimizer

      bfgs        use BFGS method (default)
      cg          use conjugated gradient
      l-bfgs-b    use L-BFGS-b method
      newton-cg   use Netown conjugated gradient
    
- step_size // set the radius of the constraining hypersphere from the starting structure

      0.1         the largest distance between the mass-weighted coordinates (default)
    
- step_tol // set the threshold for the radius on a hypersphere from the starting structure 

      1e-2        the smallest distance between the mass-weighted coordinates (default) 
    
- maxit // the maximum number of geometry optimization iterations

      30 (defult)
    
- mep_maxit // the maximum number of mep steps

      10 (default)
    
- rmsd_grad // convergence threshold for rmsd of gradients

      1e-4 (default)
    
- rmsd_step // convergence threshold for rmsd of structure changes

      1e-3 (default)
    
- max_grad // convergence threshold for max of gradients

      3e-4 (default)
    
- max_step // convergence threshold for max of structure changes

      2e-3 (default)
    
- istate // choose the state for single state optimization

      1 (default)
    
- jstate // choose the second state for conical intersection optimization

      2 (default)

- states // ordered consecutive response roots for BaekA multistate MECI

      empty (default); use 1,2 or 1,2,3,... with meci_search=auto or baeka
    
- energy_shift // convergence threshold for electronic energy changes

      1e-6 (default)
    
- energy_gap // convergence threshold for the energy gap changes

      1e-5 (default)
    
- meci-search // choose the algorithm for conical intersection optimization

      auto       default: two-state searches use auglag, with BaekA held back
                 as a rescue on the recovery budget if auglag does not meet
                 the criteria; multistate searches select BaekA directly.
                 Backends other than oqp resolve auto to auglag
      penalty    use the modified penalty method
      ubp        use the update branching plane method
      auglag     use the branching-plane projection with a least-squares
                 Lagrange multiplier; the reported objective is an augmented
                 Lagrangian value rather than the ratio used by ubp
      hybrid     use the penalty function then swith to ubp after energy gap is lower than the threshold
      baeka      use the additive Baek adaptive penalty for two or more states

- mecp-search // choose the algorithm for spin-crossing (MECP) optimization

      auto       (default) sqp when lib=oqp, auglag on the backends that
                 bring their own optimizer
      sqp        sequential quadratic programming: solves the KKT equations
                 of the constrained problem for the step and the multiplier
                 together, so the multiplier is a result rather than a formula
                 and there is no penalty parameter, i.e. gap_sigma is unused.
                 Works in delocalized internal coordinates with the native
                 model Hessian. Requires lib=oqp because it replaces the outer
                 optimizer with its own trust-region step control. coordsys=cart
                 is honoured; the other settings, including tric, select DLC,
                 because a dense KKT system needs a non-redundant basis.
                 Converges tighter and in fewer steps than auglag on the cases
                 tested so far, which is why auto selects it natively.
      auglag     least-squares multiplier plus gradient projection; with
                 gap_sigma=1 the converged form is the plain Bearpark gradient
                 projection. The gap term and the projected mean gradient are
                 orthogonal, so both must vanish separately and the stationary
                 point is a true crossing. Selected by auto on the scipy and
                 geometric backends, which supply their own optimizer
      penalty    Levine-Martinez smooth penalty
      quad       legacy fixed-weight quadratic penalty.  Its stationary point
                 balances the mean gradient against the gap term, leaving a
                 residual gap of order 1/gap_weight, so it generally cannot
                 satisfy energy_gap.  Kept only to reproduce earlier runs

- gap_sigma // strength of the gap term in the auglag objectives

      10.0 (default), must be positive; 1.0 reproduces the plain Bearpark
      projection, larger values reach the seam faster

- pen_sigma // set the sigma in the penalty function

      1.0 (defaut)
    
- pen_alpha // set the alpha in the penalty function

      0.0 (default)
    
- pen_incre // set the incremental factor in the penalty function

      1.0 (default); multiplicative legacy control used by runtype=tci,
      not by meci_search=baeka

- pen_delta // BaekA small additive delta-beta update

      0.025 (default)

- pen_jump // BaekA large additive beta schedule

      10,10,25,25,100,100,1000,1000,3000 (default)

- gap_weight // crossing-search penalty scaling

      1.0 (default); fixed at 1.0 for meci_search=baeka
    
- init_scf // do initial SCF iteration during geometry optimization

      True       do initial SCF iterations in every optimization step
      False      do not do initial SCF iterations after the first optimization step

### [geometric]

Legacy optional section controlling geomeTRIC when a traditional `.inp` file
sets `[optimize] lib=geometric`.  It is not part of the concise `.oqp` grammar;
install it with `pip install 'OpenQP[geometric]'` when needed.

- coordsys // coordinate system passed to geomeTRIC

      tric       translation-rotation internal coordinates (default)

- trust // initial trust radius in Angstrom

      0.1 (default)

- tmax // maximum trust radius in Angstrom

      0.3 (default)

- convergence_set // geomeTRIC convergence preset

      GAU (default)

- prefix // prefix for geomeTRIC output files

      geometric (default)

- hessian // initial Hessian option passed to geomeTRIC

      never (default for minima and constrained optimization)
      first (recommended/automatically used for ts and irc if hessian=never)

- irc_direction // IRC branch direction for runtype=irc

      forward (default)
      backward

- constraints_file // geomeTRIC constraints file for constrained optimization

      relative paths are resolved relative to the OpenQP input file

- enforce // geomeTRIC constraint enforcement value

      0.0 (default)

- conmethod // geomeTRIC constraint method

      0 (default)

### [hess]
hess section handle hessian and frequence calculations

- type // set the type of hessian

      numerical   compute hessian numerically by evaluating the gradient with a small displacemet(default)

- state // set the state for frequency calculation
    
      0 (defaut)

- dx // set a small displacement for numerical hessian

      0.01 (default) move each coordiante forward and backward with the displacement. The unit is in Bohr
    
- nproc // set the number of subprocess for numerical hessian

      1 (default)  the total number of running CPU is nproc * OMP_NUM_THREADS       
    
- read // read computed hessian data

      False        read the .hess.json file to retreive the computed frequency data
    
- restart // restart hessian calculation

      False        read the strach data in _num_hess folder to continue the calculation
