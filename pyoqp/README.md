# Python Wrapper for Open Quantum Platform

2024 Aug 5

## Features

- Native geometry optimization (`[optimize] lib=oqp`) for minima, TS, MECI,
  MECP, TCI, NEB, IRC, and MEP paths, with SciPy and geomeTRIC available as
  optional backends
- Native Fortran initial guesses: `hcore`, `huckel`, `modhuckel`, `minao`, and `sap` (no external quantum-chemistry package required at runtime)
- Energy, gradient, state-overlap, and hop-driver interfaces for
  nonadiabatic molecular dynamics
- Native FSSH NAMD (`runtype=namd`), SOC-NAMD, and ESPF/OpenMM QM/MM NAMD

## Prerequisite

- python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.10.0
- dftd4 >= 3.5.0   Optional, check `setup.py`

## Installation

PyOQP is installed together with the OpenQP native library by the top-level
pip install (see the main [README](../README.md) for details and build options):

      git clone https://github.com/Open-Quantum-Platform/openqp.git
      cd openqp
      pip install .

No environment variables are required afterwards: the installed package
locates its own native library and data files, so do not set `OPENQP_ROOT`.

Only the manual cmake/ninja development flow ("Detailed Compile" in the main
README, where the native library stays in the source tree and PyOQP is
installed from this directory) still needs `OPENQP_ROOT` pointing at that tree.

## Usage

- optional environment variable

      export OMP_NUM_THREADS=4

- test pyoqp

      openqp --run-tests other: Run tests from the 'other' folder in examples
      openqp --run-tests all: Run all tests from all folders in examples
      openqp --run-tests path_to_folder: Run tests from a specific folder

Can be used to run all tests in a specific folder.

Results, including log files and `test_report.txt`, will be stored in the current path in the `oqp_test_tmp_{date}_{time}` folder.

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
    istate=0
    jstate=1
    energy_shift=1e-6
    energy_gap=1e-5
    meci_search=penalty
    pen_sigma=1.0
    pen_alpha=0.0
    pen_incre=1.0
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
    trivial=True
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
    
- energy_shift // convergence threshold for electronic energy changes

      1e-6 (default)
    
- energy_gap // convergence threshold for the energy gap changes

      1e-5 (default)
    
- meci-search // choose the algorithm for conical intersection optimization

      penalty    use the modified pentaly method (default)
      ubp        use the update branching plane method
      hybrid     use the penalty function then swith to ubp after energy gap is lower than the threshold
    
- pen_sigma // set the sigma in the penalty function

      1.0 (defaut)
    
- pen_alpha // set the alpha in the penalty function

      0.0 (default)
    
- pen_incre // set the incremental factor in the penalty function

      1.0 (default)
    
- init_scf // do initial SCF iteration during geometry optimization

      True       do initial SCF iterations in every optimization step
      False      do not do initial SCF iterations after the first optimization step

### [geometric]

geometric section controls the geomeTRIC optimizer backend when [optimize]lib=geometric

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
