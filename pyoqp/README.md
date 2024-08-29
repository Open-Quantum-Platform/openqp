# Python Wrapper for Open Quantum Platform

2024 Aug 5

## Features

- Geometry optimization for arbitrary state local minimum, MECI, MEP with Scipy library
- Geometry optimization for arbitrary state local minimum, MECI, TS with DL-FIND library
- Energy gradient calculations interface for nonadiabatic molecular dynamics

## Prerequisite

- python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.10.0
- libdlfind >= 0.0.3
- dftd4 >= 3.5.0   Optional, check `setup.py`

## Installation

- download the source files

      git https://qchemlab.knu.ac.kr/open-quantum-package/modules.git

- compile OQP

      following instructions in /modules/README.md

- install PyOQP

      cd ./modules/pyoqp/pyoqp
      pip install .

## Usage

- required environment variable

      export OPENQP_ROOT=/path/to/oqp
      export OMP_NUM_THREADS=4
      export LD_LIBRARY_PATH=/other/dependent/lib:$LD_LIBRARY_PATH

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
    lib=scipy
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

[dlfind]
    printl=2
    icoord=0
    iopt=3
    ims=0

[hess]
    type=numerical
    state=0
    dx=0.01
    nproc=1
    read=False
    restart=False

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
      nac        non-adiabatic coupling (not available yet)
      soc        spin-orbit coupling (not available yet)
      optimize   local minimum geometry optimization
      meci       minimum energy conical intersection optimization
      mep        minimum energy path calculation
      ts         transition state optimization
      neb        nudge elasted band calculation (not avaialble yet)
 
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

      huckel     huckel guess (default)
      hcore      hcore guess
      model      read orbital from molden
      json       load data from json

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

      not available yet
    
- soc // compute the spin-orbit coupling for two given states

      not available yet
    
- export // save the computed data to text files

      True        save the energies, gradients
      False       do not save data

### [optimize]

optimize section handle the geometry optimization

- lib // choose the optimization library

      scipy       use scipy.optimize library (default)
      dlfind      use DL-FIND library

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

### [dlfind]

dlfind section handle the DL-FIND library for geometry optimization

- printl // set the DL-FIND printing level

      2 (default)
    
- icoord // choose the coordinate

      0          Cartesian
      1          hybrid delocalized internal coordinates, primitive internal coordinate scheme
      2          hybrid delocalized internal coordinates, total connection scheme
      3          delocalized internal coordinates, primitive internal coordinate scheme
      4          delocalized internal coordinates, total connection scheme
      10–14      Lagrange–Newton conical intersection search, with 2nd digit referring to above options

- iopt // choose the optimization job

      0          steepest descent
      1          Polak-Ribiere conjugate gradient w/ automatic restart
      2          Polak-Ribiere conjugate gradient w/ restart every 10 steps
      3          L-BGFS (default)
      9          P-RFO, for transition state searches, require [input]runtype=ts
    
- ims // set the multistate gradient calculations

      0          single-state calculation (default)
      1          conical intersection optimization with penalty function algorithm, require [input]runtype=meci
      2          conical intersection optimization with gradient projection algorithm, require [input]runtype=meci
      3          conical intersection optimization with Lagrange–Newton algorithm, require [input]runtype=meci

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
