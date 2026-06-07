## Open Quantum Platform: OpenQP

Open Quantum Platform ([OpenQP](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)) is a quantum chemical platform featuring cutting-edge capabilities like [Mixed-Reference Spin-Flip (MRSF)-TDDFT](https://doi.org/10.1021/acs.jpclett.3c02296) with an emphasis on open-source ecosystem.

### Current Functionality

| Area | What OpenQP can do now | Notes |
| --- | --- | --- |
| Electronic structure | <small>HF, DFT, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, UMRSF-TDDFT, and MRSF-EKT ground- and excited-state calculations</small> | <small>MRSF-TDDFT includes DTCAM-series exchange-correlation functionals; UMRSF-TDDFT provides MRSF excitation energies from a UHF reference (energy-only); MRSF-EKT supports IP/EA analysis.</small> |
| Derivative properties | <small>Energies, analytic gradients, numerical Hessians, and native HF/DFT analytic Hessians</small> | <small>Native CPHF/CPKS response, integral-derivative kernels, and Hessian assembly cover closed-shell and open-shell (UHF/ROHF) references, ECPs, and range-separated (CAM/LRC) functionals.</small> |
| Vibrational analysis | <small>Frequencies, normal-mode eigenvector printout, thermochemistry, and native IR/Raman intensity assembly</small> | <small>IR/Raman intensities use native OpenQP dipole, CPHF polarizability, and vibrational-intensity kernels for closed- and open-shell references.</small> |
| Molecular symmetry | <small>Point-group detection and symmetry labels for MOs, excited states, and vibrational modes</small> | <small>Petite-list symmetry reductions accelerate integral, XC, gradient, and response evaluations.</small> |
| DFT integration grids | <small>Lebedev grids plus SG-1, SG-2, and SG-3 pruned grids with per-element DE2 radial quadrature</small> | <small>OpenMP-parallel XC integration with optimized numerical kernels.</small> |
| Nonadiabatic dynamics data | <small>MRSF-TDDFT nonadiabatic couplings and NACME-oriented workflows</small> | <small>NAC support uses the TLF-based MRSF-TDDFT machinery.</small> |
| Geometry/path optimization | <small>Minima, transition states, MECI/MECP, constrained optimization, IRC, and NEB-style workflows</small> | <small>SciPy and geomeTRIC backends are available depending on the requested job.</small> |
| Initial guesses and SCF stability | <small>Native initial guesses (`hcore`, `huckel`, `modhuckel`, `minao`, `sap`), JSON restart/auto guesses, and OpenTrustRegion SCF stabilization</small> | <small>Native guess paths avoid external runtime dependencies and are suitable for MPI execution; Davidson auto-restart and MINRES/AUTO Z-vector fallbacks harden excited-state and response convergence.</small> |
| Integrations | <small>LibXC, basis_set_exchange, libecpint, DFT-D4, PyRAI2MD, and OpenqpView</small> | <small>External viewers and workflow integrations are available around the native OpenQP electronic-structure engines.</small> |
| Performance and deployment | <small>OpenMP/MPI execution, BLAS/LAPACK optimization, source builds, pip installs, and Docker images</small> | <small>MPI requires an MPI implementation such as OpenMPI.</small> |

### Upcoming Features
- **Efficient electrostatic embedding QM/MM** by [ESPF QM/MM](https://doi.org/10.1063/5.0133646)
- **Spin-Orbit Coupling** by [**Relativistic** MRSF-TDDFT](https://doi.org/10.1021/acs.jctc.2c01036)

### Quickstart

- **pip install openqp** or
- **Ready to Use Docker Image** of [openqp/openqp](https://github.com/Open-Quantum-Platform/openqp/wiki/OpenQP_Docker_Image) or
- **Building from Source Files Using the Instructions Below.**

#### Requirements

- **GCC, G++, Gfortran**: Version >= 8
- **CMake**: Version >= 3.25
- **cffi**: Perform pip install cffi
- **ninja** (optional)
- **MPI Library**: OpenMPI For MPI Support. Consult detailed documentation for other MPI libraries

#### Download the Source Files

```bash
git clone https://github.com/Open-Quantum-Platform/openqp.git
```

#### Pip install
```bash
cd openqp
pip install .
```
This is the recommended source install path. It builds and installs the OpenQP Python package and native library together, so setting `OPENQP_ROOT` is not required for normal `openqp` command-line use after installation. Native initial guesses such as `huckel`, `modhuckel`, `minao`, and `sap` work without external guess-generation dependencies.

or 
#### Detailed Compile

##### macOS with Homebrew GCC and native BLAS

On macOS, use the Homebrew GCC toolchain for C/C++/Fortran and the native Accelerate BLAS/LAPACK provider. Replace `gcc-15`, `g++-15`, and `gfortran-15` with the installed Homebrew GCC version if needed.

```bash
cd openqp
cmake -B build -G Ninja \
  -DUSE_LIBINT=OFF \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON \
  -DLINALG_LIB=auto \
  -DLINALG_LIB_INT64=OFF
ninja -C build install
cd pyoqp
pip install .
```

##### OpenMP Support

```bash
cd openqp
cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON
ninja -C build install
cd pyoqp
pip install .
```

##### OpenMP and MPI Support

```bash
cd openqp
cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON -DENABLE_MPI=ON
ninja -C build install
cd pyoqp
pip install .
```

##### OpenMP and MPI Support using make

```bash
cd openqp
cmake -B build -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON -DENABLE_MPI=ON
make -C build install
cd pyoqp
pip install .
```

- Use `-DUSE_LIBINT=ON` to replace the default ERI based on Rys Quadrature with `libint`.
- The Linux/manual examples above use ILP64 BLAS/LAPACK (`-DLINALG_LIB_INT64=ON`, the default). The macOS example uses native Accelerate BLAS/LAPACK with `-DLINALG_LIB_INT64=OFF`.

#### Environmental Settings

```bash
export OPENQP_ROOT=/path/to/openqp                           # Path to the Root of openqp
export OMP_NUM_THREADS=4                                     # The number of cores to be used for OpenMP runs
export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
```

**Special Environmental Settings for MKL Math Library:**

```bash
export MKL_INTERFACE_LAYER="@_MKL_INTERFACE_LAYER@"
export MKL_THREADING_LAYER=SEQUENTIAL
```

#### Test

```bash
openqp --run_tests all     # Run all tests from all folders in examples
```

#### Run

For OpenMP or sequential run:

```bash
openqp any_example_file.inp
```

For OpenMP and MPI run:

```bash
mpirun -np number_of_mpi openqp any_example_file.inp
```

### Detailed Documentation

For more in-depth information, visit:
- [OpenQP Documentation](https://github.com/Open-Quantum-Platform/openqp/wiki)

### Input Generator
Easily create input files for OpenQP using our [Web-based Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/).

### OpenqpView
Inspect OpenQP calculation outputs directly in the browser with [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/). Recent OpenqpView development added a GitHub Pages deployment, full-periodic-table molecule rendering, local file/drop/paste loading, WebGL auto-spin controls, and support for OpenQP log, JSON, Molden, cube, and XYZ data. Files and pasted text are processed locally in the browser and are not uploaded to a server.

### Citing OpenQP
If you use OpenQP in your research, please cite the following papers:

- **Mironov V, Komarov K, Li J, Gerasimov I, Mazaheri M, Park W, Lashkaripour A, Oh M, Nakata H, Ishimura K, Huix-Rotllant M, Lee S, and Choi CH.** "OpenQP: A Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-source Ecosystem" [Journal of Chemical Theory and Computation, 2024](https://doi.org/10.1021/acs.jctc.4c01117)
- **Park W, Komarov K, Lee S, and Choi CH.** "Mixed-Reference Spin-Flip Time-Dependent Density Functional Theory: Multireference Advantages with the Practicality of Linear Response Theory." [The Journal of Physical Chemistry Letters. 2023 Sep 28;14(39):8896-908.](https://doi.org/10.1021/acs.jpclett.3c02296)
- **Lee S, Filatov M, Lee S, and Choi CH.** "Eliminating Spin-Contamination of Spin-Flip Time-Dependent Density Functional Theory Within Linear Response Formalism by the Use of Zeroth-Order Mixed-Reference (MR) Reduced Density Matrix." [The Journal of Chemical Physics, vol. 149, no. 10, 2018.](https://doi.org/10.1063/1.5044202)
- **Lee S, Kim EE, Nakata H, Lee S, and Choi CH.** "Efficient Implementations of Analytic Energy Gradient for Mixed-Reference Spin-Flip Time-Dependent Density Functional Theory (MRSF-TDDFT)." [The Journal of Chemical Physics, vol. 150, no. 18, 2019.](https://doi.org/10.1063/1.5086895)

### Contributors

- **Cheol Ho Choi**, Kyungpook National University, South Korea, [cheolho.choi@gmail.com](mailto:cheolho.choi@gmail.com), [https://www.openqp.org](https://www.openqp.org)
- **Seunghoon Lee**, Seoul National University, South Korea, [seunghoonlee89@gmail.com](mailto:seunghoonlee89@gmail.com)
- **Vladimir Mironov**, [vladimir.a.mironov@gmail.com](mailto:vladimir.a.mironov@gmail.com)
- **Konstantin Komarov**, [constlike@gmail.com](mailto:constlike@gmail.com)
- **Jingbai Li**, Hoffmann Institute of Advanced Materials, China, [lijingbai2009@gmail.com](mailto:lijingbai2009@gmail.com)
- **Igor Gerasimov**, [i.s.ger@yandex.ru](mailto:i.s.ger@yandex.ru)
- **Hiroya Nakata**, Fukui Institute for Fundamental Chemistry, Japan, [nakata.hiro07@gmail.com](mailto:nakata.hiro07@gmail.com)
- **Mohsen Mazaherifar**, Kyungpook National University, South Korea, [moh.mazaheri@gmail.com](mailto:moh.mazaheri@gmail.com)
### Legal Notice

See the separate LICENSE file.
