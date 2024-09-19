## Open Quantum Platform: OpenQP
**August 23, 2024**

Open Quantum Platform (OpenQP) is a quantum chemical platform featuring cutting-edge capabilities like [Mixed-Reference Spin-Flip (MRSF)-TDDFT](https://doi.org/10.1021/acs.jpclett.3c02296) with an emphasis on open-source ecosystem.

### Key Features

- **Autonomous Modules of Quantum Chemistry Theories for Easy Interoperability**
- **Ground and Excited State Properties** by [MRSF-TDDFT](https://doi.org/10.1021/acs.jpclett.3c02296)
- **Nonadiabatic Coupling** based on [TLF Technology](https://doi.org/10.1021/acs.jpclett.1c00932) using **MRSF-TDDFT**
- **New Exchange-Correlation Functionals** of [**DTCAM** series](https://doi.org/10.1021/acs.jctc.4c00640) for MRSF-TDDFT
- **Ground State Properties** by HF and DFT theories
- **Geometry Optimization, Transition State Search, and Conical Intersection Search** by SciPy and [DL-Find](https://github.com/digital-chemistry-laboratory/libdlfind)
- [LibXC](https://gitlab.com/libxc/libxc) Integration to support a variety of exchange-correlation functionals
- **Support for [Molden](https://www.theochem.ru.nl/molden/) File Format** for visualization, compatible with many graphic software tools
- [DFT-D4 Dispersion Correction](https://dftd4.readthedocs.io/en/latest/)
- **OpenMP and MPI Parallelization** and **BLAS/LAPACK Optimization** for high performance

### Upcoming Features

- **Spin-Orbit Coupling** by [**Relativistic** MRSF-TDDFT](https://doi.org/10.1021/acs.jctc.2c01036)
- **Ionization Potential/Electron Affinity** by [**EKT**-MRSF-TDDFT](https://doi.org/10.1021/acs.jpclett.1c02494)

### Quickstart

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

#### Compile

##### OpenMP Support

```bash
cd openqp
cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF
ninja -C build install
cd pyoqp
pip install .
```

##### OpenMP and MPI Support

```bash
cd openqp
cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF -DENABLE_MPI=ON
ninja -C build install
cd pyoqp
pip install .
```

##### OpenMP and MPI Support using make

```bash
cd openqp
cmake -B build -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF -DENABLE_MPI=ON
make -C build install
cd pyoqp
pip install .
```

- Use `-DUSE_LIBINT=ON` to replace the default ERI based on Rys Quadrature with `libint`.
- Use `-DLINALG_LIB_INT64=OFF` to ensure compatibility with third-party software like libdlfind compiled with 32-bit BLAS.

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

### Citing OpenQP
If you use OpenQP in your research, please cite the following papers:

- **Mironov V, Komarov K, Li J, Gerasimov I, Mazaheri M, Park W, Lashkaripour A, Oh M, Nakata H, Ishimura K, Huix-Rotllant M, Lee S, and Choi CH.** "OpenQP: A Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-source Ecosystem" [Submitted](https://doi.org/10.26434/chemrxiv-2024-k846p)
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

