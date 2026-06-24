## Open Quantum Platform: OpenQP

Open Quantum Platform ([OpenQP](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)) is a quantum chemical platform featuring cutting-edge capabilities like [Mixed-Reference Spin-Flip (MRSF)-TDDFT](https://doi.org/10.1021/acs.jpclett.3c02296) with an emphasis on open-source ecosystem.

### Web tools

- [Try OpenQP Web](https://app.openqp.org/) - browser app for choosing example workflows, previewing structures, and downloading OpenQP input files. Calculations are prepared locally in the browser and are not executed online yet.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/) - existing browser input generator with a 3D molecular editor.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) - browser-based viewer for OpenQP logs, JSON, Molden, cube, and XYZ data.

### Capabilities

**Electronic structure**
- HF and DFT (LibXC functionals) with RHF, ROHF, and UHF references
- TDHF/TDDFT, SF-TDDFT, and **MRSF-TDDFT** including the DTCAM-series functionals
- **UMRSF-TDDFT** excitation energies from a UHF reference (energy-only)
- **MRSF-EKT** ionization potentials and electron affinities with Dyson orbitals and pole strengths

**Properties & spectroscopy**
- Analytic gradients; native analytic HF/DFT Hessians covering open-shell (UHF/ROHF) references, ECPs, and range-separated (CAM/LRC) functionals
- Vibrational frequencies, thermochemistry, and native **IR and Raman intensities**
- **NMR chemical shieldings** with CGO and GIAO (London-orbital) gauge treatments
- MRSF-TDDFT **nonadiabatic couplings** (NACME) for dynamics workflows
- **Spin-orbit couplings (SOC)** between MRSF-TDDFT states from one- and two-electron contributions ([Relativistic MRSF-TDDFT](https://doi.org/10.1021/acs.jctc.2c01036))
- **Implicit solvation**: energy-only ddPCM continuum solvent via the ddX library

**Geometry & reaction paths**
- Minima, transition states, MECI/MECP, constrained optimization, IRC, and NEB
- Built-in native optimizer (`lib=oqp`, TRIC/DLC + restricted-step RFO), with geomeTRIC and SciPy available as optional backends. DL-FIND is no longer used in the Python startup/dispatch path.

**Dynamics & QM/MM**
- Native fewest-switches surface hopping (`runtype=namd`) for gas-phase MRSF-TDDFT internal conversion.
- SOC-NAMD for intersystem crossing, including SHARC-like spin-adiabatic propagation and an MCH-basis SOC mode with exact active-root MCH gradients (`[md] soc_basis=mch`).
- ESPF electrostatic QM/MM with OpenMM, including PME periodic electrostatics, smooth ESPF grid forces, and QM/MM NAMD/SOC-NAMD dispatch.
- Overlap-based MRSF state tracking, finite-time NAC/TDC propagation, and SOC-QM/MM regression guards for hop bookkeeping and restart/reproducibility-sensitive state.

**Reliability & performance**
- Point-group **symmetry**: detection, MO/state/mode labels, and petite-list reductions accelerating integrals, XC, gradients, and response
- DFT grids: Lebedev plus SG-0/SG-1/SG-2/SG-3 pruned grids with per-element DE2 radial quadrature; OpenMP-parallel XC kernels
- Native initial guesses (`hcore`, `huckel`, `modhuckel`, `minao`, `sap`), JSON restart, OpenTrustRegion (TRAH) SCF stabilization, Davidson auto-restart, and MINRES/AUTO Z-vector fallbacks
- OpenMP/MPI parallelism, ILP64 BLAS/LAPACK, pip installs, and Docker images

**Ecosystem**
- LibXC, basis_set_exchange, libecpint, DFT-D4, PyRAI2MD, [OpenQP Web](https://app.openqp.org/), the [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/), and the browser-based [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) viewer

### Upcoming Features
- Full analytic spin-adiabatic SOC gradients, requiring MCH derivative-coupling vectors and SOC-gradient matrix elements.
- **Scalar-relativistic (X2C) framework** extending the relativistic MRSF-TDDFT treatment

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

Two representative configurations are shown below; both build the native library and install the Python package. See **Common options** for the remaining build switches.

##### Linux (GCC + OpenMP)

```bash
cd openqp
cmake -B build -G Ninja \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=. \
  -DENABLE_OPENMP=ON
ninja -C build install
cd pyoqp
pip install .
```

##### macOS (Homebrew GCC + Accelerate BLAS)

Use the Homebrew GCC toolchain and the native Accelerate BLAS/LAPACK provider (32-bit integer interface). Adjust the `-15` suffix to your installed Homebrew GCC.

```bash
cd openqp
cmake -B build -G Ninja \
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

##### Common options

| Option | Default | Description |
|--------|---------|-------------|
| `-DENABLE_MPI=ON` | OFF | Build with MPI; pair with an MPI Fortran wrapper, e.g. `-DCMAKE_Fortran_COMPILER=mpif90`. |
| `-DENABLE_OPENMP=ON` | ON | OpenMP threading. |
| `-DUSE_LIBINT=ON` | OFF | Use `libint` for ERIs instead of the built-in Rys-quadrature engine. |
| `-DLINALG_LIB=<vendor>` | auto | BLAS/LAPACK provider (`auto`, `MKL`, `OpenBLAS`, `netlib`, …). |
| `-DLINALG_LIB_INT64=ON` | ON | Use the default 64-bit-integer BLAS/LAPACK interface. |
| `-DENABLE_OPENTRAH=OFF` | ON | Skip the external OpenTrustRegion (OpenTRAH) solver; the TRAH SCF converger then uses the built-in native implementation only. |
| `-DOQP_REUSE_EXTERNALS=OFF` | ON | Disable the automatic reusable cache for bundled external dependency builds. |

Bundled external dependencies are cached automatically across fresh OpenQP build directories. The cache uses the platform cache directory (`~/Library/Caches/openqp/externals` on macOS, `$XDG_CACHE_HOME/openqp/externals` or `~/.cache/openqp/externals` elsewhere) and is keyed by operating system, architecture, CMake generator, compiler IDs/versions, compiler flags, build type, linkage mode, BLAS/LAPACK choice, BLAS/LAPACK integer ABI, and bundled external versions, so Fortran modules and libraries are reused only by compatible builds. Advanced users can relocate the cache with `-DOQP_EXTERNALS_ROOT=<path>`. The cache is not safe for concurrent use: do not run two fresh OpenQP builds with the same configuration at the same time, or configure one of them with `-DOQP_REUSE_EXTERNALS=OFF`.

LP64 BLAS/LAPACK (`-DLINALG_LIB_INT64=OFF`) is a macOS-only exception intended for a consistent native Accelerate stack. OpenBLAS and MKL builds should use ILP64 for normal use.

To build without Ninja, drop `-G Ninja` and replace `ninja` with `make`.

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

##### Setting the OpenMP thread count

The number of OpenMP threads (per process / MPI rank) can be set directly from
the input file or the command line, instead of relying on `OMP_NUM_THREADS`:

```ini
[input]
omp_threads=16        # OpenMP threads per process / MPI rank
```

```bash
openqp any_example_file.inp --omp 16
```

Precedence (highest first): `--omp` → input `omp_threads` → `OMP_NUM_THREADS` →
built-in default. The value is applied before the OpenMP runtime initializes, and
is honored on the programmatic `Runner(input_dict=...)` path as well. If OpenQP
was built without OpenMP (`-DENABLE_OPENMP=OFF`), the request is ignored with a
warning and the run is serial. MPI rank count remains launcher-controlled
(`mpirun -np ...`); `--nompi` disables MPI.

##### NAMD, SOC-NAMD, and QM/MM

Surface-hopping dynamics use `runtype=namd` with an `[md]` section. For
gas-phase internal conversion, set `soc=false`; for SOC dynamics, set
`soc=true`. The recommended SOC-QM/MM production mode is the MCH-basis SOC path:

```ini
[input]
runtype=namd
qmmm_flag=true

[md]
soc=true
soc_basis=mch       # recommended SOC-QM/MM production mode
dt=0.5
nstep=200
```

The SHARC-like spin-adiabatic path remains available with
`soc_basis=adiabatic`. Optional diagnostic corrections are controlled with
`soc_du_dt_corr` and `soc_tdc_grad_corr`. QM/MM single-point energy, gradient,
and optimization jobs stay on the normal `Runner` path; only legacy
ground-state OpenMM molecular dynamics uses `runtype=md` with `qmmm_flag=true`.

#### Performance and threading

OpenQP parallelizes the two-electron integral / Fock build with OpenMP. A few
defaults and knobs keep it efficient:

- **Sequential BLAS inside the integral build (automatic).** A threaded BLAS
  keeps its own worker pool that would oversubscribe the cores against the
  OpenMP integral threads. OpenQP forces BLAS to a single thread for the
  duration of the integral build (via the BLAS library's own runtime setter, so
  it cannot be re-introduced by a stray `OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`)
  and restores full BLAS threading afterward for diagonalization and other
  BLAS-heavy phases. The frontend also sets conservative `*_NUM_THREADS=1` and
  `OMP_STACKSIZE` defaults (via `setdefault`, so an explicit environment wins).
  Measured ~1.5x faster at 24-28 threads on a dedicated node; near-ideal scaling
  to the performance-core count.

- **Low-memory Fock accumulator (automatic for large, sparse systems).** The
  Fock build normally keeps one Fock copy per thread
  (`fockdim x nfocks x nthreads`); for large bases on many cores this can reach
  several GB. OpenQP can instead accumulate into a single shared Fock with atomic
  updates (~`Nthreads`x less memory). It auto-engages only when the replicated
  buffers would exceed a cap **and** the density is sparse (so atomic contention
  is negligible). Controls:
  - `OQP_FOCK_MEM_MB` — replicated-memory cap in MB before switching (default `4096`).
  - `OQP_FOCK_SPARSITY` — significant shell-pair fraction below which the density
    is considered sparse (default `0.5`).
  - `OQP_FOCK_ATOMIC=1`/`0` — force the shared-atomic / replicated mode regardless.

### Detailed Documentation

For more in-depth information, visit:
- [OpenQP Documentation](https://github.com/Open-Quantum-Platform/openqp/wiki)

### OpenQP Web
Prepare OpenQP inputs in the browser at [app.openqp.org](https://app.openqp.org/). The web app provides example workflows, local structure preview, downloadable `.inp` and `.xyz` files, tutorials, and links to the existing OpenQP browser tools. It does not run OpenQP jobs online, require login, or expose API keys in the browser.

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
- **Vladimir Makhnev**, [VladimirMakhnev](https://github.com/VladimirMakhnev)
- **Alireza Lashkaripour**, [Alireza Lashkaripour](https://github.com/Alireza-Lashkaripour)
### Legal Notice

See the separate LICENSE file.
