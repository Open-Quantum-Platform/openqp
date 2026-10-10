$ErrorActionPreference = 'Stop'
if ($env:STUDIO_BUILD_ENGINE -ne '1') { throw 'Engine build must be explicitly selected' }
# One release job owns this checkout and cache; cap both native and Rust builds.
(Get-Process -Id $PID).ProcessorAffinity = 15
$env:CMAKE_BUILD_PARALLEL_LEVEL = '4'
$env:CARGO_BUILD_JOBS = '4'
$env:OMP_NUM_THREADS = '4'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$oneapi = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$vcvars = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
$envLines = & cmd.exe /d /s /c "call `"$vcvars`" >nul && call `"$oneapi`" --force >nul && set"
if ($LASTEXITCODE -ne 0) { throw 'Compiler environment failed' }
foreach ($line in $envLines) {
  if ($line -match '^([^=]+)=(.*)$') {
    [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
  }
}
$python = 'C:\Program Files\Python311\python.exe'
$config = Get-Content release-candidate.json -Raw | ConvertFrom-Json
& $python tools/checkout_engine.py --ref $config.engine.commit
if ($LASTEXITCODE -ne 0) { throw 'Pinned gateway checkout failed' }
& $python -m venv .cache/engine-venv
$venvPython = "$PWD\.cache\engine-venv\Scripts\python.exe"
Copy-Item $venvPython "$PWD\.cache\engine-venv\Scripts\python3.exe" -Force
$env:PATH = "$PWD\.cache\engine-venv\Scripts;$env:PATH"
$env:CC = 'icx'; $env:CXX = 'icx'; $env:FC = 'ifx'
$env:CMAKE_GENERATOR = 'Ninja'
# Separate from the gateway runner's shared native cache.
$cache = 'C:\oqp-studio-024\externals'
$build = "C:\oqp-studio-024\build-$env:CI_JOB_ID"
$temp = "C:\oqp-studio-024\temp-$env:CI_JOB_ID"
New-Item -ItemType Directory -Force -Path $cache,$build,$temp | Out-Null
$env:TEMP = $temp; $env:TMP = $temp
& $venvPython -m pip install --upgrade pip cmake ninja 'numpy<2.2' cffi scikit-build-core pyinstaller==6.19.0
if ($LASTEXITCODE -ne 0) { throw 'Build requirements failed' }
& $venvPython -m pip install -v ./openqp --no-build-isolation --no-deps `
  --config-settings=build-dir=$build `
  --config-settings=cmake.define.LINALG_LIB=Intel10_64ilp `
  --config-settings=cmake.define.ENABLE_MPI=OFF `
  --config-settings=cmake.define.ENABLE_OPENMP=ON `
  --config-settings=cmake.define.USE_LIBINT=OFF `
  --config-settings=cmake.define.ENABLE_OPENTRAH=OFF `
  '--config-settings=cmake.define.OQP_DFTD4_PATCH_EXECUTABLE=C:\Program Files\Git\usr\bin\patch.exe' `
  --config-settings=cmake.define.OQP_REUSE_EXTERNALS=ON `
  --config-settings=cmake.define.OQP_EXTERNALS_ROOT=$cache
if ($LASTEXITCODE -ne 0) { throw 'Native MKL ILP64 engine build failed' }
& $venvPython -m pip install -r openqp/pyoqp/requirements.txt 'numpy<2.2' basis_set_exchange
if ($LASTEXITCODE -ne 0) { throw 'Runtime requirements failed' }
& $venvPython tools/release/freeze_engine.py --source openqp --output .cache/engine-frozen/openqp
if ($LASTEXITCODE -ne 0) { throw 'Engine freezing failed' }
& $venvPython tools/release/qualify_engine.py .cache/engine-frozen/openqp .cache/engine-qualification
if ($LASTEXITCODE -ne 0) { throw 'NMR/ACID qualification failed' }
& $python -m venv .cache/package-venv
& .cache/package-venv/Scripts/python.exe tools/release/package.py --engine .cache/engine-frozen/openqp
if ($LASTEXITCODE -ne 0) { throw 'Installer packaging failed' }
