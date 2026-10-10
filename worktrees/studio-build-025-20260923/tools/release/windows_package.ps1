$ErrorActionPreference = 'Stop'
if ($env:STUDIO_RELEASE_CANDIDATE -ne '1') { throw 'Candidate packaging must be explicitly selected' }
(Get-Process -Id $PID).ProcessorAffinity = 15
$env:CMAKE_BUILD_PARALLEL_LEVEL = '4'
$env:CARGO_BUILD_JOBS = '4'
$env:OMP_NUM_THREADS = '4'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$vcvars = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
$envLines = & cmd.exe /d /s /c "call `"$vcvars`" >nul && set"
if ($LASTEXITCODE -ne 0) { throw 'Compiler environment failed' }
foreach ($line in $envLines) {
  if ($line -match '^([^=]+)=(.*)$') {
    [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
  }
}
& 'C:\Program Files\Python311\python.exe' tools/release/ci_package.py
if ($LASTEXITCODE -ne 0) { throw 'Installer packaging failed' }
