[CmdletBinding()]
param(
    [string]$PackagePath,
    [string]$Msys2Root = $env:MSYS2_ROOT,
    [switch]$SkipPythonDeps,
    [switch]$NoDesktopLauncher
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-OpenQPMessage {
    param([string]$Message)
    Write-Host "[OpenQP] $Message"
}

function Find-Msys2Root {
    param([string]$PreferredRoot)

    $candidates = @()
    if ($PreferredRoot) {
        $candidates += $PreferredRoot
    }
    $candidates += "C:\msys64"
    $candidates += "C:\msys2"

    foreach ($candidate in ($candidates | Select-Object -Unique)) {
        if (-not $candidate) {
            continue
        }
        $root = [System.IO.Path]::GetFullPath($candidate)
        $bash = Join-Path $root "usr\bin\bash.exe"
        $python = Join-Path $root "ucrt64\bin\python.exe"
        if ((Test-Path -LiteralPath $bash) -and (Test-Path -LiteralPath $python)) {
            return $root
        }
    }

    throw "MSYS2 UCRT64 was not found. Install MSYS2 from https://www.msys2.org/ or rerun with -Msys2Root C:\msys64."
}

function Find-OpenQPPackage {
    param([string]$RequestedPackagePath)

    if ($RequestedPackagePath) {
        $resolved = Resolve-Path -LiteralPath $RequestedPackagePath
        return Get-Item -LiteralPath $resolved.Path
    }

    $searchRoots = @(
        $PSScriptRoot,
        (Join-Path $PSScriptRoot "dist\msys2-ucrt64\packages"),
        (Get-Location).Path,
        (Join-Path (Get-Location).Path "dist\msys2-ucrt64\packages")
    )

    $matches = foreach ($root in ($searchRoots | Select-Object -Unique)) {
        if (Test-Path -LiteralPath $root) {
            Get-ChildItem -LiteralPath $root -File -ErrorAction SilentlyContinue |
                Where-Object {
                    $_.Name -like "mingw-w64-ucrt-x86_64-openqp-*.pkg.tar.*" -or
                    $_.Name -like "openqp-msys2-ucrt64.pkg.tar.*"
                }
        }
    }

    $package = $matches | Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if (-not $package) {
        throw "Could not find the OpenQP .pkg.tar.* package next to this installer. Rerun with -PackagePath C:\path\to\package.pkg.tar.zst."
    }

    return $package
}

function Invoke-Ucrt64 {
    param(
        [string]$BashPath,
        [string]$Command,
        [string]$Description
    )

    Write-OpenQPMessage $Description
    $env:MSYSTEM = "UCRT64"
    $env:CHERE_INVOKING = "1"
    & $BashPath -lc $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Description failed."
    }
}

function New-DesktopLauncher {
    param(
        [string]$BashPath
    )

    $desktop = [Environment]::GetFolderPath("DesktopDirectory")
    if (-not $desktop) {
        return
    }

    $launcher = Join-Path $desktop "OpenQP UCRT64.cmd"
    $escapedBash = $BashPath.Replace('"', '""')
    $content = @"
@echo off
set "MSYSTEM=UCRT64"
set "CHERE_INVOKING=1"
"$escapedBash" -lc "cd ~; echo OpenQP is ready.; echo Try: openqp --help; exec bash"
"@
    Set-Content -LiteralPath $launcher -Value $content -Encoding ASCII
    Write-OpenQPMessage "Created launcher: $launcher"
}

Write-OpenQPMessage "Starting MSYS2 UCRT64 installer."
$root = Find-Msys2Root -PreferredRoot $Msys2Root
$bash = Join-Path $root "usr\bin\bash.exe"
$package = Find-OpenQPPackage -RequestedPackagePath $PackagePath

Write-OpenQPMessage "Using MSYS2 root: $root"
Write-OpenQPMessage "Using package: $($package.FullName)"

$env:OPENQP_PKG_WIN = $package.FullName

Invoke-Ucrt64 -BashPath $bash -Description "Installing pacman support packages" -Command @'
pacman --needed --noconfirm -S mingw-w64-ucrt-x86_64-python-pip
'@

Invoke-Ucrt64 -BashPath $bash -Description "Installing OpenQP package" -Command @'
pkg=$(cygpath -u "$OPENQP_PKG_WIN")
pacman -U --noconfirm "$pkg"
'@

if (-not $SkipPythonDeps) {
    Invoke-Ucrt64 -BashPath $bash -Description "Installing temporary Python runtime dependencies" -Command @'
python -m pip install --user --upgrade basis_set_exchange geometric
'@
}

Invoke-Ucrt64 -BashPath $bash -Description "Checking OpenQP import and launcher" -Command @'
python -c "import oqp; print(oqp.oqp_root)"
command -v openqp
'@

if (-not $NoDesktopLauncher) {
    New-DesktopLauncher -BashPath $bash
}

Write-OpenQPMessage "Done. Open a UCRT64 shell and run: openqp --help"
