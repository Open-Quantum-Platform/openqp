[CmdletBinding()]
param(
    [string]$PackagePath,
    [string]$PayloadDir = "dist\msys2-ucrt64\msi-payload",
    [string]$OutputDir = "dist\msys2-ucrt64\msi",
    [string]$ProductVersion
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$payload = Join-Path $repoRoot $PayloadDir
$output = Join-Path $repoRoot $OutputDir
$wxs = Join-Path $PSScriptRoot "msi\Product.wxs"

if (-not $ProductVersion) {
    $pyproject = Get-Content -LiteralPath (Join-Path $repoRoot "pyproject.toml") -Raw
    $match = [regex]::Match($pyproject, '(?m)^version\s*=\s*"([^"]+)"')
    if (-not $match.Success) {
        throw "Could not read project version from pyproject.toml."
    }
    $ProductVersion = $match.Groups[1].Value
}

if (-not $PackagePath) {
    $package = Get-ChildItem -LiteralPath (Join-Path $repoRoot "dist\msys2-ucrt64\packages") -File -Filter "mingw-w64-ucrt-x86_64-openqp-*.pkg.tar.*" |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
} else {
    $package = Get-Item -LiteralPath $PackagePath
}

if (-not $package) {
    throw "OpenQP MSYS2 package was not found. Build the pacman package first or pass -PackagePath."
}

New-Item -ItemType Directory -Force -Path $payload, $output | Out-Null
Copy-Item -LiteralPath (Join-Path $PSScriptRoot "install_openqp_msys2.cmd") -Destination (Join-Path $payload "install_openqp_msys2.cmd") -Force
Copy-Item -LiteralPath (Join-Path $PSScriptRoot "install_openqp_msys2.ps1") -Destination (Join-Path $payload "install_openqp_msys2.ps1") -Force
Copy-Item -LiteralPath $package.FullName -Destination (Join-Path $payload "openqp-msys2-ucrt64.pkg.tar.zst") -Force

@"
OpenQP MSYS2 Installer

This MSI installs the OpenQP MSYS2 installer bundle. It does not install MSYS2
itself. After the MSI finishes, run "Install OpenQP for MSYS2" from the Start
Menu or Desktop.

Prerequisite:
- MSYS2 UCRT64 installed at C:\msys64, or run the installer script manually with
  -Msys2Root.
"@ | Set-Content -LiteralPath (Join-Path $payload "README.txt") -Encoding ASCII

$wix = Get-Command wix -ErrorAction SilentlyContinue
if (-not $wix) {
    throw "WiX Toolset was not found. Install with: dotnet tool install --global wix"
}

& $wix.Source build $wxs `
    -arch x64 `
    -d "PayloadDir=$payload" `
    -d "ProductVersion=$ProductVersion" `
    -o (Join-Path $output "OpenQP-MSYS2-UCRT64-$ProductVersion.msi")

if ($LASTEXITCODE -ne 0) {
    throw "WiX MSI build failed."
}

Write-Host "Created MSI:" (Join-Path $output "OpenQP-MSYS2-UCRT64-$ProductVersion.msi")
