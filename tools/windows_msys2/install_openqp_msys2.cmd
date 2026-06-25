@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "INSTALLER=%SCRIPT_DIR%install_openqp_msys2.ps1"

if not exist "%INSTALLER%" (
  echo OpenQP installer script was not found:
  echo   %INSTALLER%
  echo.
  pause
  exit /b 1
)

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%INSTALLER%" %*
set "STATUS=%ERRORLEVEL%"

echo.
if not "%STATUS%"=="0" (
  echo OpenQP installation failed with exit code %STATUS%.
  pause
  exit /b %STATUS%
)

echo OpenQP installation finished.
pause
exit /b 0
