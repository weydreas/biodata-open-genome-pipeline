@echo off
setlocal
chcp 65001 >nul 2>&1
set "ROOT=%~dp0"
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"
set "CACHE=%ROOT%\.pipeline_repo"

where git >nul 2>&1
if errorlevel 1 (
  echo ERROR: git was not found in PATH.
  pause
  exit /b 1
)

where py >nul 2>&1
if not errorlevel 1 (
  set "PY=py"
) else (
  where python >nul 2>&1
  if errorlevel 1 (
    echo ERROR: Python 3 was not found.
    pause
    exit /b 1
  )
  set "PY=python"
)

if not exist "%CACHE%\.git" (
  echo Initializing Biodata genome pipeline code cache...
  git clone --depth 1 --branch main https://github.com/weydreas/biodata-open-genome-pipeline.git "%CACHE%"
  if errorlevel 1 (
    echo ERROR: GitHub clone failed. For the private repository, complete Git Credential Manager authentication and retry.
    pause
    exit /b 1
  )
)

powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%CACHE%\src\bootstrap_update.ps1" -Root "%ROOT%" -PythonCommand "%PY%"
if errorlevel 1 (
  echo ERROR: bootstrap update failed.
  pause
  exit /b 1
)

echo.
echo Bootstrap complete. Future launches through start.cmd will check GitHub and reference bundles automatically.
pause
endlocal
