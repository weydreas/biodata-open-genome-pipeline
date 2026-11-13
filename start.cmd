@echo off
setlocal
chcp 65001 >nul 2>&1

set "DANTE_ROOT=%~dp0"
rem %~dp0 always ends with a backslash; remove it before passing as a quoted argument.
if "%DANTE_ROOT:~-1%"=="\" set "DANTE_ROOT=%DANTE_ROOT:~0,-1%"

set "DANTE_SRC=%DANTE_ROOT%\src"
set "DANTE_LAYOUT_MANAGER=%DANTE_SRC%\layout_manager.py"
set "DANTE_PATIENT_MANAGER=%DANTE_SRC%\patient_manager.py"
set "DANTE_PATIENT_DISCOVERY=%DANTE_SRC%\patient_folder_discovery.py"
set "DANTE_PS1=%DANTE_SRC%\dante_wgs_pipeline.ps1"

if not exist "%DANTE_LAYOUT_MANAGER%" (
  echo ERROR: src\layout_manager.py was not found.
  echo Extract the full v9.0.0 archive over the pipeline root.
  pause
  exit /b 1
)
if not exist "%DANTE_PATIENT_MANAGER%" (
  echo ERROR: src\patient_manager.py was not found.
  pause
  exit /b 1
)
if not exist "%DANTE_PS1%" (
  echo ERROR: src\dante_wgs_pipeline.ps1 was not found.
  pause
  exit /b 1
)

where py >nul 2>&1
if not errorlevel 1 (
  set "DANTE_PYTHON=py"
) else (
  where python >nul 2>&1
  if errorlevel 1 (
    echo ERROR: Python 3 was not found. Install Python or make py/python available in PATH.
    pause
    exit /b 1
  )
  set "DANTE_PYTHON=python"
)

set "DANTE_BOOTSTRAP=%DANTE_SRC%\bootstrap_update.ps1"
if exist "%DANTE_BOOTSTRAP%" (
  powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%DANTE_BOOTSTRAP%" -Root "%DANTE_ROOT%" -PythonCommand "%DANTE_PYTHON%"
  if errorlevel 1 (
    echo WARNING: pipeline/reference auto-update did not complete. Continuing with local files.
  )
) else (
  echo WARNING: src\bootstrap_update.ps1 not found; automatic GitHub/reference update is disabled.
)

rem Final report items 25/26 are deterministic Stage-24 PDF renderers and do not
rem require the OpenAI Python SDK or any API preflight.

rem v9.0.0: retain the one-time v8.6.0 flat-layout migration before patient selection.
"%DANTE_PYTHON%" "%DANTE_LAYOUT_MANAGER%" migrate --root "%DANTE_ROOT%"
if errorlevel 1 (
  echo Code-layout migration failed. Patient data were not processed.
  pause
  endlocal
  exit /b 1
)

rem Adopt manually created immediate folders under patients\ before showing the menu.
rem This only creates missing patient metadata/standard subfolders; existing files are
rem never moved or overwritten. Source genome files still belong in patients\NAME\input\.
if exist "%DANTE_PATIENT_DISCOVERY%" (
  "%DANTE_PYTHON%" "%DANTE_PATIENT_DISCOVERY%" --root "%DANTE_ROOT%"
  if errorlevel 1 (
    echo Patient-folder discovery failed. No pipeline action was started.
    pause
    endlocal
    exit /b 1
  )
) else (
  echo WARNING: src\patient_folder_discovery.py not found; manually created patient folders without patient.json will not be auto-adopted.
)

rem Select patient BEFORE any READY/INCOMPLETE status is calculated.
"%DANTE_PYTHON%" "%DANTE_PATIENT_MANAGER%" select --root "%DANTE_ROOT%"
set "PATIENT_EXIT=%ERRORLEVEL%"
if "%PATIENT_EXIT%"=="10" (
  endlocal
  exit /b 0
)
if not "%PATIENT_EXIT%"=="0" (
  echo Patient selection/migration failed. No pipeline action was started.
  pause
  endlocal
  exit /b 1
)

powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -Command "$tokens = $null; $parseErrors = $null; [void][System.Management.Automation.Language.Parser]::ParseFile($env:DANTE_PS1, [ref]$tokens, [ref]$parseErrors); if ($parseErrors.Count -gt 0) { $parseErrors | ForEach-Object { Write-Host ('PowerShell parser error: ' + $_.Message) -ForegroundColor Red }; exit 2 }"
if errorlevel 1 (
  echo Pipeline syntax preflight failed. No pipeline action was started.
  pause
  endlocal
  exit /b 1
)

rem Interactive contract: each menu selection runs in a fresh PowerShell process.
rem Successful actions immediately reopen the menu with freshly recalculated
rem READY/STALE states. Menu item 0 returns sentinel exit code 10.
:PIPELINE_MENU
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%DANTE_PS1%" -PerformanceProfile max -NoPause
set "PIPELINE_EXIT=%ERRORLEVEL%"
if "%PIPELINE_EXIT%"=="10" (
  endlocal
  exit /b 0
)
if not "%PIPELINE_EXIT%"=="0" (
  echo Pipeline action failed with exit code %PIPELINE_EXIT%.
  pause
  endlocal
  exit /b %PIPELINE_EXIT%
)
goto PIPELINE_MENU
