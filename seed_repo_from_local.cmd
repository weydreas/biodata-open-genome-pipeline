@echo off
setlocal EnableExtensions
chcp 65001 >nul 2>&1

set "ROOT=%~dp0"
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"
set "CACHE=%ROOT%\.pipeline_repo"
set "REPO=https://github.com/weydreas/biodata-open-genome-pipeline.git"

where git >nul 2>&1
if errorlevel 1 (
  echo ERROR: git was not found in PATH.
  pause
  exit /b 1
)

if not exist "%CACHE%\.git" (
  echo Cloning private pipeline repository...
  git clone "%REPO%" "%CACHE%"
  if errorlevel 1 (
    echo ERROR: clone failed. Complete Git Credential Manager authentication and rerun this file.
    pause
    exit /b 1
  )
) else (
  echo Updating existing local repository cache...
  git -C "%CACHE%" fetch origin main
  if errorlevel 1 goto :git_error
  git -C "%CACHE%" reset --hard origin/main
  if errorlevel 1 goto :git_error
)

for %%D in (src tests docs release) do (
  if not exist "%CACHE%\%%D" mkdir "%CACHE%\%%D"
)

robocopy "%ROOT%\src" "%CACHE%\src" /E /PURGE /XD __pycache__ /XF *.pyc >nul
if errorlevel 8 goto :copy_error
robocopy "%ROOT%\tests" "%CACHE%\tests" /E /PURGE /XD __pycache__ /XF *.pyc >nul
if errorlevel 8 goto :copy_error
robocopy "%ROOT%\docs" "%CACHE%\docs" /E /PURGE >nul
if errorlevel 8 goto :copy_error
robocopy "%ROOT%\release" "%CACHE%\release" /E /PURGE >nul
if errorlevel 8 goto :copy_error
copy /Y "%ROOT%\start.cmd" "%CACHE%\start.cmd" >nul
copy /Y "%ROOT%\bootstrap_git.cmd" "%CACHE%\bootstrap_git.cmd" >nul
copy /Y "%ROOT%\seed_repo_from_local.cmd" "%CACHE%\seed_repo_from_local.cmd" >nul
if exist "%ROOT%\.gitignore" copy /Y "%ROOT%\.gitignore" "%CACHE%\.gitignore" >nul
if exist "%ROOT%\README.md" copy /Y "%ROOT%\README.md" "%CACHE%\README.md" >nul

rem Never stage local patient/reference payloads even if somebody later copies them into the cache.
git -C "%CACHE%" add -- .gitignore README.md start.cmd bootstrap_git.cmd seed_repo_from_local.cmd release src tests docs
if errorlevel 1 goto :git_error

git -C "%CACHE%" diff --cached --quiet
if not errorlevel 1 (
  echo Repository already matches the local mutable source tree.
  goto :done
)

for /f "delims=" %%A in ('git -C "%CACHE%" config user.name 2^>nul') do set "GIT_NAME=%%A"
if not defined GIT_NAME git -C "%CACHE%" config user.name "Biodata Pipeline Bootstrap"
for /f "delims=" %%A in ('git -C "%CACHE%" config user.email 2^>nul') do set "GIT_EMAIL=%%A"
if not defined GIT_EMAIL git -C "%CACHE%" config user.email "biodata-pipeline@users.noreply.github.com"

git -C "%CACHE%" commit -m "Seed v9.0.0-hotfix5 mutable pipeline source"
if errorlevel 1 goto :git_error
git -C "%CACHE%" push origin main
if errorlevel 1 (
  echo.
  echo ERROR: push failed. If Git Credential Manager opens a browser, authenticate and rerun this file.
  pause
  exit /b 1
)

echo Full mutable source was pushed to GitHub main.

:done
echo Future start.cmd launches will fetch code updates automatically.
pause
exit /b 0

:copy_error
echo ERROR: robocopy failed.
pause
exit /b 1

:git_error
echo ERROR: git operation failed.
pause
exit /b 1
