@echo off
setlocal
set "PYTHONUTF8=1"
set "OMP_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
if not exist "%~dp0.venv\Scripts\python.exe" goto missing_env
"%~dp0.venv\Scripts\python.exe" -X utf8 "%~dp0run_local.py" %*
set "TILT_EXIT=%ERRORLEVEL%"
pause
exit /b %TILT_EXIT%
:missing_env
echo Run setup_windows.bat first.
pause
exit /b 1
