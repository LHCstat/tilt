@echo off
setlocal
cd /d "%~dp0"
set "PYTHONUTF8=1"
set "OMP_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
set "TILT_BOOTSTRAP=py -3.11"
%TILT_BOOTSTRAP% -c "import sys,struct; assert sys.version_info[:2]==(3,11) and struct.calcsize('P')==8" >nul 2>&1
if errorlevel 1 set "TILT_BOOTSTRAP=python"
%TILT_BOOTSTRAP% -c "import sys,struct; assert sys.version_info[:2]==(3,11) and struct.calcsize('P')==8" >nul 2>&1
if errorlevel 1 goto missing_python
if not exist ".venv\Scripts\python.exe" %TILT_BOOTSTRAP% -m venv .venv
if errorlevel 1 goto failed
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto failed
".venv\Scripts\python.exe" -m unittest discover -s . -p "test_*.py" -v
if errorlevel 1 goto failed
echo Installation and tests completed. Run demo_windows.bat or run_windows.bat.
if /I not "%~1"=="--no-pause" pause
exit /b 0
:missing_python
echo Install Python 3.11 64-bit from python.org, including the Python launcher.
goto failed
:failed
echo Setup failed. Read the error above and README.md.
if /I not "%~1"=="--no-pause" pause
exit /b 1
