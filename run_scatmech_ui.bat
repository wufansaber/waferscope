@echo off
set "CONDA_BAT=D:\anaconda\condabin\conda.bat"
set "NIST2_PYTHON=D:\anaconda\envs\NIST2\python.exe"

if exist "%CONDA_BAT%" (
    call "%CONDA_BAT%" run -n NIST2 python "%~dp0main.py"
) else if exist "%NIST2_PYTHON%" (
    "%NIST2_PYTHON%" "%~dp0main.py"
) else (
    echo NIST2 Python environment was not found.
    echo Expected: %NIST2_PYTHON%
    exit /b 1
)
