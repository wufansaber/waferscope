@echo off
setlocal

set "PYTHON_EXE=D:\anaconda\envs\NIST2\python.exe"

if not exist "%PYTHON_EXE%" (
    echo NIST2 Python was not found at: %PYTHON_EXE%
    exit /b 1
)

cd /d "%~dp0"
"%PYTHON_EXE%" build_release.py
if errorlevel 1 exit /b %errorlevel%

echo.
echo Build completed:
echo %~dp0dist\YYYY-MM-DD\WaferScope\WaferScope.exe
endlocal
