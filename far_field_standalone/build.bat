@echo off
setlocal

set "PYTHON_EXE=D:\anaconda\envs\NIST2\python.exe"

if not exist "%PYTHON_EXE%" (
    echo NIST2 Python was not found at: %PYTHON_EXE%
    exit /b 1
)

cd /d "%~dp0"
"%PYTHON_EXE%" -m PyInstaller FarField_Simulation.spec --clean --noconfirm
if errorlevel 1 exit /b %errorlevel%
copy /Y "USER_GUIDE_CN.md" "dist\FarField_Simulation\USER_GUIDE_CN.md" >nul
if errorlevel 1 exit /b %errorlevel%

echo.
echo Build completed:
echo %~dp0dist\FarField_Simulation\FarField_Simulation.exe
endlocal
