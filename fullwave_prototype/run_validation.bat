@echo off
setlocal
cd /d "%~dp0.."
if not exist "build\ngsolve_env\Scripts\python.exe" (
  echo Create the isolated Python environment described in fullwave_prototype\README.md first.
  exit /b 1
)
"build\ngsolve_env\Scripts\python.exe" -m fullwave_prototype.run %*
exit /b %errorlevel%
