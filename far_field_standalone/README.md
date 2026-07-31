# Far-field Standalone UI

This directory contains the standalone launcher and PyInstaller configuration
for the particle-on-wafer far-field simulation UI.

Chinese user documentation: [USER_GUIDE_CN.md](USER_GUIDE_CN.md)

The simulator supports SCATMECH multilayer film stacks and optional
Rayleigh-Rice microroughness backgrounds in addition to the
Bobbert-Vlieger particle model.

## Run from source

Use the NIST2 environment from the repository root:

```powershell
D:\anaconda\envs\NIST2\python.exe far_field_standalone\main.py
```

## Build the distributable directory

Run `build.bat`. The build uses PyInstaller onedir mode and produces:

```text
dist\FarField_Simulation\
|-- FarField_Simulation.exe
`-- _internal\
```

Keep the executable and `_internal` directory together when distributing or
moving the application. The executable is not a one-file package and will not
run if it is copied away from its `_internal` dependencies.
