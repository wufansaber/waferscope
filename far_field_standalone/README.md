# WaferScope

This directory contains the standalone launcher and PyInstaller configuration
for the particle-on-wafer far-field simulation UI.

Chinese user documentation: [USER_GUIDE_CN.md](USER_GUIDE_CN.md)

Independent 20x image-plane simulation: [IMAGING_MODEL.md](IMAGING_MODEL.md).
The default 20 mm objectives use 400 mm tube lenses. This mode preserves
SCATMECH Jones phase and computes three separate vector image fields.
Optional independent detectors integrate pixel areas and simulate electron
counts, shot/read noise and full-well saturation. See the image-model document
for detector assumptions and units. SCATMECH does not provide full near fields;
the verified native API limitation is recorded in [NEAR_FIELD_STATUS.md](NEAR_FIELD_STATUS.md).
The separate experimental FEM tab provides near/far/image calculations using
NGSolve. Windows releases now bundle its own solver executable alongside the GUI.

The simulator supports SCATMECH multilayer film stacks and optional
Rayleigh-Rice microroughness backgrounds in addition to the
Bobbert-Vlieger particle model.

## Run from source

Use the NIST2 environment from the repository root:

```powershell
D:\anaconda\envs\NIST2\python.exe far_field_standalone\main.py
```

## Build the distributable directory

Run `build.bat` after installing PyInstaller in both the NIST2 and isolated
NGSolve environments. The tested solver packager is PyInstaller 6.16.0.
`build_release.py --gui-python ... --fem-python ...` supports other interpreter paths.
The build uses two PyInstaller onedir runtimes and produces a dated ZIP plus SHA-256:

```text
dist\YYYY-MM-DD\WaferScope\
|-- WaferScope.exe
|-- _internal\
|-- fem\WaferScope_FEM.exe
|-- fem\_internal\
|-- s4\S4.exe
|-- s4\openblas.dll
|-- s4\source\
|-- examples\grating\
|-- docs\
|-- licenses\
`-- README_CN.md
```

Keep the executable and `_internal` directory together when distributing or
moving the application. The executable is not a one-file package and will not
run if it is copied away from its `_internal` dependencies. Keep `fem` alongside
the GUI as well, along with `s4` for periodic calculations. The frozen application writes results to
`%LOCALAPPDATA%\FarField_Simulation\results`, not the installation directory.
Release usage: [RELEASE_README_CN.md](RELEASE_README_CN.md).
Packaging checks: [PACKAGING_VALIDATION_CN.md](PACKAGING_VALIDATION_CN.md).

The S4 Grating (periodic) page includes parameter presets and complete saved
contact-sphere/grating demonstrations. Apply parameters prepares a new run;
Open example loads the original computed near/far/image/camera data. The
bundled FEM worker also hosts the isolated S4 process protocol, so no separate
Python environment is needed. These mixed-scene examples are exploratory,
not certified isolated-defect or finite-grating solutions.
