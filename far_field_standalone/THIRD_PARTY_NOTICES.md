# Third-party components

This distribution bundles third-party shared libraries and Python modules,
plus an explicitly modified S4 native solver. Python package licenses and version inventories are copied from the installed
distributions into `licenses/gui` and `licenses/fem`. License terms remain those
of the original projects; this notice does not replace them.

Main upstream projects and source locations:

- Python: https://www.python.org/downloads/source/
- PyQt5: https://www.riverbankcomputing.com/software/pyqt/download5
- Qt 5: https://download.qt.io/archive/qt/
- pySCATMECH: https://github.com/usnistgov/pySCATMECH
- SCATMECH: https://github.com/usnistgov/SCATMECH
- NGSolve 6.2.2607: https://github.com/NGSolve/ngsolve
- Netgen 6.2.2607: https://github.com/NGSolve/netgen
- Open CASCADE: https://github.com/Open-Cascade-SAS/OCCT
- Netgen OCCT wheel 7.8.1: https://pypi.org/project/netgen-occt/7.8.1/
- OpenBLAS: https://github.com/OpenMathLib/OpenBLAS
- NumPy: https://github.com/numpy/numpy
- SciPy: https://github.com/scipy/scipy
- Matplotlib: https://github.com/matplotlib/matplotlib
- psutil: https://github.com/giampaolo/psutil
- PyInstaller: https://github.com/pyinstaller/pyinstaller
- S4 (GPL v2-or-later), pinned source plus documented project adapter:
  https://github.com/victorliu/S4
- Lua 5.2.4: https://www.lua.org/license.html

The modified S4 runtime is a separate command-line process in `s4/`.
Its license texts, modified corresponding source, pristine S4/Lua archives,
patches and build instructions are included under `s4/licenses` and `s4/source`.
It does not become an unmodified upstream build by being packaged here.

NGSolve/Netgen/OCCT libraries remain separate DLLs, not statically linked into
application code. The directories must remain intact for normal operation.
Applicable LGPL rights, including replacement/debugging of modified libraries,
are not restricted by this application. Source application and build recipes
are distributed in the same GitHub repository as this release.
