# Third-party components

This distribution bundles unmodified third-party shared libraries and Python
modules. Their licenses and version inventories are copied from the installed
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

NGSolve/Netgen/OCCT libraries remain separate DLLs, not statically linked into
application code. The directories must remain intact for normal operation.
Applicable LGPL rights, including replacement/debugging of modified libraries,
are not restricted by this application. Source application and build recipes
are distributed in the same GitHub repository as this release.
