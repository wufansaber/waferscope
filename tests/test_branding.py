import os
from pathlib import Path


def test_main_window_displays_dfi_name():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from scatmech_ui.app import ScatmechStandaloneWindow
    app = QApplication.instance() or QApplication([])
    window = ScatmechStandaloneWindow()
    assert window.windowTitle() == 'DFI Simulation'
    window.close()
    app.processEvents()


def test_specs_and_current_manual_use_dfi_names():
    root = Path(__file__).resolve().parents[1]
    assert 'name="DFI_Simulation"' in (root/'far_field_standalone/FarField_Simulation.spec').read_text()
    assert "name='DFI_FEM'" in (root/'far_field_standalone/FarField_FEM.spec').read_text()
    manual = (root/'far_field_standalone/OPERATION_MANUAL_CN.md').read_text(encoding='utf-8')
    assert manual.startswith('# DFI Simulation')
    assert 'DFI_Simulation.exe' in manual and 'DFI_FEM.exe' in manual
