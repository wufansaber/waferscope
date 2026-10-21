import os
from pathlib import Path


def test_main_window_displays_waferscope_name():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from scatmech_ui.app import ScatmechStandaloneWindow
    app = QApplication.instance() or QApplication([])
    window = ScatmechStandaloneWindow()
    assert window.windowTitle() == 'WaferScope'
    window.close()
    app.processEvents()


def test_specs_and_current_manual_use_waferscope_names():
    root = Path(__file__).resolve().parents[1]
    assert 'name="WaferScope"' in (root/'far_field_standalone/WaferScope.spec').read_text()
    assert "name='WaferScope_FEM'" in (root/'far_field_standalone/WaferScope_FEM.spec').read_text()
    manual = (root/'far_field_standalone/OPERATION_MANUAL_CN.md').read_text(encoding='utf-8')
    assert manual.startswith('# WaferScope')
    assert 'WaferScope.exe' in manual and 'WaferScope_FEM.exe' in manual


def test_platform_entry_keeps_legacy_class_import():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from main_window import APTKSimulationUI, WaferScopeMainWindow
    assert APTKSimulationUI is WaferScopeMainWindow
    app = QApplication.instance() or QApplication([])
    window = WaferScopeMainWindow()
    assert window.windowTitle() == 'WaferScope'
    window.close()
    app.processEvents()


def test_migration_links_use_current_repository_and_archive_name():
    root = Path(__file__).resolve().parents[1]
    readme = (root/'README.md').read_text(encoding='utf-8')
    assert readme.startswith('# WaferScope')
    assert 'https://github.com/wufansaber/waferscope' in readme
    guide = (root/'docs/PROJECT_MIGRATION_CN.md').read_text(encoding='utf-8')
    assert 'waferscope-2026-10-21' in guide
    assert 'WaferScope_Migration_20261021.zip' in guide
