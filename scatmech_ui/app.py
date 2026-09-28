"""Standalone application entry points for the SCATMECH interface."""

import sys
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget

from .backend import calculate_far_field
from .parameters import SimulationParameters
from .widgets import ScatmechSimulationWidget


class ScatmechStandaloneWindow(QMainWindow):
    """Standalone wrapper retained for direct execution."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("pySCATMECH Particle-on-Wafer Far-Field Simulator")
        self.setMinimumSize(1200, 720)
        self._resize_for_screen()
        self.simulation_widget = ScatmechSimulationWidget()
        from fullwave_ui.widget import FullwaveWidget
        self.fullwave_widget = FullwaveWidget()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.simulation_widget, 'SCATMECH')
        self.tabs.addTab(self.fullwave_widget, 'Full-wave FEM (experimental)')
        self.setCentralWidget(self.tabs)

    def _resize_for_screen(self):
        screen = QApplication.primaryScreen()
        if screen is None:
            self.resize(1680, 960)
            return
        available = screen.availableGeometry()
        width = min(1680, max(1200, int(available.width() * 0.92)))
        height = min(980, max(720, int(available.height() * 0.92)))
        self.resize(width, height)

    def closeEvent(self, event):
        if not self.fullwave_widget.wait_for_shutdown():
            event.ignore()
            return
        self.simulation_widget.wait_for_shutdown()
        super().closeEvent(event)


def run_headless_test(output_directory):
    parameters = SimulationParameters(angular_step_deg=12.0)
    result = calculate_far_field(parameters)
    app = QApplication.instance() or QApplication([])
    window = ScatmechSimulationWidget()
    window.result = result
    window.render_results(result)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    for key, canvas in window.canvases.items():
        canvas.figure.savefig(str(output_directory / (key + ".png")), dpi=120)
    print("Headless UI render test completed:", output_directory.resolve())


def main():
    if '--fem-package-test' in sys.argv:
        from fullwave_ui.package_smoke import run
        run(sys.argv[sys.argv.index('--fem-package-test')+1])
        return
    if "--headless-test" in sys.argv:
        index = sys.argv.index("--headless-test")
        output = (
            sys.argv[index + 1]
            if index + 1 < len(sys.argv)
            else "ui_test_results"
        )
        run_headless_test(output)
        return
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    window = ScatmechStandaloneWindow()
    window.show()
    sys.exit(app.exec_())
