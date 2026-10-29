"""Standalone application entry points for the SCATMECH interface."""

import json
import sys
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QTabWidget

from .backend import calculate_far_field
from .parameters import SimulationParameters
from .widgets import ScatmechSimulationWidget


class ScatmechStandaloneWindow(QMainWindow):
    """Standalone wrapper retained for direct execution."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("WaferScope")
        self.setMinimumSize(1200, 720)
        self._resize_for_screen()
        self.simulation_widget = ScatmechSimulationWidget()
        from fullwave_ui.widget import FullwaveWidget
        self.fullwave_widget = FullwaveWidget()
        from grating_ui.widget import GratingWidget
        self.grating_widget = GratingWidget()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.simulation_widget, 'SCATMECH')
        self.tabs.addTab(self.fullwave_widget, 'Full-wave FEM (experimental)')
        self.tabs.addTab(self.grating_widget, 'S4 Grating (periodic)')
        self.setCentralWidget(self.tabs)
        self.simulation_widget.roughness_group.match_fem_button.clicked.connect(
            self.use_fem_haze_inputs)
        self.simulation_widget.roughness_group.match_fem_button.show()

    def use_fem_haze_inputs(self):
        """Match a completed statistical-PSD FEM result, not editable form values."""
        target = self.simulation_widget
        source_dir = self.fullwave_widget.loaded_source
        result_dir = self.fullwave_widget.loaded_result
        if source_dir is None or result_dir is None:
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'Run or load a complete FEM result first. Form values may differ from a saved result.')
            return
        if target.coating_group.layers():
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'Remove SCATMECH coating layers first. FEM currently uses a homogeneous wafer.')
            return
        try:
            report = json.loads((source_dir/'report.json').read_text(encoding='utf-8'))
            image_report = json.loads((result_dir/'imaging_report.json').read_text(encoding='utf-8'))
            physical = report['parameters']
            roughness = report['roughness']
            imaging = image_report['settings']
        except (OSError, ValueError, KeyError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'The displayed FEM result is incomplete: {}'.format(error))
            return
        if roughness['model'] not in ('gaussian', 'abc'):
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'The displayed FEM result uses {} roughness without a matching '
                'stationary PSD. Select Gaussian or ABC on the FEM page and '
                'calculate a new result.'.format(roughness['model']))
            return
        if (roughness['rms_nm'] <= 0 or physical['beam_waist_um'] != 0 or
                physical.get('beam_profile', 'plane') != 'plane'):
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'This comparison requires nonzero roughness and plane-wave FEM illumination.')
            return
        irradiance = 1000*imaging['irradiance_w_um2']
        irradiance_control = target.imaging_group.irradiance
        if not irradiance_control.minimum() <= irradiance <= irradiance_control.maximum():
            QMessageBox.warning(self, 'Cannot match haze inputs',
                'FEM beam-normal irradiance is outside the SCATMECH input range.')
            return
        target.wavelength.setValue(physical['wavelength_nm'])
        target.incidence.setValue(physical['incidence_deg'])
        target.polarization.setCurrentIndex(target.polarization.findData(
            physical['polarization']))
        target.wafer_group.n_edit.setText(str(physical['substrate_n']))
        target.wafer_group.k_edit.setText(str(physical['substrate_k']))
        rough = target.roughness_group
        rough.enabled.setChecked(True)
        rough.haze_only.setChecked(True)
        rough.psd.setCurrentIndex(rough.psd.findData(roughness['model']))
        rough.mode.setCurrentIndex(rough.mode.findData('top_interface'))
        rough.rms_nm.setValue(roughness['rms_nm'])
        rough.correlation_length_um.setValue(roughness['correlation_length_um'])
        rough.exponent.setValue(roughness.get('psd_exponent', 3.))
        patch_width = roughness['width_um']
        rough.footprint_length_um.setValue(patch_width)
        rough.footprint_width_um.setValue(patch_width)
        target.angular_step.setValue(report['angular_step_deg'])
        target.objective_focal_length.setValue(imaging['focal_length_mm'])
        target.imaging_group.magnification.setValue(imaging['magnification'])
        for index in (1, 2, 3):
            getattr(target, 'objective_na_ch'+str(index)).setValue(
                imaging['numerical_apertures'][index-1])
            getattr(target, 'objective_tilt_ch'+str(index)).setValue(
                imaging['tilt_degrees'][index-1])
        irradiance_control.setValue(irradiance)
        target.result = None
        target.far_field_result = None
        target.export_button.setEnabled(False)
        for key in ('3d_intensity', '2d_intensity'):
            canvas = target.canvases[key]
            canvas.figure.clear()
            canvas.draw_idle()
        target.status.setText(
            'Matched displayed {} FEM result ({}) on nominal inputs and '
            'rough area {:.6g} um^2. Run SCATMECH far-field now. Statistical ensemble '
            'PSD and one FEM height realization can still differ in angular shape.'.format(
                roughness['model'].upper(), result_dir.name, patch_width**2))

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
        if not self.grating_widget.wait_for_shutdown():
            event.ignore()
            return
        if not self.fullwave_widget.wait_for_shutdown():
            event.ignore()
            return
        self.simulation_widget.wait_for_shutdown()
        super().closeEvent(event)


def run_headless_test(output_directory):
    parameters = SimulationParameters(angular_step_deg=12.0)
    from .surface_model import build_scattering_model
    sphere_gap_um = float(build_scattering_model(parameters).model.getParameters()['delta'])
    if sphere_gap_um != 0:
        raise RuntimeError('SCATMECH sphere must contact the top wafer interface')
    result = calculate_far_field(parameters)
    app = QApplication.instance() or QApplication([])
    window = ScatmechSimulationWidget()
    window.result = result
    window.render_results(result)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    from result_archive import save_result_archive, load_result_archive
    from .parameters import ObjectiveCollectionParameters
    from .objective_collection import calculate_objective_collection_example
    native_path = output_directory/'native_far_field_data.npz'
    save_result_archive(native_path, result, legacy_arrays={
        key: value for key, value in result.items() if hasattr(value, 'dtype')})
    collected = calculate_objective_collection_example(ObjectiveCollectionParameters(), result)
    collection_path = output_directory/'far_field_data.npz'
    save_result_archive(collection_path, collected, source=result, legacy_arrays={
        key: value for key, value in collected.items() if hasattr(value, 'dtype')})
    restored, source = load_result_archive(collection_path)
    if len(restored['channel_pupils']) != 3 or 'mueller_raw' not in source:
        raise RuntimeError('Packaged raw result export is incomplete')
    from .imaging import ImagingParameters, calculate_imaging, export_imaging
    from matplotlib.figure import Figure
    image = calculate_imaging(SimulationParameters(wavelength_nm=266, particle_diameter_nm=10),
        ObjectiveCollectionParameters(channel_numerical_apertures=(.4, .4, 0)),
        ImagingParameters(pupil_samples=32))
    if any(channel['field_reference'] != 'top-wafer plane z=0' for channel in image['channels']):
        raise RuntimeError('Packaged coherent field reference is incorrect')
    export_imaging(image, output_directory/'coherent', Figure())
    (output_directory/'raw_export_smoke.json').write_text(json.dumps(dict(
        status='passed', channels=3, native_mueller=True, source_snapshot=True,
        sphere_gap_um=sphere_gap_um, coherent_wafer_reference=True)), encoding='utf-8')
    for key, canvas in window.canvases.items():
        canvas.figure.savefig(str(output_directory / (key + ".png")), dpi=120)
    print("Headless UI render test completed:", output_directory.resolve())


def main():
    if '--grating-package-test' in sys.argv:
        output = Path(sys.argv[sys.argv.index('--grating-package-test')+1])
        try:
            from grating_ui.package_smoke import run
            run(output)
        except Exception:
            import traceback
            output.mkdir(parents=True, exist_ok=True)
            (output/'grating_package_error.txt').write_text(traceback.format_exc(), encoding='utf-8')
            raise SystemExit(1)
        return
    if '--camera-package-test' in sys.argv:
        from fullwave_ui.camera_smoke import run
        index = sys.argv.index('--camera-package-test')
        run(sys.argv[index+1], sys.argv[index+2])
        return
    if '--haze-package-test' in sys.argv:
        from fullwave_ui.ensemble_smoke import run
        run(sys.argv[sys.argv.index('--haze-package-test')+1])
        return
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
    app.setApplicationName('WaferScope')
    app.setApplicationDisplayName('WaferScope')
    window = ScatmechStandaloneWindow()
    window.show()
    sys.exit(app.exec_())
