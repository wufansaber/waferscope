"""Reusable Qt controls for the SCATMECH simulation page."""

import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from PyQt5.QtCore import Qt, QThread, pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .backend import SCATMECH_IMPORT_ERROR, calculate_objective_collection_example
from .coordinate_geometry import (
    incidence_metadata_from_result,
    incidence_projection_text,
)
from .materials import PARTICLE_MATERIALS, WAFER_MATERIALS
from .parameters import (
    ObjectiveCollectionParameters,
    ParticleSizeSweepParameters,
    SimulationParameters,
)
from .plotting import DetectorCanvas, ImagingCanvas, PlotCanvas, render_results
from .size_sweep import export_particle_size_sweep
from .wafer_controls import CoatingStackGroup, RoughnessGroup
from .worker import ParticleSizeSweepWorker, SimulationWorker
from .worker import ImagingWorker
from .imaging_controls import ImagingGroup
from .imaging import export_imaging
from .detector_controls import DetectorGroup


class MaterialGroup(QGroupBox):
    def __init__(self, title, materials, include_thickness=False, include_size=False):
        super().__init__(title)
        self.materials = materials
        layout = QFormLayout(self)
        layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        layout.setContentsMargins(8, 10, 8, 6)
        layout.setHorizontalSpacing(6)
        layout.setVerticalSpacing(5)
        self.material_combo = QComboBox()
        self.material_combo.addItems(materials.keys())
        self.n_edit = QLineEdit()
        self.k_edit = QLineEdit()
        validator = QDoubleValidator(0.0, 100.0, 8, self)
        self.n_edit.setValidator(validator)
        self.k_edit.setValidator(validator)
        layout.addRow("Material", self.material_combo)
        layout.addRow("Refractive index n", self.n_edit)
        layout.addRow("Extinction coefficient k", self.k_edit)
        if include_thickness:
            self.thickness = QDoubleSpinBox()
            self.thickness.setRange(0.001, 5000.0)
            self.thickness.setDecimals(6)
            self.thickness.setValue(725.0)
            self.thickness.setSuffix(" um")
            layout.addRow("Thickness", self.thickness)
        if include_size:
            self.diameter = QDoubleSpinBox()
            self.diameter.setRange(1.0, 10000.0)
            self.diameter.setDecimals(3)
            self.diameter.setValue(100.0)
            self.diameter.setSuffix(" nm")
            layout.addRow("Sphere diameter", self.diameter)
        self.material_combo.currentTextChanged.connect(self._apply_material)
        self._apply_material(self.material_combo.currentText())

    def _apply_material(self, name):
        value = self.materials[name]
        self.n_edit.setText("{:.8g}".format(value.real))
        self.k_edit.setText("{:.8g}".format(value.imag))


class ControlsScrollArea(QScrollArea):
    """Keep the controls panel fitted to the resizable sidebar viewport."""

    def resizeEvent(self, event):
        super().resizeEvent(event)
        panel = self.widget()
        if panel is not None:
            panel.setFixedWidth(max(0, self.viewport().width()))


class ScatmechSimulationWidget(QWidget):
    """Embeddable SCATMECH simulation page used by the APTK main window."""

    def __init__(self):
        super().__init__()
        self.result = None
        self.far_field_result = None
        self.thread = None
        self.worker = None
        self._build_ui()

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(6)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(5)
        root.addWidget(splitter)

        controls_scroll = ControlsScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        controls_scroll.setMinimumWidth(430)
        controls_scroll.setMaximumWidth(680)
        controls = QWidget()
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(5, 5, 5, 5)
        controls_layout.setSpacing(6)

        self.wafer_group = MaterialGroup(
            "Wafer", WAFER_MATERIALS, include_thickness=True
        )
        self.particle_group = MaterialGroup(
            "Spherical particle", PARTICLE_MATERIALS, include_size=True
        )
        self.coating_group = CoatingStackGroup()
        self.roughness_group = RoughnessGroup()
        controls_layout.addWidget(self.wafer_group)
        controls_layout.addWidget(self.coating_group)
        controls_layout.addWidget(self.roughness_group)
        controls_layout.addWidget(self.particle_group)

        laser_group = QGroupBox("Laser and incidence")
        laser_layout = QFormLayout(laser_group)
        laser_layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        laser_layout.setContentsMargins(8, 10, 8, 6)
        laser_layout.setHorizontalSpacing(6)
        laser_layout.setVerticalSpacing(5)
        self.wavelength = QDoubleSpinBox()
        self.wavelength.setRange(100.0, 2000.0)
        self.wavelength.setDecimals(3)
        self.wavelength.setValue(355.0)
        self.wavelength.setSuffix(" nm")
        self.incidence = QDoubleSpinBox()
        self.incidence.setRange(0.0, 89.9)
        self.incidence.setDecimals(3)
        self.incidence.setValue(75.0)
        self.incidence.setSuffix(" deg")
        self.incidence_projection = QLabel(
            incidence_projection_text(self.incidence.value())
        )
        self.incidence_projection.setWordWrap(True)
        self.incidence_projection.setToolTip(
            "The laser source is on the -KY side. The incident propagation "
            "vector lies in the KY-KZ plane and projects toward +KY."
        )
        self.incidence.valueChanged.connect(self._update_incidence_projection_label)
        self.polarization = QComboBox()
        self.polarization.addItem("P polarized", "p")
        self.polarization.addItem("S polarized", "s")
        self.polarization.addItem("Unpolarized", "unpolarized")
        laser_layout.addRow("Wavelength", self.wavelength)
        laser_layout.addRow("Incidence angle", self.incidence)
        laser_layout.addRow("In-plane projection", self.incidence_projection)
        laser_layout.addRow("Polarization", self.polarization)
        controls_layout.addWidget(laser_group)

        numerical_group = QGroupBox("Numerical and display settings")
        numerical_layout = QFormLayout(numerical_group)
        numerical_layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        numerical_layout.setContentsMargins(8, 10, 8, 6)
        numerical_layout.setHorizontalSpacing(6)
        numerical_layout.setVerticalSpacing(5)
        self.angular_step = QDoubleSpinBox()
        self.angular_step.setRange(2.0, 15.0)
        self.angular_step.setDecimals(1)
        self.angular_step.setSingleStep(1.0)
        self.angular_step.setValue(4.0)
        self.angular_step.setSuffix(" deg")
        self.dynamic_range = QDoubleSpinBox()
        self.dynamic_range.setRange(10.0, 80.0)
        self.dynamic_range.setValue(35.0)
        self.dynamic_range.setSuffix(" dB")
        numerical_layout.addRow("Angular step", self.angular_step)
        numerical_layout.addRow("Display range", self.dynamic_range)
        controls_layout.addWidget(numerical_group)

        objective_group = QGroupBox("Objective collection example")
        objective_layout = QFormLayout(objective_group)
        objective_layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        objective_layout.setContentsMargins(8, 10, 8, 6)
        objective_layout.setHorizontalSpacing(6)
        objective_layout.setVerticalSpacing(5)
        self.objective_focal_length = QDoubleSpinBox()
        self.objective_focal_length.setRange(0.1, 1000.0)
        self.objective_focal_length.setDecimals(3)
        self.objective_focal_length.setValue(20.0)
        self.objective_focal_length.setSuffix(" mm")
        self.objective_tilt_ch1 = QDoubleSpinBox()
        self.objective_tilt_ch1.setRange(-89.9, 89.9)
        self.objective_tilt_ch1.setDecimals(3)
        self.objective_tilt_ch1.setSingleStep(5.0)
        self.objective_tilt_ch1.setValue(45.0)
        self.objective_tilt_ch1.setSuffix(" deg")
        self.objective_na_ch1 = self._create_objective_na_spinbox()
        self.objective_tilt_ch2 = QDoubleSpinBox()
        self.objective_tilt_ch2.setRange(-89.9, 89.9)
        self.objective_tilt_ch2.setDecimals(3)
        self.objective_tilt_ch2.setSingleStep(5.0)
        self.objective_tilt_ch2.setValue(-45.0)
        self.objective_tilt_ch2.setSuffix(" deg")
        self.objective_na_ch2 = self._create_objective_na_spinbox()
        self.objective_tilt_ch3 = QDoubleSpinBox()
        self.objective_tilt_ch3.setRange(-89.9, 89.9)
        self.objective_tilt_ch3.setDecimals(3)
        self.objective_tilt_ch3.setSingleStep(5.0)
        self.objective_tilt_ch3.setValue(0.0)
        self.objective_tilt_ch3.setSuffix(" deg")
        self.objective_na_ch3 = self._create_objective_na_spinbox()
        objective_layout.addRow("Focal length", self.objective_focal_length)
        objective_layout.addRow("Channel 1 NA", self.objective_na_ch1)
        objective_layout.addRow("Channel 1 tilt", self.objective_tilt_ch1)
        objective_layout.addRow("Channel 2 NA", self.objective_na_ch2)
        objective_layout.addRow("Channel 2 tilt", self.objective_tilt_ch2)
        objective_layout.addRow("Channel 3 NA", self.objective_na_ch3)
        objective_layout.addRow("Channel 3 tilt", self.objective_tilt_ch3)
        controls_layout.addWidget(objective_group)

        self.imaging_group = ImagingGroup(self.objective_focal_length)
        self.imaging_group.run_button.clicked.connect(self.start_imaging)
        controls_layout.addWidget(self.imaging_group)
        self.detector_group = DetectorGroup()
        controls_layout.addWidget(self.detector_group)

        sweep_group = QGroupBox("Particle-size sweep")
        sweep_layout = QFormLayout(sweep_group)
        sweep_layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        sweep_layout.setContentsMargins(8, 10, 8, 6)
        sweep_layout.setHorizontalSpacing(6)
        sweep_layout.setVerticalSpacing(5)
        self.sweep_start_diameter = self._create_sweep_spinbox(30.0)
        self.sweep_end_diameter = self._create_sweep_spinbox(200.0)
        self.sweep_step_diameter = self._create_sweep_spinbox(10.0)
        sweep_layout.addRow("Start diameter", self.sweep_start_diameter)
        sweep_layout.addRow("End diameter", self.sweep_end_diameter)
        sweep_layout.addRow("Diameter step", self.sweep_step_diameter)
        self.sweep_button = QPushButton("Run particle-size sweep")
        self.sweep_button.setMinimumHeight(30)
        self.sweep_button.clicked.connect(self.start_particle_size_sweep)
        sweep_layout.addRow(self.sweep_button)
        controls_layout.addWidget(sweep_group)

        self.run_button = QPushButton("Run far-field simulation")
        self.run_button.setMinimumHeight(32)
        self.run_button.clicked.connect(self.start_simulation)
        self.objective_example_button = QPushButton("Run ideal objective example")
        self.objective_example_button.setMinimumHeight(30)
        self.objective_example_button.clicked.connect(
            self.run_objective_collection_example
        )
        self.export_button = QPushButton("Export figures and data")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self.export_results)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.status = QLabel("Ready")
        self.status.setWordWrap(True)
        controls_layout.addWidget(self.run_button)
        controls_layout.addWidget(self.objective_example_button)
        controls_layout.addWidget(self.export_button)
        controls_layout.addWidget(self.progress)
        controls_layout.addWidget(self.status)
        controls_layout.addStretch()
        controls_scroll.setWidget(controls)
        splitter.addWidget(controls_scroll)

        self.tabs = QTabWidget()
        self.canvases = {}
        for key, title in (
            ("3d_intensity", "3-D far-field / pupil"),
            ("2d_intensity", "2-D far-field / pupil"),
            ("size_sweep", "Particle-size sweep"),
            ("image_planes", "Three image planes"),
            ("detectors", "Detectors"),
        ):
            canvas_type = {"image_planes": ImagingCanvas, "detectors": DetectorCanvas}.get(key, PlotCanvas)
            canvas = canvas_type()
            self.canvases[key] = canvas
            self.tabs.addTab(canvas, title)
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes((540, 1100))
        self.controls_scroll = controls_scroll
        self.controls_panel = controls
        self.main_splitter = splitter

        self.setStyleSheet(
            "QGroupBox { font-weight: bold; margin-top: 6px; }"
            "QGroupBox::title { subcontrol-origin: margin; left: 6px; }"
            "QPushButton { padding: 5px; }"
            "QLabel, QCheckBox { font-size: 11px; }"
        )
        if SCATMECH_IMPORT_ERROR is not None:
            self.run_button.setEnabled(False)
            self.sweep_button.setEnabled(False)
            self.imaging_group.run_button.setEnabled(False)
            self.status.setText(
                "pySCATMECH is not installed in the active interpreter. "
                "The native particle simulation is disabled, but the ideal "
                "objective example is available."
            )

    def _create_objective_na_spinbox(self):
        spinbox = QDoubleSpinBox()
        spinbox.setRange(0.0, 0.99)
        spinbox.setDecimals(3)
        spinbox.setSingleStep(0.05)
        spinbox.setValue(0.5)
        spinbox.setToolTip("Set NA to 0 to disable this collection channel.")
        return spinbox

    def _update_incidence_projection_label(self, incidence_deg):
        self.incidence_projection.setText(incidence_projection_text(incidence_deg))

    def _create_sweep_spinbox(self, value):
        spinbox = QDoubleSpinBox()
        spinbox.setRange(0.1, 10000.0)
        spinbox.setDecimals(3)
        spinbox.setSingleStep(10.0)
        spinbox.setValue(value)
        spinbox.setSuffix(" nm")
        return spinbox

    def parameters(self):
        roughness_values = self.roughness_group.values()
        return SimulationParameters(
            wafer_material=self.wafer_group.material_combo.currentText(),
            wafer_n=float(self.wafer_group.n_edit.text()),
            wafer_k=float(self.wafer_group.k_edit.text()),
            wafer_thickness_um=self.wafer_group.thickness.value(),
            coating_layers=self.coating_group.layers(),
            particle_material=self.particle_group.material_combo.currentText(),
            particle_n=float(self.particle_group.n_edit.text()),
            particle_k=float(self.particle_group.k_edit.text()),
            particle_diameter_nm=self.particle_group.diameter.value(),
            wavelength_nm=self.wavelength.value(),
            incidence_deg=self.incidence.value(),
            polarization=self.polarization.currentData(),
            angular_step_deg=self.angular_step.value(),
            display_floor_db=-self.dynamic_range.value(),
            **roughness_values,
        )

    def objective_collection_parameters(self):
        return ObjectiveCollectionParameters(
            channel_numerical_apertures=(
                self.objective_na_ch1.value(),
                self.objective_na_ch2.value(),
                self.objective_na_ch3.value(),
            ),
            focal_length_mm=self.objective_focal_length.value(),
            tilt_degrees=(
                self.objective_tilt_ch1.value(),
                self.objective_tilt_ch2.value(),
                self.objective_tilt_ch3.value(),
            ),
            display_floor_db=-self.dynamic_range.value(),
        )

    def particle_size_sweep_parameters(self):
        return ParticleSizeSweepParameters(
            start_diameter_nm=self.sweep_start_diameter.value(),
            end_diameter_nm=self.sweep_end_diameter.value(),
            step_diameter_nm=self.sweep_step_diameter.value(),
        )

    @pyqtSlot()
    def start_simulation(self):
        try:
            parameters = self.parameters()
            parameters.validate()
        except Exception as error:
            QMessageBox.warning(self, "Invalid parameters", str(error))
            return

        self._start_worker(
            SimulationWorker(parameters),
            self._simulation_completed,
            "Building SCATMECH model...",
        )

    @pyqtSlot()
    def start_particle_size_sweep(self):
        try:
            simulation_parameters = self.parameters()
            objective_parameters = self.objective_collection_parameters()
            sweep_parameters = self.particle_size_sweep_parameters()
            simulation_parameters.validate()
            objective_parameters.validate()
            sweep_parameters.validate()
        except Exception as error:
            QMessageBox.warning(self, "Invalid sweep parameters", str(error))
            return

        self._start_worker(
            ParticleSizeSweepWorker(
                simulation_parameters,
                objective_parameters,
                sweep_parameters,
            ),
            self._sweep_completed,
            "Running particle-size sweep...",
        )

    @pyqtSlot()
    def start_imaging(self):
        try:
            simulation = self.parameters()
            objective = self.objective_collection_parameters()
            imaging = replace(self.imaging_group.parameters(), detector=self.detector_group.parameters())
            simulation.validate()
            objective.validate()
            imaging.validate()
            if simulation.roughness_enabled:
                raise ValueError("Disable roughness for coherent particle imaging. "
                                 "Statistical BRDF alone does not define a complex field.")
        except Exception as error:
            QMessageBox.warning(self, "Invalid imaging parameters", str(error))
            return
        self._start_worker(ImagingWorker(simulation, objective, imaging),
                           self._imaging_completed, "Sampling complex Jones pupils...")

    @pyqtSlot(object)
    def _imaging_completed(self, result):
        self.result = result
        self.render_results(result)
        self.export_button.setEnabled(True)
        self.progress.setValue(100)
        powers = "/".join("{:.4g}".format(c["power_w"]) for c in result["channels"])
        warning = ""
        if any(c["admits_specular"] for c in result["channels"]):
            warning += " WARNING: a collection cone admits specular reflection; a dark-field stop is required."
        if any(c["edge_fraction"] > 0.02 for c in result["channels"]):
            warning += " Image edge energy is high; increase pupil samples to enlarge the field."
        if result["imaging_parameters"].detector.enabled:
            counts = "/".join("{:.4g}".format(c["detector"]["signal_e"]) for c in result["channels"])
            warning += " Expected particle signal: {} e-.".format(counts)
            if any(c["detector"]["saturated_pixels"] or c["detector"]["expected_overfull_pixels"]
                   for c in result["channels"]):
                warning += " WARNING: detector saturation; reported SNR is pre-saturation only."
        self.status.setText("Particle-only image powers: {} W. Ideal aplanatic optics; "
                            "near-field and background leakage are not included.{}".format(powers, warning))

    def _start_worker(self, worker, completed_slot, status_text):
        self._set_calculation_controls_enabled(False)
        self.export_button.setEnabled(False)
        self.progress.setValue(0)
        self.status.setText(status_text)
        self.thread = QThread(self)
        self.worker = worker
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.progress.setValue)
        self.worker.completed.connect(completed_slot)
        self.worker.failed.connect(self._simulation_failed)
        self.worker.completed.connect(self.thread.quit)
        self.worker.failed.connect(self.thread.quit)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self._thread_finished)
        self.thread.start()

    def _set_calculation_controls_enabled(self, enabled):
        native_enabled = enabled and SCATMECH_IMPORT_ERROR is None
        self.run_button.setEnabled(native_enabled)
        self.sweep_button.setEnabled(native_enabled)
        self.imaging_group.run_button.setEnabled(native_enabled)
        self.objective_example_button.setEnabled(enabled)

    @pyqtSlot()
    def _thread_finished(self):
        self.worker = None
        self.thread = None
        self._set_calculation_controls_enabled(True)

    @pyqtSlot(object)
    def _simulation_completed(self, result):
        self.result = result
        self.far_field_result = result
        self.render_results(result)
        self.export_button.setEnabled(True)
        self.progress.setValue(100)
        self.status.setText(
            "Completed. Effective DSC range: {:.3e} to {:.3e} um^2/sr\n{}\n{}".format(
                np.min(result["intensity"]),
                np.max(result["intensity"]),
                result["thickness_note"],
                result["model_note"],
            )
        )

    @pyqtSlot(str)
    def _simulation_failed(self, details):
        self.progress.setValue(0)
        self.status.setText("Simulation failed.")
        QMessageBox.critical(self, "Simulation failed", details)

    @pyqtSlot(object)
    def _sweep_completed(self, result):
        self.result = result
        self.render_results(result)
        self.export_button.setEnabled(True)
        self.progress.setValue(100)
        values = result["channel_collected_dsc_um2"]
        self.status.setText(
            "Particle-size sweep completed: {} points, {:.3g} to {:.3g} nm. "
            "Collected effective DSC range: {:.3e} to {:.3e} um^2.\n{}".format(
                len(result["diameters_nm"]),
                result["diameters_nm"][0],
                result["diameters_nm"][-1],
                np.min(values),
                np.max(values),
                result["model_note"],
            )
        )

    def run_objective_collection_example(self):
        if self.thread is not None and self.thread.isRunning():
            QMessageBox.information(
                self,
                "Simulation is running",
                "Wait for the active far-field simulation to finish first.",
            )
            return

        try:
            parameters = self.objective_collection_parameters()
            result = calculate_objective_collection_example(
                parameters,
                self.far_field_result,
            )
        except Exception as error:
            QMessageBox.warning(self, "Invalid objective example", str(error))
            return

        self.result = result
        self.render_results(result)
        self.export_button.setEnabled(True)
        self.progress.setValue(100)
        if self.far_field_result is None:
            source_status = "No SCATMECH far-field result is available; using the synthetic dark-field example."
        else:
            source_status = "Using the latest SCATMECH far-field result."
        self.status.setText(
            "Objective collection completed. {}\n{}".format(
                source_status,
                result["thickness_note"],
            )
        )

    def render_results(self, result):
        render_results(self.canvases, result)
        if result.get("plot_kind") == "coherent_imaging":
            self.tabs.setCurrentWidget(self.canvases["image_planes"])
        elif result.get("plot_kind") == "particle_size_sweep":
            self.tabs.setCurrentWidget(self.canvases["size_sweep"])
        else:
            self.tabs.setCurrentWidget(self.canvases["3d_intensity"])

    def export_results(self):
        if self.result is None:
            return
        directory = QFileDialog.getExistingDirectory(self, "Select export directory")
        if not directory:
            return
        directory = Path(directory)
        if self.result.get("plot_kind") == "coherent_imaging":
            export_imaging(self.result, directory, self.canvases["image_planes"].figure,
                           self.canvases["detectors"].figure)
            QMessageBox.information(self, "Export complete", str(directory))
            return
        if self.result.get("plot_kind") == "particle_size_sweep":
            export_particle_size_sweep(
                self.result,
                directory,
                self.canvases["size_sweep"].figure,
            )
            QMessageBox.information(self, "Export complete", str(directory))
            return
        names = {
            "3d_intensity": "far_field_3d_intensity_with_polarization.png",
            "2d_intensity": "far_field_2d_intensity_with_polarization.png",
        }
        for key, filename in names.items():
            self.canvases[key].figure.savefig(
                str(directory / filename), dpi=180, facecolor="white"
            )
        arrays = {
            key: value
            for key, value in self.result.items()
            if isinstance(value, np.ndarray)
        }
        np.savez_compressed(str(directory / "far_field_data.npz"), **arrays)
        with (directory / "simulation_parameters.json").open(
            "w", encoding="utf-8"
        ) as stream:
            payload = asdict(self.result["parameters"])
            if hasattr(self.result["parameters"], "illuminated_area_um2"):
                payload["illuminated_area_um2"] = self.result[
                    "parameters"
                ].illuminated_area_um2
            payload["wafer_stack_note"] = self.result.get("thickness_note", "")
            payload["scatmech_model_note"] = self.result.get("model_note", "")
            coordinate_convention = incidence_metadata_from_result(self.result)
            if coordinate_convention:
                payload["coordinate_convention"] = coordinate_convention
            json.dump(payload, stream, indent=2)
        QMessageBox.information(self, "Export complete", str(directory))

    def wait_for_shutdown(self):
        """Wait for an active native SCATMECH calculation before app exit."""
        if self.thread is not None and self.thread.isRunning():
            self.status.setText("Waiting for the active SCATMECH calculation...")
            self.thread.quit()
            self.thread.wait()
