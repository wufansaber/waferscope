"""Reusable Qt controls for the SCATMECH simulation page."""

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from PyQt5.QtCore import QThread, pyqtSlot
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
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .backend import SCATMECH_IMPORT_ERROR
from .materials import PARTICLE_MATERIALS, WAFER_MATERIALS
from .parameters import SimulationParameters
from .plotting import PlotCanvas, render_results
from .worker import SimulationWorker


class MaterialGroup(QGroupBox):
    def __init__(self, title, materials, include_thickness=False, include_size=False):
        super().__init__(title)
        self.materials = materials
        layout = QFormLayout(self)
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


class ScatmechSimulationWidget(QWidget):
    """Embeddable SCATMECH simulation page used by the APTK main window."""

    def __init__(self):
        super().__init__()
        self.result = None
        self.thread = None
        self.worker = None
        self._build_ui()

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(6)

        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setMinimumWidth(285)
        controls_scroll.setMaximumWidth(325)
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
        controls_layout.addWidget(self.wafer_group)
        controls_layout.addWidget(self.particle_group)

        laser_group = QGroupBox("Laser and incidence")
        laser_layout = QFormLayout(laser_group)
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
        self.polarization = QComboBox()
        self.polarization.addItem("P polarized", "p")
        self.polarization.addItem("S polarized", "s")
        self.polarization.addItem("Unpolarized", "unpolarized")
        laser_layout.addRow("Wavelength", self.wavelength)
        laser_layout.addRow("Incidence angle", self.incidence)
        laser_layout.addRow("Polarization", self.polarization)
        controls_layout.addWidget(laser_group)

        numerical_group = QGroupBox("Numerical and display settings")
        numerical_layout = QFormLayout(numerical_group)
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

        self.run_button = QPushButton("Run far-field simulation")
        self.run_button.setMinimumHeight(32)
        self.run_button.clicked.connect(self.start_simulation)
        self.export_button = QPushButton("Export figures and data")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self.export_results)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.status = QLabel("Ready")
        self.status.setWordWrap(True)
        controls_layout.addWidget(self.run_button)
        controls_layout.addWidget(self.export_button)
        controls_layout.addWidget(self.progress)
        controls_layout.addWidget(self.status)
        controls_layout.addStretch()
        controls_scroll.setWidget(controls)
        root.addWidget(controls_scroll)

        self.tabs = QTabWidget()
        self.canvases = {}
        for key, title in (
            ("3d_intensity", "3-D intensity + polarization"),
            ("2d_intensity", "2-D intensity + polarization"),
        ):
            canvas = PlotCanvas()
            self.canvases[key] = canvas
            self.tabs.addTab(canvas, title)
        root.addWidget(self.tabs, stretch=1)

        self.setStyleSheet(
            "QGroupBox { font-weight: bold; margin-top: 6px; }"
            "QGroupBox::title { subcontrol-origin: margin; left: 6px; }"
            "QPushButton { padding: 5px; }"
            "QLabel, QCheckBox { font-size: 11px; }"
        )
        if SCATMECH_IMPORT_ERROR is not None:
            self.run_button.setEnabled(False)
            self.status.setText(
                "pySCATMECH is not installed in the active interpreter. "
                "The integrated UI is available, but simulation is disabled."
            )

    def parameters(self):
        return SimulationParameters(
            wafer_material=self.wafer_group.material_combo.currentText(),
            wafer_n=float(self.wafer_group.n_edit.text()),
            wafer_k=float(self.wafer_group.k_edit.text()),
            wafer_thickness_um=self.wafer_group.thickness.value(),
            particle_material=self.particle_group.material_combo.currentText(),
            particle_n=float(self.particle_group.n_edit.text()),
            particle_k=float(self.particle_group.k_edit.text()),
            particle_diameter_nm=self.particle_group.diameter.value(),
            wavelength_nm=self.wavelength.value(),
            incidence_deg=self.incidence.value(),
            polarization=self.polarization.currentData(),
            angular_step_deg=self.angular_step.value(),
            display_floor_db=-self.dynamic_range.value(),
        )

    @pyqtSlot()
    def start_simulation(self):
        try:
            parameters = self.parameters()
            parameters.validate()
        except Exception as error:
            QMessageBox.warning(self, "Invalid parameters", str(error))
            return

        self.run_button.setEnabled(False)
        self.export_button.setEnabled(False)
        self.progress.setValue(0)
        self.status.setText("Building SCATMECH model...")
        self.thread = QThread(self)
        self.worker = SimulationWorker(parameters)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.progress.setValue)
        self.worker.completed.connect(self._simulation_completed)
        self.worker.failed.connect(self._simulation_failed)
        self.worker.completed.connect(self.thread.quit)
        self.worker.failed.connect(self.thread.quit)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self._thread_finished)
        self.thread.start()

    @pyqtSlot()
    def _thread_finished(self):
        self.worker = None
        self.thread = None

    @pyqtSlot(object)
    def _simulation_completed(self, result):
        self.result = result
        self.render_results(result)
        self.run_button.setEnabled(True)
        self.export_button.setEnabled(True)
        self.progress.setValue(100)
        self.status.setText(
            "Completed. DSC range: {:.3e} to {:.3e} um^2/sr\n{}".format(
                np.min(result["intensity"]),
                np.max(result["intensity"]),
                result["thickness_note"],
            )
        )

    @pyqtSlot(str)
    def _simulation_failed(self, details):
        self.run_button.setEnabled(True)
        self.progress.setValue(0)
        self.status.setText("Simulation failed.")
        QMessageBox.critical(self, "Simulation failed", details)

    def render_results(self, result):
        render_results(self.canvases, result)

    def export_results(self):
        if self.result is None:
            return
        directory = QFileDialog.getExistingDirectory(self, "Select export directory")
        if not directory:
            return
        directory = Path(directory)
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
            json.dump(asdict(self.result["parameters"]), stream, indent=2)
        QMessageBox.information(self, "Export complete", str(directory))

    def wait_for_shutdown(self):
        """Wait for an active native SCATMECH calculation before app exit."""
        if self.thread is not None and self.thread.isRunning():
            self.status.setText("Waiting for the active SCATMECH calculation...")
            self.thread.quit()
            self.thread.wait()
