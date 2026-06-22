"""PyQt5 front end for pySCATMECH particle-on-wafer far-field simulation."""

import json
import sys
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.cm as cm
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from PyQt5.QtCore import QObject, QThread, Qt, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
try:
    from pySCATMECH.fresnel import Film, FilmStack, OpticalFunction
    from pySCATMECH.local import Local_BRDF_Model
    from pySCATMECH.mueller import Polarization
    SCATMECH_IMPORT_ERROR = None
except ModuleNotFoundError as error:
    if error.name not in ("pySCATMECH", "SCATPY"):
        raise
    Film = FilmStack = OpticalFunction = Local_BRDF_Model = Polarization = None
    SCATMECH_IMPORT_ERROR = error


# Reference values near 355 nm. Every n and k field remains user-editable.
WAFER_MATERIALS = {
    "Silicon (Si)": 5.57 + 0.387j,
    "Fused silica (SiO2)": 1.476 + 0.0j,
    "Sapphire (Al2O3)": 1.79 + 0.0j,
    "Silicon carbide (SiC)": 2.75 + 0.02j,
    "Gallium arsenide (GaAs)": 3.65 + 1.80j,
    "Custom": 1.50 + 0.0j,
}

PARTICLE_MATERIALS = {
    "PSL": 1.59 + 0.0j,
    "Silica (SiO2)": 1.476 + 0.0j,
    "Silicon (Si)": 5.57 + 0.387j,
    "Alumina (Al2O3)": 1.79 + 0.0j,
    "Aluminum (Al)": 0.49 + 4.84j,
    "Gold (Au)": 1.48 + 1.86j,
    "Custom": 1.50 + 0.0j,
}

PARTICLE_INTENSITY_CMAP = LinearSegmentedColormap.from_list(
    "particle_intensity",
    [
        (0.00, "#050b18"),
        (0.08, "#151746"),
        (0.18, "#d7191c"),
        (0.48, "#ff6d00"),
        (0.82, "#fff21a"),
        (1.00, "#fffde8"),
    ],
)


@dataclass(frozen=True)
class SimulationParameters:
    wafer_material: str = "Silicon (Si)"
    wafer_n: float = 5.57
    wafer_k: float = 0.387
    wafer_thickness_um: float = 725.0
    particle_material: str = "PSL"
    particle_n: float = 1.59
    particle_k: float = 0.0
    particle_diameter_nm: float = 100.0
    wavelength_nm: float = 355.0
    incidence_deg: float = 75.0
    polarization: str = "p"
    angular_step_deg: float = 4.0
    display_floor_db: float = -35.0

    def validate(self):
        if self.wafer_n <= 0 or self.particle_n <= 0:
            raise ValueError("Refractive-index real parts must be positive.")
        if self.wafer_k < 0 or self.particle_k < 0:
            raise ValueError("Extinction coefficients cannot be negative.")
        if self.wafer_thickness_um <= 0:
            raise ValueError("Wafer thickness must be positive.")
        if self.particle_diameter_nm <= 0 or self.wavelength_nm <= 0:
            raise ValueError("Particle diameter and wavelength must be positive.")
        if not 0 <= self.incidence_deg < 90:
            raise ValueError("Incidence angle must be in [0, 90) degrees.")
        if not 2.0 <= self.angular_step_deg <= 15.0:
            raise ValueError("Angular step must be between 2 and 15 degrees.")


def _build_particle_model(parameters):
    """Build a finite-film model where numerically meaningful."""
    wavelength_um = parameters.wavelength_nm / 1000.0
    wafer_index = complex(parameters.wafer_n, parameters.wafer_k)
    particle_index = complex(parameters.particle_n, parameters.particle_k)
    thickness = parameters.wafer_thickness_um
    optical_depth = 4.0 * np.pi * parameters.wafer_k * thickness / wavelength_um
    coherent_limit = 100.0 * wavelength_um

    model_parameters = {
        "lambda": wavelength_um,
        "radius": parameters.particle_diameter_nm / 2000.0,
        "sphere": particle_index,
        "type": 0,
    }

    if thickness <= coherent_limit and optical_depth < 60.0:
        wafer_stack = FilmStack(
            [Film(OpticalFunction(wafer_index), thickness=thickness)]
        )
        model_parameters["substrate"] = 1.0 + 0.0j
        model_parameters["stack"] = wafer_stack
        thickness_note = (
            "Finite coherent wafer layer: {:.6g} um on air backing."
        ).format(thickness)
    else:
        model_parameters["substrate"] = wafer_index
        thickness_note = (
            "Wafer is optically/geometrically thick at this wavelength; "
            "the backside is negligible and the wafer is modeled as semi-infinite."
        )

    model = Local_BRDF_Model("Bobbert_Vlieger_BRDF_Model", model_parameters)
    return model, thickness_note


def calculate_far_field(parameters, progress_callback=None):
    """Return intensity and full Stokes-derived polarization on a hemisphere."""
    if SCATMECH_IMPORT_ERROR is not None:
        raise RuntimeError(
            "pySCATMECH is unavailable in the active Python interpreter. "
            "Select a supported environment containing pySCATMECH and SCATPY."
        ) from SCATMECH_IMPORT_ERROR
    parameters.validate()
    model, thickness_note = _build_particle_model(parameters)
    step = parameters.angular_step_deg
    theta_deg = np.arange(0.0, 90.0, step)
    if theta_deg[-1] < 89.0:
        theta_deg = np.append(theta_deg, 89.0)
    phi_count = max(24, int(round(360.0 / step)))
    phi = np.linspace(0.0, 2.0 * np.pi, phi_count + 1)
    theta = np.deg2rad(theta_deg)
    phi_grid, theta_grid = np.meshgrid(phi, theta)

    shape = theta_grid.shape
    intensity = np.empty(shape, dtype=float)
    q = np.empty(shape, dtype=float)
    u = np.empty(shape, dtype=float)
    v = np.empty(shape, dtype=float)
    incident = Polarization(
        {"p": "p", "s": "s", "unpolarized": "u"}[parameters.polarization]
    )
    incidence = np.deg2rad(parameters.incidence_deg)

    for row in range(shape[0]):
        for column in range(shape[1]):
            mueller = model.MuellerDSC(
                incidence,
                theta_grid[row, column],
                phi_grid[row, column],
                coords="psps",
            )
            output = mueller @ incident
            intensity[row, column] = max(float(output[0]), 0.0)
            q[row, column] = float(output[1])
            u[row, column] = float(output[2])
            v[row, column] = float(output[3])
        if progress_callback is not None:
            progress_callback(int(round(100.0 * (row + 1) / shape[0])))

    if not np.isfinite(intensity).all() or np.nanmax(intensity) <= 0:
        raise RuntimeError("SCATMECH returned invalid far-field intensity values.")

    sin_theta = np.sin(theta_grid)
    cos_theta = np.cos(theta_grid)
    cos_phi = np.cos(phi_grid)
    sin_phi = np.sin(phi_grid)
    x = sin_theta * cos_phi
    y = sin_theta * sin_phi
    z = cos_theta

    safe_intensity = np.maximum(intensity, np.nanmax(intensity) * 1e-15)
    relative = safe_intensity / np.nanmax(safe_intensity)
    intensity_db = 10.0 * np.log10(relative)
    dolp = np.clip(np.sqrt(q**2 + u**2) / safe_intensity, 0.0, 1.0)
    docp = np.clip(v / safe_intensity, -1.0, 1.0)

    # Q is Is-Ip in SCATMECH's s-p basis. alpha=0 is s; alpha=pi/2 is p.
    alpha = 0.5 * np.arctan2(u, q)
    s_x = -sin_phi
    s_y = cos_phi
    s_z = np.zeros_like(s_x)
    p_x = cos_theta * cos_phi
    p_y = cos_theta * sin_phi
    p_z = -sin_theta
    e_x = np.cos(alpha) * s_x + np.sin(alpha) * p_x
    e_y = np.cos(alpha) * s_y + np.sin(alpha) * p_y
    e_z = np.cos(alpha) * s_z + np.sin(alpha) * p_z

    return {
        "parameters": parameters,
        "thickness_note": thickness_note,
        "theta": theta_grid,
        "phi": phi_grid,
        "x": x,
        "y": y,
        "z": z,
        "intensity": intensity,
        "intensity_db": intensity_db,
        "q": q,
        "u": u,
        "v": v,
        "dolp": dolp,
        "docp": docp,
        "e_x": e_x,
        "e_y": e_y,
        "e_z": e_z,
    }


class SimulationWorker(QObject):
    progress = pyqtSignal(int)
    completed = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, parameters):
        super().__init__()
        self.parameters = parameters

    @pyqtSlot()
    def run(self):
        try:
            result = calculate_far_field(self.parameters, self.progress.emit)
            self.completed.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self):
        self.figure = Figure(figsize=(8, 6), constrained_layout=True)
        super().__init__(self.figure)


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

    def _sparse_indices(self, result):
        rows = np.arange(2, result["x"].shape[0] - 1, max(1, result["x"].shape[0] // 9))
        columns = np.arange(0, result["x"].shape[1] - 1, max(1, result["x"].shape[1] // 18))
        return np.ix_(rows, columns)

    def render_results(self, result):
        self._plot_3d_intensity(result)
        self._plot_2d_intensity(result)

    def _plot_3d_intensity(self, result):
        canvas = self.canvases["3d_intensity"]
        canvas.figure.clear()
        axis = canvas.figure.add_subplot(111, projection="3d")
        floor = result["parameters"].display_floor_db
        norm = Normalize(floor, 0.0)
        colors = cm.turbo(norm(np.clip(result["intensity_db"], floor, 0.0)))
        axis.plot_surface(
            result["x"], result["y"], result["z"], facecolors=colors,
            linewidth=0, antialiased=True, shade=False,
        )
        self._add_3d_polarization_axes(axis, result)
        self._draw_3d_incidence_geometry(axis, result)
        self._format_3d_axis(
            axis, "Far-field intensity with polarization directions"
        )
        colorbar = canvas.figure.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cm.turbo), ax=axis, shrink=0.70
        )
        colorbar.set_label("Relative DSC (dB)")
        canvas.draw_idle()

    def _add_3d_polarization_axes(self, axis, result):
        """Overlay Stokes polarization directions on the intensity surface."""
        index = self._sparse_indices(result)
        length = 0.12
        segments = []
        for x, y, z, ex, ey, ez, dolp in zip(
            result["x"][index].ravel(), result["y"][index].ravel(),
            result["z"][index].ravel(), result["e_x"][index].ravel(),
            result["e_y"][index].ravel(), result["e_z"][index].ravel(),
            result["dolp"][index].ravel(),
        ):
            half = 0.5 * length * max(0.25, dolp)
            segments.append(
                [(x - half * ex, y - half * ey, z - half * ez),
                 (x + half * ex, y + half * ey, z + half * ez)]
            )
        axis.add_collection3d(
            Line3DCollection(segments, colors="#1236b5", linewidths=1.4)
        )

    def _plot_2d_intensity(self, result):
        canvas = self.canvases["2d_intensity"]
        canvas.figure.clear()
        axis = canvas.figure.add_subplot(111)
        relative = result["intensity"] / np.max(result["intensity"])
        display_intensity = np.clip(relative, 0.0, 1.0) ** 0.80
        image = axis.pcolormesh(
            result["x"], result["y"], display_intensity,
            cmap=PARTICLE_INTENSITY_CMAP, vmin=0.0, vmax=1.0,
            shading="gouraud",
        )
        self._add_2d_polarization_axes(axis, result)
        self._format_2d_axis(
            axis, "Particle Scattering\nIntensity with Polarization Directions"
        )
        colorbar = canvas.figure.colorbar(image, ax=axis)
        colorbar.set_label("Normalized particle scattering intensity")
        canvas.draw_idle()

    def _add_2d_polarization_axes(self, axis, result):
        """Overlay projected Stokes polarization axes on a pupil plot."""
        index = self._sparse_indices(result)
        px = result["x"][index].ravel()
        py = result["y"][index].ravel()
        ex = result["e_x"][index].ravel()
        ey = result["e_y"][index].ravel()
        dolp = result["dolp"][index].ravel()
        projected_norm = np.hypot(ex, ey)
        valid = projected_norm > 1e-8
        segments = []
        for x, y, vx, vy, degree, magnitude in zip(
            px[valid], py[valid], ex[valid], ey[valid], dolp[valid],
            projected_norm[valid],
        ):
            half = 0.06 * max(0.25, degree)
            vx, vy = vx / magnitude, vy / magnitude
            segments.append(
                [(x - half * vx, y - half * vy),
                 (x + half * vx, y + half * vy)]
            )
        polarization_lines = LineCollection(
            segments, colors="#3157d5", linewidths=1.45
        )
        polarization_lines.set_clip_path(
            Circle((0.0, 0.0), 0.995, transform=axis.transData)
        )
        axis.add_collection(polarization_lines)

    @staticmethod
    def _draw_3d_incidence_geometry(axis, result):
        """Draw wafer, surface normal, laser ray, and incidence-angle arc."""
        incidence_deg = result["parameters"].incidence_deg
        incidence = np.deg2rad(incidence_deg)
        # All incidence geometry is anchored at the physical wafer origin.
        diagram_y = 0.0

        # A translucent wafer plane extending beyond the unit pupil.
        plane_x, plane_y = np.meshgrid(
            np.linspace(-1.45, 1.10, 2),
            np.linspace(-1.10, 1.10, 2),
        )
        axis.plot_surface(
            plane_x,
            plane_y,
            np.zeros_like(plane_x),
            color="#8d9499",
            alpha=0.16,
            shade=False,
            zorder=0,
        )
        axis.text(0.90, -0.95, 0.015, "wafer plane", color="#50565a")

        # Outward surface normal.
        axis.quiver(
            0.0, diagram_y, 0.0,
            0.0, 0.0, 0.72,
            color="#666666",
            linewidth=2.0,
            arrow_length_ratio=0.10,
            zorder=20,
        )
        axis.text(
            0.03, diagram_y, 0.76, "surface normal",
            color="#4f5559", zorder=20,
        )

        # Source lies above and to the left. The arrow points toward the wafer.
        source_distance = 1.40
        source_x = -source_distance * np.sin(incidence)
        source_z = source_distance * np.cos(incidence)
        axis.quiver(
            source_x, diagram_y, source_z,
            -source_x, 0.0, -source_z,
            color="#133e9a",
            linewidth=3.2,
            arrow_length_ratio=0.10,
            zorder=20,
        )
        axis.text(
            source_x - 0.02,
            diagram_y,
            source_z + 0.07,
            "incident laser",
            color="#133e9a",
            weight="bold",
            ha="center",
            zorder=20,
        )
        axis.scatter(
            [0.0], [0.0], [0.0],
            color="#133e9a", s=28, depthshade=False, zorder=21,
        )
        axis.text(
            0.04, 0.02, 0.03, "incidence point (0, 0, 0)",
            color="#133e9a", fontsize=8, zorder=21,
        )

        # Incidence angle is measured between the source ray and +z normal.
        arc_angle = np.linspace(0.0, incidence, 80)
        arc_radius = 0.30
        arc_x = -arc_radius * np.sin(arc_angle)
        arc_z = arc_radius * np.cos(arc_angle)
        axis.plot(
            arc_x,
            np.full_like(arc_x, diagram_y),
            arc_z,
            color="#133e9a",
            linewidth=1.8,
            zorder=20,
        )
        label_angle = incidence * 0.52
        axis.text(
            -0.38 * np.sin(label_angle),
            diagram_y,
            0.38 * np.cos(label_angle),
            "{:.1f} deg".format(incidence_deg),
            color="#133e9a",
            ha="center",
            va="center",
            zorder=20,
        )

    @staticmethod
    def _format_3d_axis(axis, title):
        axis.set(
            xlim=(-1.50, 1.10), ylim=(-1.25, 1.10), zlim=(0.0, 1.12),
            xlabel=r"$k_x/k$", ylabel=r"$k_y/k$", zlabel=r"$k_z/k$",
        )
        axis.set_title(title, fontsize=11, pad=4)
        axis.tick_params(labelsize=8, pad=1)
        axis.xaxis.label.set_size(9)
        axis.yaxis.label.set_size(9)
        axis.zaxis.label.set_size(9)
        axis.set_box_aspect((1, 1, 0.7))
        axis.view_init(elev=28, azim=-55)

    @staticmethod
    def _format_2d_axis(axis, title):
        horizon = np.linspace(0.0, 2.0 * np.pi, 361)
        axis.plot(np.cos(horizon), np.sin(horizon), color="white", linewidth=1.0)
        axis.text(
            -0.97, 0.94, "Particle", color="white", fontsize=13,
            weight="bold", ha="left", va="top", zorder=20,
        )
        axis.annotate(
            "",
            xy=(0.90, -1.20),
            xytext=(-0.90, -1.20),
            arrowprops=dict(arrowstyle="->", color="#123b91", lw=2.8),
            annotation_clip=False,
        )
        axis.text(
            0.0, -1.13, "Laser Incidence Direction",
            color="#123b91", fontsize=10, weight="bold",
            ha="center", va="bottom",
        )
        axis.set(
            xlim=(-1.05, 1.05), ylim=(-1.30, 1.05), aspect="equal",
            xlabel="", ylabel=r"$k_y/k$",
        )
        axis.set_title(title, fontsize=11, pad=4)
        axis.tick_params(labelsize=8, pad=1)
        axis.yaxis.label.set_size(9)
        axis.set_facecolor("#111111")

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
            key: value for key, value in self.result.items()
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


class ScatmechStandaloneWindow(QMainWindow):
    """Standalone wrapper retained for direct execution of this module."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("pySCATMECH Particle-on-Wafer Far-Field Simulator")
        self.resize(1500, 920)
        self.simulation_widget = ScatmechSimulationWidget()
        self.setCentralWidget(self.simulation_widget)

    def closeEvent(self, event):
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
    if "--headless-test" in sys.argv:
        index = sys.argv.index("--headless-test")
        output = sys.argv[index + 1] if index + 1 < len(sys.argv) else "ui_test_results"
        run_headless_test(output)
        return
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    window = ScatmechStandaloneWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
