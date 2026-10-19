"""Qt controls for multilayer coatings and SCATMECH microroughness."""

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .materials import COATING_MATERIALS
from .parameters import CoatingLayerParameters


ROUGHNESS_MODE_LABELS = (
    ("Top exposed interface", "top_interface"),
    ("All interfaces correlated", "correlated_interfaces"),
    ("All interfaces uncorrelated", "uncorrelated_interfaces"),
)


class CoatingStackGroup(QGroupBox):
    """Editable coating list ordered from the wafer toward the exposed surface."""

    def __init__(self):
        super().__init__("Wafer coating stack")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 10, 8, 6)
        layout.setSpacing(5)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ("#", "Material", "n", "k", "Thickness (nm)")
        )
        self.table.setMinimumHeight(132)
        self.table.setMaximumHeight(210)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setToolTip(
            "Rows are ordered from the wafer substrate toward the exposed surface."
        )
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Fixed)
        self.table.setColumnWidth(0, 26)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.setColumnWidth(2, 58)
        self.table.setColumnWidth(3, 58)
        self.table.setColumnWidth(4, 96)
        layout.addWidget(self.table)

        buttons = QHBoxLayout()
        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")
        self.up_button = QPushButton("Up")
        self.down_button = QPushButton("Down")
        for button in (
            self.add_button,
            self.remove_button,
            self.up_button,
            self.down_button,
        ):
            buttons.addWidget(button)
        self.add_button.clicked.connect(self.add_layer)
        self.remove_button.clicked.connect(self.remove_selected_layer)
        self.up_button.clicked.connect(lambda: self.move_selected_layer(-1))
        self.down_button.clicked.connect(lambda: self.move_selected_layer(1))
        layout.addLayout(buttons)

    def add_layer(self, layer=None):
        if self.table.rowCount() >= 32:
            return
        layer = layer or CoatingLayerParameters()
        row = self.table.rowCount()
        self.table.insertRow(row)
        self._set_row(row, layer)
        self.table.selectRow(row)

    def remove_selected_layer(self):
        row = self.table.currentRow()
        if row >= 0:
            self.table.removeRow(row)
            self._renumber_rows()

    def move_selected_layer(self, offset):
        row = self.table.currentRow()
        target = row + offset
        if row < 0 or target < 0 or target >= self.table.rowCount():
            return
        layers = list(self.layers())
        layers[row], layers[target] = layers[target], layers[row]
        self.set_layers(layers)
        self.table.selectRow(target)

    def layers(self):
        result = []
        for row in range(self.table.rowCount()):
            material = self.table.cellWidget(row, 1)
            n_edit = self.table.cellWidget(row, 2)
            k_edit = self.table.cellWidget(row, 3)
            thickness = self.table.cellWidget(row, 4)
            result.append(
                CoatingLayerParameters(
                    material=material.currentText(),
                    n=float(n_edit.text()),
                    k=float(k_edit.text()),
                    thickness_nm=thickness.value(),
                )
            )
        return tuple(result)

    def set_layers(self, layers):
        self.table.setRowCount(0)
        for layer in layers:
            self.add_layer(layer)

    def _set_row(self, row, layer):
        number = QTableWidgetItem(str(row + 1))
        number.setToolTip("1 is nearest the wafer; the last row is the top layer.")
        self.table.setItem(row, 0, number)

        material = QComboBox()
        material.addItems(COATING_MATERIALS.keys())
        material.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        material.setMinimumContentsLength(7)
        index = material.findText(layer.material)
        material.setCurrentIndex(index if index >= 0 else material.findText("Custom"))

        n_edit = QLineEdit("{:.8g}".format(layer.n))
        k_edit = QLineEdit("{:.8g}".format(layer.k))
        n_edit.setMaximumWidth(58)
        k_edit.setMaximumWidth(58)
        validator = QDoubleValidator(0.0, 100.0, 8, self)
        n_edit.setValidator(validator)
        k_edit.setValidator(validator)

        thickness = QDoubleSpinBox()
        thickness.setRange(0.001, 1_000_000.0)
        thickness.setDecimals(6)
        thickness.setValue(layer.thickness_nm)
        thickness.setMaximumWidth(96)

        def apply_material(name):
            value = COATING_MATERIALS[name]
            n_edit.setText("{:.8g}".format(value.real))
            k_edit.setText("{:.8g}".format(value.imag))

        material.currentTextChanged.connect(apply_material)
        self.table.setCellWidget(row, 1, material)
        self.table.setCellWidget(row, 2, n_edit)
        self.table.setCellWidget(row, 3, k_edit)
        self.table.setCellWidget(row, 4, thickness)

    def _renumber_rows(self):
        for row in range(self.table.rowCount()):
            self.table.setItem(row, 0, QTableWidgetItem(str(row + 1)))


class RoughnessGroup(QGroupBox):
    """Controls for SCATMECH Rayleigh-Rice roughness and beam footprint."""

    def __init__(self):
        super().__init__("Surface and interface roughness")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 10, 8, 6)
        layout.setSpacing(5)

        self.enabled = QCheckBox("Include roughness background")
        layout.addWidget(self.enabled)
        self.match_fem_button = QPushButton("Match displayed FEM haze result")
        self.match_fem_button.setToolTip(
            "Copy parameters from the displayed FEM result. Gaussian FEM roughness "
            "is required; a legacy surface has no equivalent Gaussian L.")
        self.match_fem_button.hide()
        layout.addWidget(self.match_fem_button)

        form = QWidget()
        form_layout = QFormLayout(form)
        form_layout.setRowWrapPolicy(QFormLayout.WrapLongRows)
        form_layout.setContentsMargins(0, 0, 0, 0)
        form_layout.setHorizontalSpacing(6)
        form_layout.setVerticalSpacing(4)
        self.mode = _labeled_combo(form_layout, "Interfaces")
        for label, value in ROUGHNESS_MODE_LABELS:
            self.mode.addItem(label, value)
        self.haze_only = QCheckBox("Haze only (no particle)")
        form_layout.addRow(self.haze_only)
        self.psd = _labeled_combo(form_layout, "PSD model")
        self.psd.addItem("ABC", "abc")
        self.psd.addItem("Gaussian", "gaussian")
        self.psd.setCurrentIndex(self.psd.findData("gaussian"))
        self.rms_nm = _labeled_spinbox(
            form_layout,
            "RMS height",
            0.001,
            100.0,
            0.2,
            " nm",
            4,
        )
        self.correlation_length_um = _labeled_spinbox(
            form_layout,
            "PSD length B / L",
            0.001,
            10000.0,
            0.2,
            " um",
            4,
        )
        self.exponent = _labeled_spinbox(
            form_layout,
            "ABC exponent C",
            2.01,
            10.0,
            3.0,
            "",
            3,
        )
        self.footprint_length_um = _labeled_spinbox(
            form_layout,
            "Scattering area length",
            0.001,
            1_000_000.0,
            10.0,
            " um",
            3,
        )
        self.footprint_width_um = _labeled_spinbox(
            form_layout,
            "Scattering area width",
            0.001,
            1_000_000.0,
            10.0,
            " um",
            3,
        )
        self.area_label = QLabel()
        self.area_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        form_layout.addRow("Effective area", self.area_label)
        layout.addWidget(form)

        self._parameter_form = form
        self.enabled.toggled.connect(form.setEnabled)
        self.enabled.toggled.connect(self._roughness_toggled)
        self.psd.currentIndexChanged.connect(self._update_psd_controls)
        self.footprint_length_um.valueChanged.connect(self._update_area)
        self.footprint_width_um.valueChanged.connect(self._update_area)
        form.setEnabled(False)
        self._update_psd_controls()
        self._update_area()

        self.mode.setToolTip(
            "Top uses the exposed interface only. Correlated and uncorrelated "
            "modes include every interface in the SCATMECH film stack."
        )
        self.rms_nm.setToolTip("Stationary RMS height before FEM patch edge taper.")
        self.correlation_length_um.setToolTip(
            "Gaussian: correlation length L in C(r)=RMS^2 exp[-(r/L)^2]. "
            "ABC: PSD roll-off parameter B; B and L are not equivalent."
        )
        self.exponent.setToolTip("ABC PSD exponent C must be greater than 2.")
        self.footprint_length_um.setToolTip(
            "Modeled wafer area length. FEM haze matching uses its finite rough patch, "
            "not the larger laser power-normalization spot."
        )
        self.footprint_width_um.setToolTip(
            "Modeled wafer area width. FEM haze matching uses its finite rough patch, "
            "not the larger laser power-normalization spot."
        )

    def values(self):
        return {
            "roughness_enabled": self.enabled.isChecked(),
            "haze_only": self.haze_only.isChecked(),
            "roughness_psd": self.psd.currentData(),
            "roughness_interface_mode": self.mode.currentData(),
            "roughness_rms_nm": self.rms_nm.value(),
            "roughness_correlation_length_um": self.correlation_length_um.value(),
            "roughness_exponent": self.exponent.value(),
            "illumination_length_um": self.footprint_length_um.value(),
            "illumination_width_um": self.footprint_width_um.value(),
        }

    def _update_area(self):
        area = self.footprint_length_um.value() * self.footprint_width_um.value()
        self.area_label.setText("{:.6g} um^2".format(area))

    def _roughness_toggled(self, enabled):
        if not enabled:
            self.haze_only.setChecked(False)

    def _update_psd_controls(self, *_):
        self.exponent.setEnabled(self.psd.currentData() == "abc")


def _labeled_combo(layout, label):
    combo = QComboBox()
    layout.addRow(label, combo)
    return combo


def _labeled_spinbox(layout, label, minimum, maximum, value, suffix, decimals):
    spinbox = QDoubleSpinBox()
    spinbox.setRange(minimum, maximum)
    spinbox.setDecimals(decimals)
    spinbox.setValue(value)
    spinbox.setSuffix(suffix)
    layout.addRow(label, spinbox)
    return spinbox
