# pages/page_wafer.py
from PyQt5.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout, QPushButton,
    QComboBox, QLabel, QCheckBox, QLineEdit
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from utils.ui_helpers import create_input_widget


class WaferPatternPage(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(15)

        col1_layout = QVBoxLayout()
        group_defect = QGroupBox("Defect Definition")
        form_defect = QFormLayout(group_defect)
        form_defect.setVerticalSpacing(15)

        btn_import = QPushButton("Import Defect ")
        form_defect.addRow(btn_import)

        combo_defect_type = QComboBox()
        combo_defect_type.addItems(["Particle", "PSL Sphere", "Scratch", "COP (Pit)", "Bridge / Short"])
        form_defect.addRow("Defect Type", combo_defect_type)

        form_defect.addRow("Defect Material (n+ik)", create_input_widget("1.5+0.0j", ""))
        form_defect.addRow("Height/Depth (Z)", create_input_widget("100", "nm"))
        form_defect.addRow("Dimensions (X, Y)", create_input_widget("100, 100", "nm"))
        form_defect.addRow("Coating Material", create_input_widget("SiO2", ""))
        form_defect.addRow("Coating Thickness", create_input_widget("10", "nm"))

        psd_parameters = QWidget()
        psd_layout = QHBoxLayout(psd_parameters)
        psd_layout.setContentsMargins(0, 0, 0, 0)
        psd_layout.setSpacing(6)

        self.psd_checkbox = QCheckBox("PSD")

        self.psd_inputs = {}
        for parameter in ("A", "B", "C"):
            psd_layout.addWidget(QLabel(parameter))
            parameter_input = QLineEdit("0.0")
            parameter_input.setPlaceholderText(parameter)
            parameter_input.setValidator(QDoubleValidator(-1e100, 1e100, 12, parameter_input))
            parameter_input.setMinimumWidth(60)
            parameter_input.setMaximumWidth(120)
            psd_layout.addWidget(parameter_input, stretch=1)
            self.psd_inputs[parameter] = parameter_input

        form_defect.addRow(self.psd_checkbox, psd_parameters)

        form_defect.addRow(QLabel(""))
        btn_generate = QPushButton("Generate 3D Model")
        btn_generate.setStyleSheet("background-color: #67C23A; color: white;")
        form_defect.addRow(btn_generate)

        col1_layout.addWidget(group_defect)
        col1_layout.addStretch()

        col2_layout = QVBoxLayout()
        group_3d = QGroupBox("3D Defect Topology Viewport")
        vbox_3d = QVBoxLayout(group_3d)

        render_viewport = QLabel(
            "3D Wafer & Defect Topology Render Viewport\n(Configure parameters on the left and click Generate)")
        render_viewport.setAlignment(Qt.AlignCenter)
        render_viewport.setStyleSheet("background-color: #0d1117; color: #409EFF; font-size: 13px; border-radius: 4px;")
        vbox_3d.addWidget(render_viewport)
        col2_layout.addWidget(group_3d)

        layout.addLayout(col1_layout, stretch=1)
        layout.addLayout(col2_layout, stretch=2)
