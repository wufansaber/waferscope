# pages/page_advanced.py
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout,
                             QPushButton, QComboBox, QLabel, QScrollArea, QFrame,
                             QTableWidget, QTableWidgetItem, QHeaderView, QSpinBox)
from PyQt5.QtCore import Qt
from utils.ui_helpers import create_input_widget


class AdvancedSettingsPage(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(15)

        # --- Left Col: Noise & Perturbations ---
        col1_layout = QVBoxLayout()
        scroll_col1 = QScrollArea()
        scroll_col1.setWidgetResizable(True)
        scroll_col1.setFrameShape(QFrame.NoFrame)
        col1_widget = QWidget()
        vbox_col1 = QVBoxLayout(col1_widget)
        vbox_col1.setSpacing(15)

        group_illum = QGroupBox("Illumination Perturbations")
        form_illum = QFormLayout(group_illum)
        form_illum.addRow("Power Stability", create_input_widget("0.86", "% (RMS)"))
        form_illum.addRow("Pointing Fluctuation", create_input_widget("20", "μm/1m"))
        form_illum.addRow("Spot Size Fluctuation", create_input_widget("250", "nm"))
        form_illum.addRow("Spot Center Drift", create_input_widget("250", "nm"))
        vbox_col1.addWidget(group_illum)

        group_cam = QGroupBox("Camera Noise")
        form_cam = QFormLayout(group_cam)
        form_cam.addRow("Readout Noise (RN)", create_input_widget("12.1", "e- @12bit"))
        form_cam.addRow("Dark Current", create_input_widget("400", "e-/s"))
        form_cam.addRow("PRNU", create_input_widget("4.43", "e-"))
        form_cam.addRow("DSNU", create_input_widget("0.03", "%"))
        vbox_col1.addWidget(group_cam)

        group_vib = QGroupBox("Mechanical Vibration Noise")
        form_vib = QFormLayout(group_vib)
        form_vib.addRow("Scan Vel. Uniformity", create_input_widget("0.02", "%"))
        form_vib.addRow("Low-Freq Res. (X)", create_input_widget("500", "nm"))
        form_vib.addRow("Low-Freq Res. (Y)", create_input_widget("500", "nm"))
        form_vib.addRow("Low-Freq Res. (Z)", create_input_widget("1000", "nm"))
        form_vib.addRow("Dynamic HF Jitter (X)", create_input_widget("500", "nm"))
        form_vib.addRow("Dynamic HF Jitter (Y)", create_input_widget("20", "nm"))
        form_vib.addRow("Dynamic HF Jitter (Z)", create_input_widget("20", "nm"))
        form_vib.addRow("Defocus Jitter Range", create_input_widget("< 20", "nm"))
        vbox_col1.addWidget(group_vib)
        vbox_col1.addStretch()

        scroll_col1.setWidget(col1_widget)
        col1_layout.addWidget(scroll_col1)

        # --- Right Col: Fourier Filter & Zernike ---
        col2_layout = QVBoxLayout()
        col2_layout.setSpacing(15)

        group_fourier = QGroupBox("Fourier Pupil Filter Configuration (Max 11 Bars)")
        hbox_fourier = QHBoxLayout(group_fourier)

        scroll_fourier = QScrollArea()
        scroll_fourier.setWidgetResizable(True)
        scroll_fourier.setFrameShape(QFrame.NoFrame)
        scroll_fourier.setFixedWidth(340)

        fourier_content = QWidget()
        form_fourier = QFormLayout(fourier_content)

        spin_bars = QSpinBox()
        spin_bars.setRange(0, 11)
        spin_bars.setValue(3)
        form_fourier.addRow("Total Bars Amount", spin_bars)

        # Bar 0
        lbl_bar0 = QLabel("▶ Bar 0 (Cylindrical)")
        lbl_bar0.setStyleSheet("font-weight: bold; color: #409EFF; padding-top: 10px;")
        form_fourier.addRow(lbl_bar0)
        form_fourier.addRow("X Position", create_input_widget("0.0", "mm"))
        form_fourier.addRow("Diameter", create_input_widget("2.0", "mm"))
        form_fourier.addRow("Length", create_input_widget("10.0", "mm"))

        # Bars 1-10
        lbl_bar1_10 = QLabel("▶ Bars 1-10 (Triangular Toothed)")
        lbl_bar1_10.setStyleSheet("font-weight: bold; color: #409EFF; padding-top: 10px;")
        form_fourier.addRow(lbl_bar1_10)

        combo_bar_id = QComboBox()
        combo_bar_id.addItems([f"Bar {i}" for i in range(1, 11)])
        form_fourier.addRow("Select Bar ID", combo_bar_id)
        form_fourier.addRow("X Position", create_input_widget("1.5", "mm"))
        form_fourier.addRow("Period", create_input_widget("0.5", "mm"))
        form_fourier.addRow("Coeff A", create_input_widget("1.0", ""))
        form_fourier.addRow("Coeff B", create_input_widget("0.5", ""))
        form_fourier.addRow("Coeff C", create_input_widget("0.2", ""))

        btn_apply_fourier = QPushButton("Apply Parameters")
        btn_apply_fourier.setStyleSheet("background-color: #E6A23C; color: white; margin-top: 10px;")
        form_fourier.addRow(btn_apply_fourier)

        scroll_fourier.setWidget(fourier_content)
        hbox_fourier.addWidget(scroll_fourier)

        vbox_fourier_view = QVBoxLayout()
        combo_fourier_channel = QComboBox()
        combo_fourier_channel.addItems([
            "View: Channel 1 (Edge -45°)",
            "View: Channel 2 (Center 0°)",
            "View: Channel 3 (Edge +45°)"
        ])
        combo_fourier_channel.setFixedHeight(32)
        combo_fourier_channel.setStyleSheet("font-weight: bold;")

        lbl_fourier_view = QLabel(
            "Pupil Plane Field Intensity & Mask\n[Currently displaying Channel 1]\n\n(Waiting for Diffraction Orders...)")
        lbl_fourier_view.setAlignment(Qt.AlignCenter)
        lbl_fourier_view.setStyleSheet(
            "background-color: #0d1117; color: #409EFF; font-size: 16px; border-radius: 4px;")

        vbox_fourier_view.addWidget(combo_fourier_channel)
        vbox_fourier_view.addWidget(lbl_fourier_view, stretch=1)
        hbox_fourier.addLayout(vbox_fourier_view, stretch=1)

        # Zernike
        group_zernike = QGroupBox("Collection System Wavefront Aberration")
        hbox_zernike = QHBoxLayout(group_zernike)
        vbox_zernike_params = QVBoxLayout()
        table_zernike = QTableWidget(36, 2)
        table_zernike.setHorizontalHeaderLabels(["Order", "Coefficient"])
        table_zernike.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table_zernike.verticalHeader().setVisible(False)
        for i in range(36):
            table_zernike.setItem(i, 0, QTableWidgetItem(f"Z{i + 1}"))
            val = "0.000" if i > 8 else str(round((9 - i) * 0.01, 3))
            table_zernike.setItem(i, 1, QTableWidgetItem(val))

        btn_apply_zernike = QPushButton("Fit Wavefront Aberration")
        btn_apply_zernike.setStyleSheet("background-color: #67C23A; color: white;")
        vbox_zernike_params.addWidget(table_zernike)
        vbox_zernike_params.addWidget(btn_apply_zernike)

        lbl_zernike_view = QLabel("2D Wavefront Aberration Map\nFitted by Zernike Polynomials")
        lbl_zernike_view.setAlignment(Qt.AlignCenter)
        lbl_zernike_view.setStyleSheet(
            "background-color: #0d1117; color: #E6A23C; font-size: 16px; border-radius: 4px;")

        hbox_zernike.addLayout(vbox_zernike_params, stretch=1)
        hbox_zernike.addWidget(lbl_zernike_view, stretch=2)

        col2_layout.addWidget(group_fourier, stretch=1)
        col2_layout.addWidget(group_zernike, stretch=1)

        layout.addLayout(col1_layout, stretch=1)
        layout.addLayout(col2_layout, stretch=2)