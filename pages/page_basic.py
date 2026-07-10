# pages/page_basic.py
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QScrollArea, QFrame, QFormLayout, QPushButton, \
    QComboBox
from utils.ui_helpers import create_input_widget


class BasicSettingsPage(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(20)

        # --- Left Col: Lighting & Laser ---
        col1_layout = QVBoxLayout()
        group_lighting = QGroupBox("Lighting & Laser Configuration")
        scroll_lighting = QScrollArea()
        scroll_lighting.setWidgetResizable(True)
        scroll_lighting.setFrameShape(QFrame.NoFrame)

        lighting_content = QWidget()
        form_lighting = QFormLayout(lighting_content)
        form_lighting.setVerticalSpacing(12)

        form_lighting.addRow("Nominal Wavelength", create_input_widget("355", "nm"))
        form_lighting.addRow("Spectral Bandwidth", create_input_widget("50", "pm (FWHM)"))
        form_lighting.addRow("Total Laser Power", create_input_widget("5", "W"))
        form_lighting.addRow("Wafer Spot Size", create_input_widget("L: 12mm / W: 2μm", ""))

        btn_import_laser = QPushButton("Import Laser Parameter List")
        btn_import_laser.setStyleSheet("background-color: #67C23A; color: white;")
        form_lighting.addRow(btn_import_laser)

        btn_import_spot = QPushButton("Import Spot Intensity Matrix")
        btn_import_spot.setStyleSheet("background-color: #67C23A; color: white;")
        form_lighting.addRow(btn_import_spot)

        combo_path = QComboBox()
        combo_path.addItems(["Oblique (Edge Path)", "Normal (Center Path)"])
        form_lighting.addRow("Illumination Path", combo_path)
        form_lighting.addRow("Oblique Angle", create_input_widget("75", "°"))
        form_lighting.addRow("Illumination NA", create_input_widget("0.24", ""))

        combo_pol_mode = QComboBox()
        combo_pol_mode.addItems(["P-Polarized", "S-Polarized", "Circular (C)"])
        form_lighting.addRow("Polarization Mode", combo_pol_mode)
        form_lighting.addRow("Polarization Ratio", create_input_widget("100:1", ""))

        scroll_lighting.setWidget(lighting_content)
        vbox_lighting = QVBoxLayout(group_lighting)
        vbox_lighting.addWidget(scroll_lighting)
        col1_layout.addWidget(group_lighting)

        # --- Right Col: Collecting + Sensor ---
        col2_layout = QVBoxLayout()
        col2_layout.setSpacing(20)

        group_collecting = QGroupBox("Collecting Configuration")
        form_collecting = QFormLayout(group_collecting)
        form_collecting.setVerticalSpacing(10)
        combo_mag = QComboBox()
        combo_mag.addItems(["40X", "20X", "10X"])
        form_collecting.addRow("Magnification", combo_mag)
        combo_fov = QComboBox()
        combo_fov.addItems(["5.2 mm", "3.8 mm", "2.5 mm"])
        form_collecting.addRow("Object FOV", combo_fov)
        form_collecting.addRow("Numerical Aperture", create_input_widget("Edge 0.65 / Ctr 0.5", "(Comb 0.997)"))
        form_collecting.addRow("Objective Eq. Focal Length", create_input_widget("200", "mm"))
        form_collecting.addRow("Collection Tilt", create_input_widget("45°", ""))

        combo_col_pol = QComboBox()
        combo_col_pol.addItems(["P", "S", "N (Unpolarized)", "C (Circular)"])
        form_collecting.addRow("Polarization (Analyzer)", combo_col_pol)

        form_collecting.addRow("Ch1 Attenuation Filter", create_input_widget("100", "%"))
        form_collecting.addRow("Ch2 Attenuation Filter", create_input_widget("100", "%"))
        form_collecting.addRow("Ch3 Attenuation Filter", create_input_widget("100", "%"))

        btn_import_pupil = QPushButton("Import Actual Pupil Matrix")
        btn_import_pupil.setStyleSheet("background-color: #67C23A; color: white;")
        form_collecting.addRow(btn_import_pupil)
        col2_layout.addWidget(group_collecting)

        group_sensor = QGroupBox("TDI Sensor Configuration")
        form_sensor = QFormLayout(group_sensor)
        form_sensor.setVerticalSpacing(10)
        form_sensor.addRow("Pixel Physical Size", create_input_widget("5 x 5", "μm"))
        form_sensor.addRow("Array Resolution", create_input_widget("8028 x 256", ""))
        self.input_scan_image_size = create_input_widget("1024 x 1024", "pixels")
        form_sensor.addRow("Scan Image Size", self.input_scan_image_size)
        form_sensor.addRow("Line Rate", create_input_widget("1000", "kHz"))
        form_sensor.addRow("Integration Stages", create_input_widget("100", "Stages"))
        form_sensor.addRow("Quantum Efficiency", create_input_widget("55", "% @ 355nm"))
        form_sensor.addRow("Full Well Capacity", create_input_widget("16,000", "e-"))
        form_sensor.addRow("Conversion Gain", create_input_widget("0.6", "SDN/e-"))
        col2_layout.addWidget(group_sensor)

        layout.addLayout(col1_layout, stretch=1)
        layout.addLayout(col2_layout, stretch=1)
