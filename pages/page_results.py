# pages/page_results.py
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel, QPushButton, QComboBox
from PyQt5.QtCore import Qt


class ResultsPage(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(15)

        # --- Top Half: Near-Field & Far-Field (1:1) ---
        top_layout = QHBoxLayout()
        top_layout.setSpacing(15)

        group_near = QGroupBox("Wafer Surface Near-Field 3D Results")
        vbox_near = QVBoxLayout(group_near)
        lbl_near = QLabel(
            "Localized 3D Electromagnetic Scattering Enhancement Field\n[Including Near-Field Distortion of Micro-Defects]\n\n(Waiting for calculation...)")
        lbl_near.setAlignment(Qt.AlignCenter)
        lbl_near.setStyleSheet("background-color: #0d1117; color: #67C23A; font-size: 16px; border-radius: 4px;")
        vbox_near.addWidget(lbl_near)
        top_layout.addWidget(group_near, stretch=1)

        group_pupil = QGroupBox("Far-Field Pupil Diffraction Intensity Results")
        vbox_pupil = QVBoxLayout(group_pupil)

        combo_res_pupil_ch = QComboBox()
        combo_res_pupil_ch.addItems([
            "Display Result: Channel 1 (Edge -45°)",
            "Display Result: Channel 2 (Center 0°)",
            "Display Result: Channel 3 (Edge +45°)"
        ])
        combo_res_pupil_ch.setFixedHeight(32)
        combo_res_pupil_ch.setStyleSheet("font-weight: bold;")

        lbl_pupil = QLabel(
            "Objective Back Focal Plane Diffraction Orders & Background\n[Channel 1 Filtered Results]\n\n(Waiting for calculation...)")
        lbl_pupil.setAlignment(Qt.AlignCenter)
        lbl_pupil.setStyleSheet("background-color: #0d1117; color: #E6A23C; font-size: 16px; border-radius: 4px;")

        vbox_pupil.addWidget(combo_res_pupil_ch)
        vbox_pupil.addWidget(lbl_pupil, stretch=1)
        top_layout.addWidget(group_pupil, stretch=1)

        layout.addLayout(top_layout, stretch=3)

        # --- Middle: TDI Imaging Flow ---
        group_imaging = QGroupBox("12-bit TDI Dynamic Overexposure Scanning Image Flow")
        vbox_imaging = QVBoxLayout(group_imaging)

        combo_res_tdi_ch = QComboBox()
        combo_res_tdi_ch.addItems([
            "Imaging Output: Combined (Fused)",
            "Imaging Output: Channel 1",
            "Imaging Output: Channel 2",
            "Imaging Output: Channel 3"
        ])
        combo_res_tdi_ch.setFixedHeight(32)
        combo_res_tdi_ch.setStyleSheet("font-weight: bold;")

        lbl_imaging = QLabel(
            "12-bit TDI Grayscale Rolling Scan Pattern\n(Supports 4095 Full-Well Truncation Red Alert / Includes Full-Link Noise)\n\n(Waiting for engine execution...)")
        lbl_imaging.setAlignment(Qt.AlignCenter)
        lbl_imaging.setStyleSheet(
            "background-color: #1a1a1a; color: #FFFFFF; font-size: 18px; border: 1px solid #444; border-radius: 4px;")

        vbox_imaging.addWidget(combo_res_tdi_ch)
        vbox_imaging.addWidget(lbl_imaging, stretch=1)

        layout.addWidget(group_imaging, stretch=4)

        # --- Bottom: SNR Metrics & Actions ---
        group_snr = QGroupBox("System-Level Defect Comprehensive Evaluation Metrics")
        group_snr.setStyleSheet("QGroupBox { background-color: #eef1f6; border: 1px solid #d1dbe5; }")
        hbox_snr = QHBoxLayout(group_snr)

        vbox_text = QVBoxLayout()
        lbl_snr_title = QLabel("Final Defect SNR:")
        lbl_snr_title.setStyleSheet("font-size: 18px; color: #333333;")
        lbl_snr_value = QLabel("-- dB")
        lbl_snr_value.setStyleSheet("font-size: 32px; font-weight: bold; color: #F56C6C;")
        lbl_haze_info = QLabel("Background Mean: -- DN  |  Background StdDev: -- DN  |  Defect Peak: -- DN")
        lbl_haze_info.setStyleSheet("font-size: 14px; color: #888888;")

        vbox_text.addWidget(lbl_snr_title)
        vbox_text.addWidget(lbl_snr_value)
        vbox_text.addWidget(lbl_haze_info)

        vbox_btns = QVBoxLayout()
        btn_run = QPushButton("▶ Run Full-Link Simulation")
        btn_run.setStyleSheet(
            "background-color: #409EFF; color: white; font-size: 16px; padding: 15px; border-radius: 6px;")
        btn_export = QPushButton("Export Simulation Results & Images")
        btn_export.setStyleSheet(
            "background-color: #909399; color: white; font-size: 14px; padding: 10px; border-radius: 4px;")

        vbox_btns.addWidget(btn_run)
        vbox_btns.addWidget(btn_export)

        hbox_snr.addLayout(vbox_text, stretch=2)
        hbox_snr.addStretch(1)
        hbox_snr.addLayout(vbox_btns, stretch=1)

        layout.addWidget(group_snr, stretch=2)