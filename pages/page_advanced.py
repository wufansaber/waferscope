# pages/page_advanced.py
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout,
                             QPushButton, QComboBox, QLabel, QScrollArea, QFrame,
                             QTableWidget, QTableWidgetItem, QHeaderView, QSpinBox,
                             QMessageBox, QSizePolicy, QLineEdit)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle
from utils.ui_helpers import create_input_widget
from utils.zernike import evaluate_zernike_wavefront


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
        scroll_fourier.setFixedWidth(300)

        fourier_content = QWidget()
        form_fourier = QFormLayout(fourier_content)

        self.spin_fourier_bars = QSpinBox()
        self.spin_fourier_bars.setRange(0, 11)
        self.spin_fourier_bars.setValue(3)
        form_fourier.addRow("Total Bars Amount", self.spin_fourier_bars)

        # Bar 0
        lbl_bar0 = QLabel("▶ Bar 0 (Cylindrical)")
        lbl_bar0.setStyleSheet("font-weight: bold; color: #409EFF; padding-top: 10px;")
        form_fourier.addRow(lbl_bar0)
        self.input_bar0_x = create_input_widget("0.0", "mm")
        self.input_bar0_diameter = create_input_widget("2.0", "mm")
        self.input_bar0_length = create_input_widget("10.0", "mm")
        form_fourier.addRow("X Position", self.input_bar0_x)
        form_fourier.addRow("Diameter", self.input_bar0_diameter)
        form_fourier.addRow("Length", self.input_bar0_length)

        # Bars 1-10
        lbl_bar1_10 = QLabel("▶ Bars 1-10 (Triangular Toothed)")
        lbl_bar1_10.setStyleSheet("font-weight: bold; color: #409EFF; padding-top: 10px;")
        form_fourier.addRow(lbl_bar1_10)

        self.combo_fourier_bar_id = QComboBox()
        self.combo_fourier_bar_id.addItems([f"Bar {i}" for i in range(1, 11)])
        form_fourier.addRow("Select Bar ID", self.combo_fourier_bar_id)
        self.input_fourier_x = create_input_widget("1.5", "mm")
        self.input_fourier_period = create_input_widget("0.5", "mm")
        self.input_fourier_a = create_input_widget("1.0", "mm")
        self.input_fourier_b = create_input_widget("0.5", "mm")
        self.input_fourier_c = create_input_widget("0.2", "mm")
        form_fourier.addRow("X Position", self.input_fourier_x)
        form_fourier.addRow("Period", self.input_fourier_period)
        form_fourier.addRow("Coeff A (tooth height)", self.input_fourier_a)
        form_fourier.addRow("Coeff B (tooth width)", self.input_fourier_b)
        form_fourier.addRow("Coeff C (bar thickness)", self.input_fourier_c)

        self.btn_apply_fourier = QPushButton("Apply Parameters")
        self.btn_apply_fourier.setStyleSheet("background-color: #E6A23C; color: white; margin-top: 10px;")
        self.btn_apply_fourier.clicked.connect(self._apply_fourier_parameters)
        form_fourier.addRow(self.btn_apply_fourier)

        scroll_fourier.setWidget(fourier_content)
        hbox_fourier.addWidget(scroll_fourier)

        vbox_fourier_view = QVBoxLayout()
        self.combo_fourier_channel = QComboBox()
        self.combo_fourier_channel.addItems([
            "View: Channel 1 (Edge -45°)",
            "View: Channel 2 (Center 0°)",
            "View: Channel 3 (Edge +45°)"
        ])
        self.combo_fourier_channel.setFixedHeight(32)
        self.combo_fourier_channel.setStyleSheet("font-weight: bold;")
        self.combo_fourier_channel.currentIndexChanged.connect(
            self._refresh_fourier_channel
        )

        self.fourier_figure = Figure(figsize=(7, 3), constrained_layout=True)
        self.fourier_canvas = FigureCanvasQTAgg(self.fourier_figure)
        self.fourier_canvas.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding
        )
        self._fourier_parameters = {
            "x_position": 1.5,
            "period": 0.5,
            "a": 1.0,
            "b": 0.5,
            "c": 0.2,
        }
        self._render_fourier_bar(self._fourier_parameters)

        vbox_fourier_view.addWidget(self.combo_fourier_channel)
        vbox_fourier_view.addWidget(self.fourier_canvas, stretch=1)
        hbox_fourier.addLayout(vbox_fourier_view, stretch=1)

        # Zernike
        group_zernike = QGroupBox("Collection System Wavefront Aberration")
        hbox_zernike = QHBoxLayout(group_zernike)
        vbox_zernike_params = QVBoxLayout()
        self.table_zernike = QTableWidget(36, 2)
        self.table_zernike.setHorizontalHeaderLabels(["Order", "Coefficient"])
        self.table_zernike.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_zernike.verticalHeader().setVisible(False)
        for i in range(36):
            self.table_zernike.setItem(i, 0, QTableWidgetItem(f"Z{i + 1}"))
            val = "0.000" if i > 8 else str(round((9 - i) * 0.01, 3))
            self.table_zernike.setItem(i, 1, QTableWidgetItem(val))

        self.btn_apply_zernike = QPushButton("Fit Wavefront Aberration")
        self.btn_apply_zernike.setStyleSheet("background-color: #67C23A; color: white;")
        self.btn_apply_zernike.clicked.connect(self._apply_zernike_coefficients)
        vbox_zernike_params.addWidget(self.table_zernike)
        vbox_zernike_params.addWidget(self.btn_apply_zernike)

        zernike_result_layout = QVBoxLayout()
        zernike_result_layout.setContentsMargins(0, 0, 0, 0)
        zernike_result_layout.setSpacing(3)
        self.lbl_zernike_metrics = QLabel(
            "Click Fit Wavefront Aberration to render the 3-D wavefront"
        )
        self.lbl_zernike_metrics.setStyleSheet(
            "color: #606266; font-size: 10px;"
        )
        self.zernike_figure = Figure(figsize=(7, 4), constrained_layout=True)
        self.zernike_canvas = FigureCanvasQTAgg(self.zernike_figure)
        self.zernike_canvas.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding
        )
        placeholder_axis = self.zernike_figure.add_subplot(111)
        placeholder_axis.set_facecolor("#0d1117")
        placeholder_axis.text(
            0.5, 0.5,
            "3-D Zernike Wavefront Aberration",
            color="#E6A23C", ha="center", va="center",
            transform=placeholder_axis.transAxes,
        )
        placeholder_axis.set_xticks([])
        placeholder_axis.set_yticks([])
        zernike_result_layout.addWidget(self.lbl_zernike_metrics)
        zernike_result_layout.addWidget(self.zernike_canvas, stretch=1)

        hbox_zernike.addLayout(vbox_zernike_params, stretch=1)
        hbox_zernike.addLayout(zernike_result_layout, stretch=2)

        col2_layout.addWidget(group_fourier, stretch=1)
        col2_layout.addWidget(group_zernike, stretch=1)

        layout.addLayout(col1_layout, stretch=1)
        layout.addLayout(col2_layout, stretch=2)

    @staticmethod
    def _input_value(input_widget, field_name):
        line_edit = input_widget.findChild(QLineEdit)
        text = line_edit.text().strip() if line_edit is not None else ""
        try:
            return float(text)
        except ValueError as error:
            raise ValueError("{} must be numeric.".format(field_name)) from error

    def _apply_fourier_parameters(self):
        try:
            parameters = {
                "x_position": self._input_value(
                    self.input_fourier_x, "X Position"
                ),
                "period": self._input_value(
                    self.input_fourier_period, "Period"
                ),
                "a": self._input_value(self.input_fourier_a, "Coeff A"),
                "b": self._input_value(self.input_fourier_b, "Coeff B"),
                "c": self._input_value(self.input_fourier_c, "Coeff C"),
            }
        except ValueError as error:
            QMessageBox.warning(self, "Invalid Fourier bar parameter", str(error))
            return

        if parameters["period"] <= 0 or parameters["a"] <= 0:
            QMessageBox.warning(
                self,
                "Invalid Fourier bar parameter",
                "Period and coefficient A must be greater than zero.",
            )
            return
        if parameters["b"] <= 0 or parameters["b"] > parameters["period"]:
            QMessageBox.warning(
                self,
                "Invalid Fourier bar parameter",
                "Coefficient B must be greater than zero and no larger than Period.",
            )
            return
        if parameters["c"] <= 0:
            QMessageBox.warning(
                self,
                "Invalid Fourier bar parameter",
                "Coefficient C must be greater than zero.",
            )
            return

        self._fourier_parameters = parameters
        self._render_fourier_bar(parameters)

    def _refresh_fourier_channel(self, _index=None):
        self._render_fourier_bar(self._fourier_parameters)

    def _render_fourier_bar(self, parameters):
        """Render the selected triangular-toothed Fourier stop geometry."""
        period = parameters["period"]
        tooth_height = parameters["a"]
        tooth_width = parameters["b"]
        bar_thickness = parameters["c"]
        tooth_count = 7
        total_length = tooth_count * period
        left = -total_length / 2.0
        right = total_length / 2.0
        blue = "#4f79c7"
        edge = "#2d5aa7"

        self.fourier_figure.clear()
        axis = self.fourier_figure.add_subplot(111)
        self.fourier_figure.patch.set_facecolor("#0d1117")
        axis.set_facecolor("#0d1117")
        axis.add_patch(
            Rectangle(
                (left, -bar_thickness / 2.0),
                total_length,
                bar_thickness,
                facecolor=blue,
                edgecolor=edge,
                linewidth=1.2,
                zorder=3,
            )
        )

        centers = [
            left + (index + 0.5) * period for index in range(tooth_count)
        ]
        for center in centers:
            axis.add_patch(
                Polygon(
                    [
                        (center - tooth_width / 2.0, bar_thickness / 2.0),
                        (center, bar_thickness / 2.0 + tooth_height),
                        (center + tooth_width / 2.0, bar_thickness / 2.0),
                    ],
                    closed=True,
                    facecolor=blue,
                    edgecolor=edge,
                    linewidth=1.0,
                    zorder=2,
                )
            )
            axis.add_patch(
                Polygon(
                    [
                        (center - tooth_width / 2.0, -bar_thickness / 2.0),
                        (center, -bar_thickness / 2.0 - tooth_height),
                        (center + tooth_width / 2.0, -bar_thickness / 2.0),
                    ],
                    closed=True,
                    facecolor=blue,
                    edgecolor=edge,
                    linewidth=1.0,
                    zorder=2,
                )
            )

        # A: tooth height.
        dimension_x = centers[-1] + tooth_width * 0.75
        axis.annotate(
            "",
            xy=(dimension_x, bar_thickness / 2.0 + tooth_height),
            xytext=(dimension_x, bar_thickness / 2.0),
            arrowprops=dict(
                arrowstyle="<->", color="#7db1ff", linewidth=1.4
            ),
        )
        axis.text(
            dimension_x + total_length * 0.025,
            bar_thickness / 2.0 + tooth_height / 2.0,
            "A = {:.3g} mm".format(tooth_height),
            color="#7db1ff",
            va="center",
            fontsize=8,
        )

        # B: base width of one triangular tooth.
        width_center = centers[-2]
        width_y = bar_thickness / 2.0 + tooth_height * 0.12
        axis.annotate(
            "",
            xy=(width_center + tooth_width / 2.0, width_y),
            xytext=(width_center - tooth_width / 2.0, width_y),
            arrowprops=dict(
                arrowstyle="<->", color="#ff5c5c", linewidth=1.4
            ),
        )
        axis.text(
            width_center,
            width_y - tooth_height * 0.11,
            "B = {:.3g} mm".format(tooth_width),
            color="#ff7373",
            ha="center",
            va="top",
            fontsize=8,
        )

        # C: thickness of the central rectangular bar.
        thickness_x = right + total_length * 0.055
        axis.annotate(
            "",
            xy=(thickness_x, bar_thickness / 2.0),
            xytext=(thickness_x, -bar_thickness / 2.0),
            arrowprops=dict(
                arrowstyle="<->", color="#ffd166", linewidth=1.4
            ),
        )
        axis.text(
            thickness_x + total_length * 0.025,
            0,
            "C = {:.3g} mm".format(bar_thickness),
            color="#ffd166",
            va="center",
            fontsize=8,
        )

        channel_angles = (-45, 0, 45)
        channel_index = max(0, self.combo_fourier_channel.currentIndex())
        axis.set_title(
            "{} Fourier Stop | Channel {} ({:+d} deg) | x = {:.3g} mm | "
            "Period = {:.3g} mm".format(
                self.combo_fourier_bar_id.currentText(),
                channel_index + 1,
                channel_angles[channel_index],
                parameters["x_position"],
                period,
            ),
            color="#dce7f7",
            fontsize=9,
            pad=4,
        )
        horizontal_padding = total_length * 0.16
        vertical_extent = bar_thickness / 2.0 + tooth_height
        axis.set_xlim(
            left - horizontal_padding, right + horizontal_padding * 1.55
        )
        axis.set_ylim(-vertical_extent * 1.18, vertical_extent * 1.18)
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        self.fourier_canvas.draw_idle()

    def _apply_zernike_coefficients(self):
        coefficients = []
        for row in range(self.table_zernike.rowCount()):
            item = self.table_zernike.item(row, 1)
            text = item.text().strip() if item is not None else "0"
            try:
                coefficients.append(float(text or 0.0))
            except ValueError:
                QMessageBox.warning(
                    self,
                    "Invalid Zernike coefficient",
                    "Z{} coefficient is not numeric: {}".format(row + 1, text),
                )
                return

        x, y, wavefront, metrics = evaluate_zernike_wavefront(coefficients)
        self.zernike_figure.clear()
        axis = self.zernike_figure.add_subplot(111, projection="3d")
        surface = axis.plot_surface(
            x,
            y,
            wavefront,
            cmap="coolwarm",
            rstride=3,
            cstride=3,
            linewidth=0,
            antialiased=True,
            shade=True,
        )
        axis.set(
            xlabel="Pupil x",
            ylabel="Pupil y",
            zlabel="WFE (waves)",
        )
        axis.set_title("3-D Zernike Wavefront", fontsize=10, pad=2)
        axis.set_box_aspect((1.0, 1.0, 0.52))
        axis.view_init(elev=28, azim=-55)
        axis.tick_params(labelsize=7, pad=0)
        colorbar = self.zernike_figure.colorbar(
            surface, ax=axis, shrink=0.68, pad=0.03
        )
        colorbar.ax.tick_params(labelsize=7)
        colorbar.set_label("WFE (waves)", fontsize=8)
        self.zernike_canvas.draw_idle()
        self.lbl_zernike_metrics.setText(
            "Noll Z1-Z{} | PV: {:.5f} waves | RMS: {:.5f} waves".format(
                len(coefficients), metrics["pv"], metrics["rms"]
            )
        )
