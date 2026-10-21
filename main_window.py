# main_window.py
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QListWidget,
    QStackedWidget, QAction, QMessageBox
)
from utils.ui_helpers import get_global_stylesheet

# Import page modules
from pages.page_basic import BasicSettingsPage
from pages.page_wafer import WaferPatternPage
from pages.page_advanced import AdvancedSettingsPage
from pages.page_results import ResultsPage
from scatmech_ui.widgets import ScatmechSimulationWidget


class WaferScopeMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("WaferScope")
        self.configure_window_size()

        # ---------------- Create Menu Bar (商业通用主菜单栏栏) ----------------
        self.init_menu_bar()

        # ---------------- Base Layout ----------------
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(6, 6, 6, 6)
        main_layout.setSpacing(8)

        # ---------------- Left Navigation Bar ----------------
        self.nav_list = QListWidget()
        self.nav_list.setFixedWidth(180)
        nav_items = [
            "Basic Settings",
            "Bare Wafer & Pattern",
            "Advanced Settings",
            "Simulation Results",
            "Far-Field"
        ]
        self.nav_list.addItems(nav_items)
        main_layout.addWidget(self.nav_list)

        # ---------------- Right Stacked Widget ----------------
        self.stacked_widget = QStackedWidget()
        main_layout.addWidget(self.stacked_widget)

        # Instantiate Pages
        self.page_basic = BasicSettingsPage()
        self.page_wafer = WaferPatternPage()
        self.page_advanced = AdvancedSettingsPage()
        self.page_results = ResultsPage()
        self.page_scatmech = ScatmechSimulationWidget()

        # Add to stack
        self.stacked_widget.addWidget(self.page_basic)
        self.stacked_widget.addWidget(self.page_wafer)
        self.stacked_widget.addWidget(self.page_advanced)
        self.stacked_widget.addWidget(self.page_results)
        self.stacked_widget.addWidget(self.page_scatmech)

        # Routing connections
        self.nav_list.currentRowChanged.connect(self.stacked_widget.setCurrentIndex)
        self.nav_list.setCurrentRow(0)

        # Apply Global QSS Stylesheet
        self.setStyleSheet(get_global_stylesheet())

    def configure_window_size(self):
        """Scale and center the application window for Full HD and 2K displays."""
        screen = QApplication.primaryScreen()
        if screen is None:
            self.resize(1450, 950)
            return

        available = screen.availableGeometry()
        width = min(1680, max(960, int(available.width() * 0.92)))
        height = min(1050, max(620, int(available.height() * 0.90)))
        width = min(width, available.width())
        height = min(height, available.height())

        self.setMinimumSize(
            min(900, available.width()),
            min(600, available.height())
        )
        self.resize(width, height)
        self.move(
            available.x() + (available.width() - width) // 2,
            available.y() + (available.height() - height) // 2
        )

    def init_menu_bar(self):
        """Initializes standard EDA software menus with professional actions."""
        menu_bar = self.menuBar()

        # 1. File Menu
        file_menu = menu_bar.addMenu("&File")

        new_proj_action = QAction("&New Project", self)
        open_proj_action = QAction("&Open Project...", self)
        save_param_action = QAction("&Save Parameters", self)
        import_param_action = QAction("&Import Parameters...", self)
        exit_action = QAction("&Exit", self)
        exit_action.triggered.connect(self.close)

        file_menu.addAction(new_proj_action)
        file_menu.addAction(open_proj_action)
        file_menu.addSeparator()
        file_menu.addAction(save_param_action)
        file_menu.addAction(import_param_action)
        file_menu.addSeparator()
        file_menu.addAction(exit_action)

        # 2. Settings Menu
        settings_menu = menu_bar.addMenu("&Settings")
        pref_action = QAction("&Preferences...", self)
        engine_config_action = QAction("&Core Engine Config", self)
        settings_menu.addAction(pref_action)
        settings_menu.addAction(engine_config_action)

        # 3. Help Menu
        help_menu = menu_bar.addMenu("&Help")
        doc_action = QAction("&Documentation (F1)", self)
        about_action = QAction("&About Simulation Platform", self)
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(doc_action)
        help_menu.addAction(about_action)

    def show_about_dialog(self):
        """Triggers standard informative industrial software dialog."""
        QMessageBox.about(
            self,
            "About WaferScope",
            "WaferScope - Wafer Defect Optical Simulation\n"
            "Version 1.7 (Modular Release)\n\n"
            "Beijing Opto-Tech Co., Ltd. All Rights Reserved."
        )

    def closeEvent(self, event):
        self.page_scatmech.wait_for_shutdown()
        super().closeEvent(event)


# Retain imports used by existing launchers and notebooks.
APTKSimulationUI = WaferScopeMainWindow
