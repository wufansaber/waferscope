# utils/ui_helpers.py
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QLabel


def create_input_widget(default_text, unit):
    """Creates a unified text layout with an input field and unit label."""
    widget = QWidget()
    layout = QHBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    line_edit = QLineEdit(default_text)
    layout.addWidget(line_edit)

    if unit:
        unit_label = QLabel(unit)
        unit_label.setFixedWidth(85)
        unit_label.setStyleSheet("color: #666666;")
        layout.addWidget(unit_label)

    return widget


def get_global_stylesheet():
    """Returns professional dark/light contrasting industrial QSS."""
    return """
    QWidget { font-size: 12px; }
    QMainWindow { background-color: #F0F2F5; }
    QMenuBar {
        background-color: #FFFFFF;
        border-bottom: 1px solid #D9D9D9;
        font-size: 12px;
    }
    QMenuBar::item {
        background-color: transparent;
        padding: 4px 9px;
        color: #333333;
    }
    QMenuBar::item:selected {
        background-color: #E6F0F9;
        color: #409EFF;
    }
    QMenu {
        background-color: #FFFFFF;
        border: 1px solid #D9D9D9;
        padding: 5px 0px;
    }
    QMenu::item {
        padding: 4px 18px;
        font-size: 12px;
    }
    QMenu::item:selected {
        background-color: #409EFF;
        color: #FFFFFF;
    }
    QListWidget {
        background-color: #FFFFFF;
        border: 1px solid #D9D9D9;
        font-size: 13px; padding: 3px; outline: 0;
    }
    QListWidget::item {
        padding: 10px 8px;
        border-bottom: 1px solid #F0F0F0;
        color: #555555;
    }
    QListWidget::item:selected {
        background-color: #E6F0F9;
        color: #333333; font-weight: bold;
        border-left: 5px solid #409EFF;
    }
    QGroupBox {
        font-size: 12px; font-weight: bold;
        border: 1px solid #B0B0B0; border-radius: 4px;
        margin-top: 10px; background-color: #FFFFFF;
    }
    QGroupBox::title {
        subcontrol-origin: margin; subcontrol-position: top left;
        padding: 0 4px; left: 8px; color: #333333;
    }
    QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
        border: 1px solid #C0C4CC; border-radius: 3px;
        padding: 3px 6px; background-color: #FAFAFA; font-weight: normal;
        min-height: 20px;
    }
    QTableWidget {
        border: 1px solid #C0C4CC; background-color: #FFFFFF; font-weight: normal;
    }
    QHeaderView::section {
        background-color: #F5F7FA; border: 1px solid #E4E7ED;
        font-weight: bold; padding: 3px;
    }
    QPushButton {
        font-weight: bold; border: none; border-radius: 3px; padding: 5px 9px;
    }
    QCheckBox { spacing: 4px; font-size: 11px; }
    QTabBar::tab { padding: 5px 10px; font-size: 11px; }
    QPushButton:hover { opacity: 0.85; }
    QScrollArea { border: none; background-color: transparent; }
    """
