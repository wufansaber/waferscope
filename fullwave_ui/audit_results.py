"""Read-only provenance-aware viewer for saved symmetry and NIST validation."""
import json
from pathlib import Path

from PyQt5.QtWidgets import QPlainTextEdit, QTabWidget


class AuditResults(QTabWidget):
    def __init__(self):
        super().__init__()
        from .results import FigureView
        self.figures = {}
        for filename, title in (('nist_haze_comparison.png', 'NIST PSD'),
            ('low_haze_comparison.png', 'Low roughness'), ('peak_decomposition.png', 'Peak / interference'),
            ('rough_mirror.png', 'Mirrored surface'), ('ensemble_confidence.png', 'Confidence intervals')):
            view = FigureView()
            self.figures[filename] = view
            self.addTab(view, title)
        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.addTab(self.details, 'Saved inputs / limits')

    def load(self, directory):
        directory = Path(directory)
        filename = 'summary_report.json' if (directory/'summary_report.json').exists() else 'nist_haze_report.json'
        report = json.loads((directory/filename).read_text(encoding='utf-8'))
        for name, view in self.figures.items():
            available = (directory/name).is_file()
            if available:
                view.load(directory/name)
            else:
                view.scene().clear()
                view.item = None
            self.setTabEnabled(self.indexOf(view), available)
        self.details.setPlainText('Saved validation: '+str(directory.resolve())+'\n\n'+
                                  json.dumps(report, indent=2, ensure_ascii=False))
        self.setCurrentIndex(next((self.indexOf(v) for v in self.figures.values() if v.item is not None),
                                 self.indexOf(self.details)))
