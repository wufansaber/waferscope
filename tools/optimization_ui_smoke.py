"""Offscreen UI check of the polynomial-order row and saved convergence results."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication, QDialog, QScrollArea
    from fullwave_ui.convergence_dialog import ConvergenceDialog
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--result')
    parser.add_argument('--default-scans', action='store_true')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    plan = json.loads(Path(args.plan).read_text(encoding='utf-8'))
    if args.default_scans:
        from fullwave_prototype.haze_convergence import default_plan
        from fullwave_prototype.parameters import Parameters
        from fullwave_prototype.rough_surface import Roughness
        plan = default_plan(Parameters(**plan['parameters']), Roughness(**plan['roughness']),
                            plan['samples'], plan['numerical_apertures'], plan['tilt_degrees'])
    dialog = ConvergenceDialog(plan)
    dialog.show()
    app.processEvents()
    if set(dialog.levels) != set(plan['scans']):
        raise RuntimeError('Convergence rows missing')
    if dialog.height() > app.primaryScreen().availableGeometry().height():
        raise RuntimeError('Dialog exceeds available screen height')
    if not dialog.grab().save(str(output/'convergence_dialog.png')):
        raise RuntimeError('Could not save dialog screenshot')
    scroll = dialog.findChild(QScrollArea)
    scroll.verticalScrollBar().setValue(scroll.verticalScrollBar().maximum())
    app.processEvents()
    if not dialog.grab().save(str(output/'convergence_dialog_bottom.png')):
        raise RuntimeError('Could not save scrolled dialog screenshot')
    dialog.accept()
    if dialog.result() != QDialog.Accepted:
        raise RuntimeError('UI plan validation failed')
    if args.result:
        from fullwave_ui.results import ResultsTabs
        results = ResultsTabs()
        results.resize(1050, 760)
        results.load_convergence(args.result)
        results.setCurrentWidget(results.convergence_tabs)
        results.show()
        app.processEvents()
        if not results.grab().save(str(output/'convergence_result.png')):
            raise RuntimeError('Could not save results screenshot')
        if 'paired_same_realization' not in results.convergence_details.toPlainText():
            raise RuntimeError('Coupled refinement metadata missing from UI')
        if (Path(args.result)/'haze_paired_errors.png').is_file():
            if results.paired_error_view.item is None:
                raise RuntimeError('Paired numerical error figure missing')
            results.convergence_tabs.setCurrentWidget(results.paired_error_view)
            app.processEvents()
            if not results.grab().save(str(output/'convergence_paired_errors.png')):
                raise RuntimeError('Could not save paired error screenshot')
    print('Convergence dialog/result smoke passed')


if __name__ == '__main__':
    main()
