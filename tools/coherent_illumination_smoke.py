"""Check finite-beam controls and real saved results in offscreen Qt."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication, QScrollArea
    from fullwave_ui.controls import SimulationControls
    from fullwave_ui.results import ResultsTabs

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline', required=True)
    parser.add_argument('--ensemble', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    controls = SimulationControls()
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setWidget(controls)
    scroll.resize(420, 720)
    scroll.show()
    app.processEvents()
    scroll.ensureWidgetVisible(controls.fields['beam-profile'], 0, 160)
    app.processEvents()
    assert controls.fields['beam-profile'].currentText() == 'gaussian'
    assert controls.fields['beam-order'].currentText() == '31'
    assert not controls.fields['beam-waist-um'].isEnabled()
    assert scroll.grab().save(str(output/'illumination_controls.png'))

    results = ResultsTabs()
    results.resize(1100, 800)
    results.show()
    for kind, folder in (('pipeline', args.pipeline), ('ensemble', args.ensemble)):
        if kind == 'pipeline':
            results.load(folder)
            assert all(view.item is not None for view in results.near_views.values())
            assert all(view.item is not None for view in results.image_views.values())
        else:
            results.load_ensemble(folder)
        results.setCurrentWidget(results.illumination_tabs)
        results.illumination_tabs.setCurrentWidget(results.illumination_view)
        app.processEvents()
        results.fit_current()
        assert results.illumination_view.item is not None
        assert 'normalized_incident_power_w' in results.illumination_details.toPlainText()
        assert results.grab().save(str(output/(kind+'_illumination.png')))
    report = json.loads((Path(args.ensemble)/'haze_ensemble_report.json').read_text(encoding='utf-8'))
    assert report['metrics']['samples'] == 2
    assert not report['metrics']['area_extrapolation_applied']
    assert report['illumination']['model'] == 'gaussian'
    results.setCurrentWidget(results.ensemble_tabs)
    results.ensemble_tabs.setCurrentWidget(results.ensemble_view)
    app.processEvents()
    results.fit_current()
    assert results.grab().save(str(output/'finite_haze_ensemble.png'))
    print('Finite coherent illumination controls, pipeline and ensemble UI smoke passed')


if __name__ == '__main__':
    main()
