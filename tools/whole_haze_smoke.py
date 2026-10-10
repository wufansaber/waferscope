"""Prepare and visually check the UI-configured whole-beam haze calculation."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.widget import FullwaveWidget
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    from fullwave_prototype.imaging_options import imaging_settings
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--result')
    parser.add_argument('--samples', type=int, default=16)
    parser.add_argument('--psd-model', choices=('gaussian', 'abc'), default='gaussian')
    parser.add_argument('--psd-length', type=float, default=.05)
    parser.add_argument('--psd-exponent', type=float, default=3.)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    widget = FullwaveWidget()
    widget.resize(1480, 900)
    widget.controls.fields['wavelength-nm'].setValue(266)
    widget.controls.fields['roughness-model'].setCurrentText(args.psd_model)
    widget.controls.fields['correlation-length-um'].setValue(args.psd_length)
    widget.controls.fields['psd-exponent'].setValue(args.psd_exponent)
    widget.controls.fields['pupil-samples'].setCurrentText('128')
    widget.controls.fields['padding'].setCurrentText('2')
    widget.ensemble_count.setValue(args.samples)
    parsed = build_parser().parse_args(widget.controls.arguments(particle_enabled=False))
    p, rough = scattering_settings(parsed, particle_enabled=False)
    settings = imaging_settings(parsed, parameters=p)
    inputs = dict(parameters=asdict(p), roughness=asdict(rough), imaging_settings=asdict(settings),
                  samples=args.samples, coverage=.99, variance_fraction=.9, taper_width_um=.025)
    (output/'whole_haze_input.json').write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    widget.show()
    app.processEvents()
    if args.result:
        assert widget.load_whole_haze_result(args.result)
        app.processEvents()
        widget.results.fit_current()
        assert widget.grab().save(str(output/'whole_haze_ui.png'))
        widget.results.whole_haze_tabs.setCurrentWidget(widget.results.whole_haze_views['whole_haze_far.png'])
        app.processEvents()
        assert widget.grab().save(str(output/'whole_haze_far_ui.png'))
    print('Whole-beam haze UI input/screenshot check passed')


if __name__ == '__main__':
    main()
