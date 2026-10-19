"""Prepare small UI hybrid inputs and check result loading, without solving inside Qt."""
import argparse
from dataclasses import asdict, replace
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
    parser.add_argument('--power', type=float)
    parser.add_argument('--disable-channel', type=int, choices=(1, 2, 3))
    parser.add_argument('--validation', nargs='+')
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    widget = FullwaveWidget()
    widget.resize(1480, 900)
    controls = widget.controls.fields
    controls['diameter-nm'].setValue(100)
    controls['particle-mesh-um'].setValue(0)
    controls['wavelength-nm'].setValue(266)
    controls['pupil-samples'].setCurrentText('32')
    controls['padding'].setCurrentText('2')
    parsed = build_parser().parse_args(widget.controls.arguments())
    p, r = scattering_settings(parsed)
    r = replace(r, taper_fraction=.125)
    inputs = dict(parameters=asdict(p), roughness=asdict(r),
        imaging_settings=asdict(imaging_settings(parsed, parameters=p)),
        coverage=.99, variance_fraction=.9, taper_width_um=.025, quadrature_order=5)
    (out/'hybrid_input.json').write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    widget.show()
    app.processEvents()
    if args.result:
        assert widget.load_result(args.result)
        assert widget.loaded_hybrid is not None
        assert widget.hybrid_reimage.isEnabled()
        assert not widget.image_current_button.isEnabled()
        assert not widget.compare_button.isEnabled()
        assert widget.validation_button.isEnabled()
        from fullwave_prototype.hybrid_run import load_sources
        from fullwave_prototype.parameters import Parameters
        _, metadata = load_sources(widget.loaded_hybrid)
        saved_p = Parameters(**metadata['parameters'])
        if args.power is not None:
            controls['illumination-power-w'].setValue(args.power)
        if args.disable_channel:
            controls['na'+str(args.disable_channel)].setValue(0)
        optical_args = build_parser().parse_args(widget.controls.arguments(
            imaging_only=True, parameters=saved_p, incidence_deg=saved_p.incidence_deg))
        optical = imaging_settings(optical_args, parameters=saved_p)
        (out/'reimage_settings.json').write_text(json.dumps(asdict(optical), indent=2), encoding='utf-8')
        app.processEvents()
        widget.results.fit_current()
        assert widget.grab().save(str(out/'hybrid_ui.png'))
    else:
        assert widget.hybrid_button.isEnabled()
        assert widget.grab().save(str(out/'hybrid_input_ui.png'))
    for index, directory in enumerate(args.validation or []):
        source_before = widget.loaded_hybrid
        assert widget.load_model_validation(directory)
        assert widget.loaded_hybrid == source_before
        widget.results.validation_tabs.setCurrentWidget(widget.results.validation_view)
        app.processEvents()
        widget.results.fit_current()
        assert widget.grab().save(str(out/'validation_ui_{}.png'.format(index+1)))
    print('Hybrid UI inputs/result routing passed', flush=True)
    widget.close()


if __name__ == '__main__':
    main()
