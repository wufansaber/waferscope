"""Check a freshly solved camera-matched example and export a compact comparison."""
import argparse
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters
from fullwave_prototype.imaging_audit import check_saved
from tools.particle100_camera import camera_case, independent_budget, monte_carlo, direct_pixels
from tools.validate_particle100 import load_image


def plot_comparison(output, arrays, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    settings = report['settings']
    detector = settings['detector']
    lx = detector['pixels_x']*detector['pixel_pitch_um']/2
    ly = detector['pixels_y']*detector['pixel_pitch_um']/2
    maximum = max(float(arrays['ch{}_image_{}_w_um2'.format(i, name)].max())
                  for i in (1, 2, 3) for name in ('with', 'without'))
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    for row, name in enumerate(('with', 'without')):
        for i in (1, 2, 3):
            axis = arrays['ch{}_image_axis_um'.format(i)]
            x, y = abs(axis) <= lx+abs(axis[1]-axis[0]), abs(axis) <= ly+abs(axis[1]-axis[0])
            values = arrays['ch{}_image_{}_w_um2'.format(i, name)][np.ix_(y, x)]
            im = axes[row, i-1].pcolormesh(axis[x], axis[y], values,
                norm=PowerNorm(.5, vmin=0, vmax=maximum), cmap='inferno', shading='auto')
            axes[row, i-1].set(xlim=(-lx, lx), ylim=(-ly, ly), xlabel='Image X (um)',
                ylabel='Image Y (um)', aspect='equal', title='Ch{} | {}'.format(i,
                    '100 nm particle + haze' if name == 'with' else 'Haze only'))
    fig.colorbar(im, ax=list(axes.ravel()), label='W/um^2 | one absolute scale (sqrt colors)')
    fig.suptitle('Full camera field | {} nm | {} W | M={} | NA={}\n'
                 'Single coherent surface | optical transmission {}'.format(
                     report['scattering_parameters']['wavelength_nm'], settings['illumination_power_w'],
                     settings['magnification'], settings['numerical_apertures'], settings['transmission']))
    fig.savefig(output/'camera_field_comparison.png', dpi=160)
    plt.close(fig)


def ui_smoke(source, output):
    from PyQt5.QtWidgets import QApplication, QScrollArea
    from PyQt5.QtGui import QFont, QFontDatabase
    from fullwave_ui.widget import FullwaveWidget
    app = QApplication.instance() or QApplication([])
    font = Path('C:/Windows/Fonts/segoeui.ttf')
    if font.is_file():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    widget = FullwaveWidget()
    widget.resize(1440, 1000)
    if not widget.load_hybrid_result(source):
        raise ValueError('The fresh hybrid result could not be loaded by the UI')
    if widget.controls.fields['diameter-nm'].value() != 100:
        raise ValueError('The result loaded but the particle input was not restored')
    widget.show()
    widget.results.setCurrentWidget(widget.results.camera_tabs)
    widget.results.camera_tabs.setCurrentWidget(widget.results.camera_components)
    widget.findChildren(QScrollArea)[0].ensureWidgetVisible(widget.controls.field_coverage.extent)
    app.processEvents()
    widget.grab().save(str(output/'ui_camera_field.png'))
    widget.close()
    return dict(loaded=True, fit=widget.controls.fields['illumination-fit'].currentText())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default='build/field_coverage_20261016/result')
    parser.add_argument('--output', default='build/field_coverage_20261016/validation')
    parser.add_argument('--ui-smoke', action='store_true')
    args = parser.parse_args()
    source, output = Path(args.source), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    arrays, report = load_image(source)
    wavelength = report['scattering_parameters']['wavelength_nm']
    p = DetectorParameters(**report['settings']['detector'])
    camera, detectors = camera_case(arrays, report, p)
    propagation, _ = check_saved(source)
    checks = [dict(channel=i+1, budget=independent_budget(d, p, wavelength),
                   monte_carlo=monte_carlo(d, p, seed=168+i),
                   pixel_integration=direct_pixels(arrays, report, d, i, p),
                   haze_pixels_positive=int((d['haze_signal_e'] > 0).sum()),
                   haze_pixels_total=int(d['haze_signal_e'].size)) for i, d in enumerate(detectors)]
    result = dict(source=str(source.resolve()), field_coverage=report['field_coverage'],
        propagation=propagation, camera=camera, checks=checks,
        source_convergence_certified=report['surface_convergence']['full_model_converged'],
        notes=['Same-surface conditional SNR, not measured detection probability.',
               'One shared synthetic surface, not a 64-surface particle/FEM ensemble.',
               'Gaussian field construction unchanged; beam dimensions and rough extent changed.'])
    if args.ui_smoke:
        result['ui'] = ui_smoke(source, output)
    plot_comparison(output, arrays, report)
    (output/'validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    passed = (report['field_coverage']['all_channels_covered'] and
              all(c['passed'] for c in propagation['channels']) and
              all(c['monte_carlo']['passed'] and c['budget']['variance_relative_max_error'] < 1e-12 for c in checks))
    print(json.dumps(dict(passed=passed, channels=[{k: d[k] for k in
        ('detection_roi_snr', 'detection_template_snr', 'saturated_pixels')} for d in detectors]), indent=2))
    if not passed:
        raise RuntimeError('Camera-field verification failed; inspect validation.json')


if __name__ == '__main__':
    main()
