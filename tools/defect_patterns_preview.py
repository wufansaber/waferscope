"""Common-scale geometry/near-field/image figure and saved-result GUI checks."""
import argparse
import json
import os
from pathlib import Path

import numpy as np


def compare(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Polygon
    from fullwave_prototype.defects import footprint_nm, draw_section
    cases = []
    for shape in ('line_x', 'line_y', 'cross'):
        folder = root/shape
        report = json.loads((folder/'scattering/report.json').read_text(encoding='utf-8'))
        with np.load(folder/'scattering/paired_fields.npz') as saved:
            data = {k: saved[k] for k in ('near_xy_u_um', 'near_xy_v_um', 'near_xy_position_um', 'near_xy_e_with')}
        with np.load(folder/'image_fields.npz') as saved:
            data.update(image=saved['ch3_image_with_w_um2'], axis=saved['ch3_image_axis_um'])
        cases.append((shape, report['parameters'], data))
    maximum = max(np.sum(abs(data['near_xy_e_with'])**2, axis=-1).max() for _, _, data in cases)
    image_max = max(data['image'].max() for _, _, data in cases)
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), constrained_layout=True)
    for j, (name, p, data) in enumerate(cases):
        ax = axes[0, j]
        ax.add_patch(Polygon(footprint_nm(p), facecolor='#d78ba6', edgecolor='#812d59'))
        ax.set(title=name+' | XY geometry', xlim=(-130, 130), ylim=(-130, 130),
               xlabel='X (nm)', ylabel='Y (nm)', aspect='equal')
        value = np.sum(abs(data['near_xy_e_with'])**2, axis=-1)
        near = axes[1, j].pcolormesh(data['near_xy_u_um']*1000, data['near_xy_v_um']*1000,
            value, shading='auto', cmap='magma', norm=PowerNorm(.5, vmin=0, vmax=maximum))
        draw_section(axes[1, j], p, 'xy', float(data['near_xy_position_um'])*1000)
        axes[1, j].set(title='FEM near XY | mid-height', xlabel='X (nm)', ylabel='Y (nm)', aspect='equal')
        im = axes[2, j].pcolormesh(data['axis'], data['axis'], data['image'], shading='auto',
            cmap='inferno', norm=PowerNorm(.5, vmin=0, vmax=image_max))
        axes[2, j].set(title='CH3 image | 20x, NA 0.5', xlabel='Image x (um)', ylabel='Image y (um)',
                       xlim=(-12, 12), ylim=(-12, 12), aspect='equal')
    fig.colorbar(near, ax=list(axes[1]), label='|E / E_inc|^2 | shared sqrt colors', shrink=.8)
    fig.colorbar(im, ax=list(axes[2]), label='W / um^2 | shared sqrt colors', shrink=.8)
    fig.suptitle('FEM PREVIEW | 266 nm, p, 75 deg from -Y | L 200, W 50, H 50, gap 10 nm\n'
                 'Smooth wafer | 1 W, plane-wave scaling over 2 x 2 um | transmission 1e-6')
    fig.savefig(root/'pattern_comparison.png', dpi=140)
    plt.close(fig)


def ui_check(root, solve=False):
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['FARFIELD_RESULTS_DIR'] = str((root/'ui_results').resolve())
    from PyQt5.QtCore import QTimer, QProcess
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication
    from scatmech_ui.app import ScatmechStandaloneWindow
    app = QApplication.instance() or QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    window = ScatmechStandaloneWindow()
    window.tabs.setCurrentIndex(1)
    widget = window.fullwave_widget
    for shape in ('line_x', 'line_y', 'cross'):
        if not widget.load_result(root/shape):
            raise RuntimeError('Cannot open '+shape+' in the UI')
        f = widget.controls.fields
        if (f['defect-shape'].currentText() != shape or f['defect-length-nm'].value() != 200
                or f['defect-width-nm'].value() != 50 or f['defect-height-nm'].value() != 50):
            raise RuntimeError('Saved pattern geometry was not restored')
    for width, height in ((1680, 980), (1200, 800)):
        window.resize(width, height)
        window.show()
        app.processEvents()
        widget.controls.fields['defect-shape'].setCurrentText('cross')
        widget.results.near_tabs.setCurrentIndex(2)
        if not window.grab().save(str(root/('pattern_ui_{}.png'.format(width)))):
            raise RuntimeError('Cannot capture pattern UI')
    result = dict(status='passed', loaded_shapes=['line_x', 'line_y', 'cross'], sizes=[1680, 1200])
    if solve:
        widget.controls.preset.setCurrentText('266 nm cross preview')
        widget.process.finished.connect(app.quit)
        watchdog = QTimer()
        watchdog.setSingleShot(True)
        watchdog.timeout.connect(app.quit)
        watchdog.start(180000)
        QTimer.singleShot(0, widget.start.click)
        app.exec_()
        watchdog.stop()
        if not widget.wait_for_shutdown() or widget.output is None:
            raise RuntimeError('UI pattern solver did not finish')
        resources = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
        report = json.loads((widget.output/'imaging_report.json').read_text(encoding='utf-8'))
        if (resources['status'] != 'completed' or resources['surviving_owned_pids']
                or widget.process.exitStatus() != QProcess.NormalExit or widget.process.exitCode()
                or report['scattering_parameters']['defect_shape'] != 'cross'
                or not all(ch['detector']['gray_valid'] for ch in report['channels'])):
            raise RuntimeError('UI pattern solve or camera acquisition failed')
        result['ui_calculate'] = dict(status='passed', result=str(widget.output),
            gray_snr=[ch['detector']['gray_snr'] for ch in report['channels']])
    window.close()
    app.processEvents()
    (root/'pattern_ui_check.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('build/defect_patterns_20261017/validated'))
    parser.add_argument('--solve-ui', action='store_true')
    args = parser.parse_args()
    compare(args.root)
    ui_check(args.root, args.solve_ui)


if __name__ == '__main__':
    main()
