"""Aspect-preserving scientific result views and saved-parameter provenance."""
import json
import math
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QGraphicsScene, QGraphicsView, QPlainTextEdit,
                            QTabWidget)


class FigureView(QGraphicsView):
    def __init__(self):
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setToolTip('Result figure')
        self.item = None

    def load(self, path):
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            raise ValueError('Cannot load figure: '+str(path))
        self.scene().clear()
        self.item = self.scene().addPixmap(pixmap)
        self.item.setTransformationMode(Qt.SmoothTransformation)
        self.scene().setSceneRect(self.item.boundingRect())
        self.fit()

    def fit(self):
        if self.item is not None:
            self.fitInView(self.item, Qt.KeepAspectRatio)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.fit()

    def wheelEvent(self, event):
        self.scale(1.2 if event.angleDelta().y() > 0 else 1/1.2,
                   1.2 if event.angleDelta().y() > 0 else 1/1.2)
        event.accept()


class ResultsTabs(QTabWidget):
    FIGURES = (('Near field', 'near_field.png'), ('Far field', 'far_field.png'),
               ('Three image planes', 'image_planes.png'),
               ('Coherent decomposition', 'image_decomposition.png'))

    def __init__(self):
        super().__init__()
        self.views = []
        self.validation_tabs = QTabWidget()
        self.validation_view = FigureView()
        self.validation_details = QPlainTextEdit()
        self.validation_details.setReadOnly(True)
        self.validation_tabs.addTab(self.validation_view, 'Diagnostics')
        self.validation_tabs.addTab(self.validation_details, 'Inputs / limits')
        self.validation_tabs.currentChanged.connect(self.fit_current)
        for title, _ in self.FIGURES:
            view = FigureView()
            self.views.append(view)
            if not self.count():
                self.near_tabs = QTabWidget()
                self.near_tabs.addTab(view, 'Overview')
                self.near_views = {}
                for plane in ('xz', 'xy', 'yz'):
                    self.near_views[plane] = FigureView()
                    index = self.near_tabs.addTab(self.near_views[plane], plane.upper())
                    self.near_tabs.setTabEnabled(index, False)
                self.near_tabs.currentChanged.connect(self.fit_current)
                self.addTab(self.near_tabs, title)
            elif title == 'Three image planes':
                self.image_tabs = QTabWidget()
                self.image_tabs.addTab(view, 'Overview')
                self.image_views = {}
                for channel in (1, 2, 3):
                    self.image_views[channel] = FigureView()
                    index = self.image_tabs.addTab(self.image_views[channel], 'Ch{}'.format(channel))
                    self.image_tabs.setTabEnabled(index, False)
                self.image_tabs.currentChanged.connect(self.fit_current)
                self.addTab(self.image_tabs, title)
            else:
                self.addTab(view, title)
        self.comparison_tabs = QTabWidget()
        self.comparison_views = {}
        for key, title in (('far', 'Far field'), ('far_log', 'Far log'),
                           ('image', 'Image plane'), ('image_log', 'Image log')):
            view = FigureView()
            self.comparison_views[key] = view
            self.comparison_tabs.addTab(view, title)
        self.comparison_totals = QPlainTextEdit()
        self.comparison_totals.setReadOnly(True)
        self.comparison_tabs.addTab(self.comparison_totals, 'Totals')
        self.addTab(self.comparison_tabs, 'Particle / haze')
        self.setTabEnabled(self.indexOf(self.comparison_tabs), False)
        self.ensemble_tabs = QTabWidget()
        self.ensemble_view = FigureView()
        self.ensemble_tabs.addTab(self.ensemble_view, 'Overview')
        self.ensemble_figures = {}
        for key, title in (('fem', 'FEM mean'), ('finite', 'Finite patch PSD'),
                           ('nist', 'NIST PSD'), ('stderr', 'Std. error')):
            view = FigureView()
            self.ensemble_figures[key] = view
            self.ensemble_tabs.addTab(view, title)
        self.ensemble_metrics = QPlainTextEdit()
        self.ensemble_metrics.setReadOnly(True)
        self.ensemble_tabs.addTab(self.ensemble_metrics, 'Statistics')
        self.sample_convergence_view = FigureView()
        self.ensemble_tabs.addTab(self.sample_convergence_view, 'Sample convergence')
        self.addTab(self.ensemble_tabs, 'Haze ensemble')
        self.setTabEnabled(self.indexOf(self.ensemble_tabs), False)
        self.convergence_tabs = QTabWidget()
        self.illumination_view = FigureView()
        self.illumination_details = QPlainTextEdit()
        self.illumination_details.setReadOnly(True)
        self.illumination_tabs = QTabWidget()
        self.illumination_tabs.addTab(self.illumination_view, 'Incident footprint')
        self.illumination_tabs.addTab(self.illumination_details, 'Power / coverage')
        self.addTab(self.illumination_tabs, 'Illumination')
        self.setTabEnabled(self.indexOf(self.illumination_tabs), False)
        self.whole_haze_tabs = QTabWidget()
        self.whole_haze_views = {}
        for filename, title in (('whole_surface.png', 'Surface / illumination'),
                                ('whole_haze_far.png', 'Mean far field'),
                                ('whole_haze_images.png', 'Mean image planes'),
                                ('whole_haze_statistics.png', 'Sample convergence')):
            view = FigureView()
            self.whole_haze_views[filename] = view
            self.whole_haze_tabs.addTab(view, title)
        self.whole_haze_details = QPlainTextEdit()
        self.whole_haze_details.setReadOnly(True)
        self.whole_haze_tabs.addTab(self.whole_haze_details, 'Model / coverage / power')
        self.whole_haze_tabs.currentChanged.connect(self.fit_current)
        self.addTab(self.whole_haze_tabs, 'Whole-beam haze (SPM)')
        self.setTabEnabled(self.indexOf(self.whole_haze_tabs), False)
        self.convergence_view = FigureView()
        self.paired_error_view = FigureView()
        self.convergence_details = QPlainTextEdit()
        self.convergence_details.setReadOnly(True)
        self.convergence_tabs.addTab(self.convergence_view, 'Energy / scan levels')
        self.convergence_tabs.addTab(self.paired_error_view, 'Paired numerical errors')
        self.convergence_tabs.addTab(self.convergence_details, 'Statistics / resources')
        self.addTab(self.convergence_tabs, 'Convergence')
        self.setTabEnabled(self.indexOf(self.convergence_tabs), False)
        self.parameters = QPlainTextEdit()
        self.parameters.setReadOnly(True)
        self.addTab(self.parameters, 'Result parameters')
        self.detector_view = FigureView()
        self.addTab(self.detector_view, 'Camera')
        self.setTabEnabled(self.indexOf(self.detector_view), False)
        self.currentChanged.connect(self.fit_current)

    def fit_current(self, *_):
        view = self.currentWidget()
        if isinstance(view, QTabWidget):
            view = view.currentWidget()
        if isinstance(view, FigureView):
            view.fit()

    def load_whole_haze(self, directory):
        directory = Path(directory)
        report = json.loads((directory/'whole_haze_report.json').read_text(encoding='utf-8'))
        if report.get('model') != 'whole_beam_first_order_surface' or report.get('status') != 'complete':
            raise ValueError('Select a complete whole-beam first-order surface result')
        for filename in self.whole_haze_views:
            if QPixmap(str(directory/filename)).isNull():
                raise ValueError('Missing whole-beam figure: '+filename)
        for filename, view in self.whole_haze_views.items():
            view.load(directory/filename)
        self.whole_haze_details.setPlainText(json.dumps(report, indent=2))
        self.setTabEnabled(self.indexOf(self.whole_haze_tabs), True)
        return report

    def load(self, directory):
        directory = Path(directory)
        report = json.loads((directory/'imaging_report.json').read_text(encoding='utf-8'))
        channels = report['channels']
        if len(channels) != 3:
            raise ValueError('Result must contain three channels')
        powers = [float(c['with']['image_power_w'])*1e12 for c in channels]
        if not all(math.isfinite(power) and power >= 0 for power in powers):
            raise ValueError('Result contains invalid image powers')
        # Validate every image before replacing the current result.
        detector = report.get('settings', {}).get('detector', {}).get('enabled', False)
        if detector and QPixmap(str(directory/'detector.png')).isNull():
            raise ValueError('Missing or invalid detector figure')
        for _, name in self.FIGURES:
            if QPixmap(str(directory/name)).isNull():
                raise ValueError('Missing or invalid result figure: '+name)
        planes = {info['plane'].lower() for info in report.get('near_field_planes', [])}
        if not planes.issubset(self.near_views):
            raise ValueError('Unknown near-field plane')
        for plane in planes:
            if QPixmap(str(directory/('near_field_'+plane+'.png'))).isNull():
                raise ValueError('Missing or invalid near-field figure: '+plane)
        image_channels = report.get('image_channel_figures', [])
        if image_channels not in ([], [1, 2, 3]):
            raise ValueError('Invalid image channel figure list')
        for channel in image_channels:
            if QPixmap(str(directory/('image_channel_{}.png'.format(channel)))).isNull():
                raise ValueError('Missing or invalid image channel figure: '+str(channel))
        comparison = None
        manifest_path = directory/'component_comparison_manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
            comparison_dir = (directory/manifest['directory']).resolve()
            comparison_dir.relative_to(directory.resolve())
            comparison_report = json.loads((comparison_dir/manifest['report']).read_text(encoding='utf-8'))
            if comparison_report.get('status') != 'complete':
                raise ValueError('Incomplete component comparison')
            figures = manifest['figures']
            for key in self.comparison_views:
                name = figures[key]
                if Path(name).name != name or QPixmap(str(comparison_dir/name)).isNull():
                    raise ValueError('Missing or invalid component comparison figure: '+key)
            totals = manifest['totals']
            if Path(totals).name != totals:
                raise ValueError('Invalid component comparison totals path')
            comparison = (comparison_dir, figures,
                          (comparison_dir/totals).read_text(encoding='utf-8'))
        for view, (_, name) in zip(self.views, self.FIGURES):
            view.load(directory/name)
        self.near_tabs.setCurrentIndex(0)
        for plane, view in self.near_views.items():
            if plane in planes:
                view.load(directory/('near_field_'+plane+'.png'))
            else:
                view.scene().clear()
                view.item = None
            self.near_tabs.setTabEnabled(self.near_tabs.indexOf(view), plane in planes)
        self.image_tabs.setCurrentIndex(0)
        for channel, view in self.image_views.items():
            if channel in image_channels:
                view.load(directory/('image_channel_{}.png'.format(channel)))
            else:
                view.scene().clear()
                view.item = None
            self.image_tabs.setTabEnabled(self.image_tabs.indexOf(view), channel in image_channels)
        if comparison is not None:
            comparison_dir, figures, totals = comparison
            for key, view in self.comparison_views.items():
                view.load(comparison_dir/figures[key])
            self.comparison_totals.setPlainText(totals)
        else:
            for view in self.comparison_views.values():
                view.scene().clear()
                view.item = None
            self.comparison_totals.clear()
        self.setTabEnabled(self.indexOf(self.comparison_tabs), comparison is not None)
        if detector:
            self.detector_view.load(directory/'detector.png')
        else:
            self.detector_view.scene().clear()
            self.detector_view.item = None
        self.setTabEnabled(self.indexOf(self.detector_view), detector)
        self.parameters.setPlainText('Loaded result: '+str(directory.resolve())+'\n\n'+
                                     json.dumps(report, indent=2, ensure_ascii=False))
        self.load_illumination(directory)
        return powers

    def load_illumination(self, directory):
        directory = Path(directory)
        available = (directory/'illumination.png').is_file() and (directory/'illumination_report.json').is_file()
        if available:
            report = json.loads((directory/'illumination_report.json').read_text(encoding='utf-8'))
            self.illumination_view.load(directory/'illumination.png')
            self.illumination_details.setPlainText(json.dumps(report, indent=2))
        else:
            self.illumination_view.scene().clear()
            self.illumination_view.item = None
            self.illumination_details.clear()
        self.setTabEnabled(self.indexOf(self.illumination_tabs), available)

    def load_validation(self, directory):
        directory = Path(directory)
        if (directory/'surface_mesh_report.json').is_file():
            filename, figure = 'surface_mesh_report.json', 'surface_mesh.png'
            model, metric = 'same_parent_spatial_taper_audit', 'coverage'
        elif (directory/'range_report.json').is_file():
            filename, figure = 'range_report.json', 'range_validation.png'
            model = 'recovered_parent_outer_spm_range_validation'
            metric = 'input_power_coverage'
        else:
            filename, figure = 'validation_report.json', 'local_validation.png'
            model = 'shared_surface_low_height_fem_spm_validation'
            metric = None
        report = json.loads((directory/filename).read_text(encoding='utf-8'))
        if report.get('model') != model or report.get('status') != 'complete_with_physical_limits':
            raise ValueError('Select a completed model check; interrupted calculations are not validated results')
        if not report.get('cases'):
            raise ValueError('Validation report has no cases')
        for case in report['cases']:
            value = case[metric] if metric else case['comparison']['complex_relative_l2']
            if not math.isfinite(value) or value < 0:
                raise ValueError('Invalid validation metric')
            if metric:
                for key in ('collected_spm_power_w', 'from_baseline_cap_complex_relative_l2',
                            'from_previous_cap_complex_relative_l2', 'angular_power_relative_changes',
                            'cap_power_w', 'complex_relative_l2', 'power_relative_changes'):
                    if key in case and (len(case[key]) != 3 or not all(
                            math.isfinite(v) and v >= 0 for v in case[key])):
                        raise ValueError('Invalid validation metric: '+key)
        self.validation_view.load(directory/figure)
        self.validation_details.setPlainText(json.dumps(report, indent=2))
        if self.indexOf(self.validation_tabs) < 0:
            self.addTab(self.validation_tabs, 'Model checks')
        return report

    def load_ensemble(self, directory):
        directory = Path(directory)
        report = json.loads((directory/'haze_ensemble_report.json').read_text(encoding='utf-8'))
        if report.get('status') != 'complete':
            raise ValueError('Haze ensemble result is incomplete')
        figures = ['haze_ensemble.png', 'haze_ensemble_fem.png',
                   'haze_ensemble_stderr.png']
        if report.get('nist_reference'):
            figures.append('haze_ensemble_nist.png')
        if report.get('finite_patch_reference'):
            figures.append('haze_ensemble_finite.png')
        for name in figures:
            if QPixmap(str(directory/name)).isNull():
                raise ValueError('Missing or invalid haze ensemble figure: '+name)
        metrics = report['metrics']
        if metrics['samples'] != len(report['seeds']) or metrics['samples'] < 1:
            raise ValueError('Invalid haze ensemble sample count')
        self.ensemble_view.load(directory/'haze_ensemble.png')
        for key, view in self.ensemble_figures.items():
            name = directory/('haze_ensemble_'+key+'.png')
            if name.is_file():
                view.load(name)
            else:
                view.scene().clear()
                view.item = None
            self.ensemble_tabs.setTabEnabled(self.ensemble_tabs.indexOf(view), name.is_file())
        self.ensemble_tabs.setCurrentWidget(self.ensemble_figures['fem'])
        lines = ['Haze-only FEM ensemble',
                 'Seeds: '+', '.join(str(seed) for seed in report['seeds']),
                 'Comparison area: {:.6g} um^2'.format(metrics['comparison_area_um2']),
                 'Patch geometric area: {:.6g} um^2'.format(metrics['patch_geometric_area_um2']),
                 'Patch squared-window effective area: {:.6g} um^2'.format(
                     metrics['patch_effective_area_um2']),
                 'FEM hemisphere DSC: {:.8g} um^2'.format(metrics['fem_mean_hemisphere_dsc_um2'])]
        if metrics.get('fem_hemisphere_standard_error_um2') is not None:
            lines.append('FEM standard error: {:.8g} um^2'.format(
                metrics['fem_hemisphere_standard_error_um2']))
        if metrics.get('nist_hemisphere_dsc_um2') is not None:
            lines += ['NIST {} PSD: {:.8g} um^2'.format(
                          metrics.get('psd_model', 'Gaussian').upper(),
                          metrics['nist_hemisphere_dsc_um2']),
                      'FEM / NIST: {}'.format(metrics['fem_to_nist_hemisphere_ratio']),
                      'Angular shape cosine: {}'.format(metrics['angular_shape_cosine'])]
        if metrics.get('finite_patch_hemisphere_dsc_um2') is not None:
            lines += ['Finite patch expected DSC: {:.8g} um^2'.format(
                          metrics['finite_patch_hemisphere_dsc_um2']),
                      'FEM / finite patch: {}'.format(
                          metrics['fem_to_finite_patch_hemisphere_ratio']),
                      'FEM / finite patch shape cosine: {}'.format(
                          metrics['finite_patch_angular_shape_cosine'])]
            if metrics.get('finite_patch_to_nist_hemisphere_ratio') is not None:
                lines += ['Finite patch / infinite NIST: {:.5g}'.format(
                              metrics['finite_patch_to_nist_hemisphere_ratio']),
                          'Finite patch / NIST shape cosine: {:.5g}'.format(
                              metrics['finite_patch_to_nist_angular_shape_cosine'])]
        lines += ['Input stationary RMS: {:.6g} nm'.format(metrics['input_stationary_rms_nm']),
                  'Mean realized patch RMS: {:.6g} nm'.format(metrics['realized_patch_rms_mean_nm']),
                  'Height-grid Nyquist PSD variance fraction: {}'.format(
                      metrics.get('nyquist_variance_fraction')),
                  '', 'Limits:']+['- '+warning for warning in report['warnings']]
        if metrics.get('observables'):
            lines += ['', 'Pointwise 95% intervals (um2):']+[
                '{}: {:.6g}, CI {}'.format(name, stat['mean_um2'], stat['ci95_um2'])
                for name, stat in metrics['observables'].items()]
        convergence_path = directory/'haze_sample_convergence.png'
        if convergence_path.is_file():
            self.sample_convergence_view.load(convergence_path)
        self.ensemble_tabs.setTabEnabled(self.ensemble_tabs.indexOf(self.sample_convergence_view),
                                         convergence_path.is_file())
        self.ensemble_metrics.setPlainText('\n'.join(lines))
        self.setTabEnabled(self.indexOf(self.ensemble_tabs), True)
        self.load_illumination(directory)
        return report

    def load_convergence(self, directory):
        directory = Path(directory)
        report = json.loads((directory/'haze_convergence_report.json').read_text(encoding='utf-8'))
        if report.get('status') not in ('complete', 'complete_with_limits'):
            raise ValueError('The convergence study is incomplete; resume its saved plan first')
        self.convergence_view.load(directory/'haze_convergence.png')
        error_path = directory/'haze_paired_errors.png'
        if error_path.is_file():
            self.paired_error_view.load(error_path)
        self.convergence_tabs.setTabEnabled(self.convergence_tabs.indexOf(self.paired_error_view),
                                            error_path.is_file())
        self.convergence_details.setPlainText(json.dumps(report, indent=2))
        self.setTabEnabled(self.indexOf(self.convergence_tabs), True)
        return report
