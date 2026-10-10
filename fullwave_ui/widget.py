"""Asynchronous FEM subprocess UI with cooperative, process-tree cancellation."""
import codecs
from dataclasses import asdict, replace
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcess, QProcessEnvironment, Qt
from PyQt5.QtWidgets import (QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
    QScrollArea, QSpinBox, QSplitter, QStyle, QVBoxLayout, QWidget, QToolButton, QMenu, QComboBox)

from .controls import SimulationControls
from .results import ResultsTabs
from .runtime_paths import application_root, default_engine, results_root
from .saved_source import reusable_source
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.parameters import Parameters

ROOT = application_root()


class FullwaveWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.read_output)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)
        self.output = None
        self.loaded_result = None
        self.loaded_source = None
        self.loaded_hybrid = None
        self.loaded_whole = None
        self.imaging_only_job = False
        self.job_kind = None
        self.active_report = None
        self.cancel_file = None
        self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
        self.controls = SimulationControls()
        self.results = ResultsTabs()
        self.interpreter = QLineEdit(str(default_engine(ROOT)))
        self.browse = QPushButton(self.style().standardIcon(QStyle.SP_DirOpenIcon), '')
        self.browse.setToolTip('Select FEM Python interpreter')
        self.browse.clicked.connect(self.choose_interpreter)
        interpreter_row = QHBoxLayout()
        interpreter_row.addWidget(self.interpreter)
        interpreter_row.addWidget(self.browse)
        self.memory = QDoubleSpinBox()
        self.memory.setRange(1, 64)
        self.memory.setValue(10)
        self.memory.setSuffix(' GiB')
        self.timeout = QDoubleSpinBox()
        self.timeout.setRange(1, 240)
        self.timeout.setValue(30)
        self.timeout.setSuffix(' min')
        self.controls.preset.currentIndexChanged.connect(self.apply_runtime_preset)
        self.controls.preset_button.clicked.connect(self.apply_runtime_preset)
        self.runtime = QWidget()
        form = QFormLayout(self.runtime)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        form.addRow('FEM engine / Python', interpreter_row)
        form.addRow('Memory limit', self.memory)
        form.addRow('Time limit', self.timeout)
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.addWidget(self.runtime)
        left_layout.addWidget(self.controls)
        left_layout.addStretch()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(350)
        scroll.setWidget(left)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(3000)
        self.log.setMaximumHeight(150)
        self.status = QLabel('FEM experimental | No result loaded')
        self.status.setWordWrap(True)
        self.loaded = QLabel('')
        self.loaded.setWordWrap(True)
        self.loaded.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.start = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Calculate')
        self.cancel = QPushButton(self.style().standardIcon(QStyle.SP_MediaStop), 'Cancel')
        self.load_button = QPushButton(self.style().standardIcon(QStyle.SP_DirOpenIcon), 'Load result')
        self.reimage_button = QPushButton(self.style().standardIcon(QStyle.SP_BrowserReload), 'Re-image source')
        self.reimage_button.setToolTip('Reuse saved scattering with new imaging, camera or laser power settings. '
                                      'Changing a finite beam footprint requires Calculate.')
        self.reimage_button.clicked.connect(self.choose_source)
        self.image_current_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Image current result')
        self.image_current_button.setToolTip('Imaging, camera and laser power controls apply. '
                                            'Particle, wafer and illumination direction/shape stay from the loaded source.')
        self.image_current_button.setEnabled(False)
        self.image_current_button.clicked.connect(self.image_current_result)
        self.compare_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Compare particle / haze')
        self.compare_button.setToolTip('Compare using the current Ch1/2/3 NA and tilt settings; NA=0 disables a channel. '
                                       'Other imaging controls also apply. Particle and wafer settings come from the loaded source.')
        self.compare_button.setEnabled(False)
        self.compare_button.clicked.connect(self.start_comparison)
        self.ensemble_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Average haze')
        self.ensemble_button.setToolTip('Solve independent Gaussian or ABC rough surfaces without a particle; '
                                       'average DSC intensities and compare with finite-patch and '
                                       'infinite NIST PSD references.')
        self.ensemble_button.clicked.connect(self.start_ensemble)
        self.whole_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Whole-beam haze (SPM)')
        self.whole_button.setToolTip('Coherent first-order rough-interface scattering over the illuminated wafer; '
                                    'not a full-wafer FEM or particle-coupled solution.')
        self.whole_button.clicked.connect(self.start_whole_haze)
        self.hybrid_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Particle + whole haze (hybrid)')
        self.hybrid_button.setToolTip('Shared mother surface; whole SPM plus local FEM minus local SPM. '
                                     'Outer/local high-order coupling and outer near fields are not solved.')
        self.hybrid_button.clicked.connect(self.start_hybrid_job)
        self.hybrid_reimage = QPushButton(self.style().standardIcon(QStyle.SP_BrowserReload), 'Re-image hybrid')
        self.hybrid_reimage.setEnabled(False)
        self.hybrid_reimage.clicked.connect(self.reimage_hybrid)
        self.validation_button = QToolButton()
        self.validation_button.setText('Model checks')
        self.validation_button.setIcon(self.style().standardIcon(QStyle.SP_FileDialogDetailedView))
        self.validation_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.validation_button.setPopupMode(QToolButton.InstantPopup)
        menu = QMenu(self.validation_button)
        self.validation_actions = {}
        for mode, title in (('range', 'Whole-surface range'), ('mesh', 'Same-surface grid / taper'),
                            ('power', 'Low-height absolute power / energy balance'),
                            ('ensemble128', 'Extend whole haze to 128 surfaces'),
                            ('local', 'Low-height FEM / SPM (P1)'),
                            ('local_p2', 'Thin-limit FEM (P2, sparse)')):
            action = menu.addAction(title)
            self.validation_actions[mode] = action
            action.triggered.connect(lambda checked=False, mode=mode: self.start_model_validation(mode))
        self.validation_button.setMenu(menu)
        menu.aboutToShow.connect(self.refresh_validation_actions)
        self.validation_button.setEnabled(False)
        self.whole_coverage = QDoubleSpinBox()
        self.whole_coverage.setRange(.9, .999)
        self.whole_coverage.setDecimals(3)
        self.whole_coverage.setSingleStep(.005)
        self.whole_coverage.setValue(.99)
        self.whole_coverage.setToolTip('Fraction of whole-beam incident power covered after squared height taper; 0.99 = 99%.')
        self.whole_variance = QDoubleSpinBox()
        self.whole_variance.setRange(.5, .999)
        self.whole_variance.setDecimals(3)
        self.whole_variance.setSingleStep(.01)
        self.whole_variance.setValue(.9)
        self.whole_variance.setToolTip('Minimum parent PSD variance inside the grid Nyquist circle. '
                                      'Unmeasured high frequencies are not restored by RMS rescaling.')
        self.controls.form.addRow('Whole-beam power coverage', self.whole_coverage)
        self.controls.form.addRow('Whole-beam PSD variance fraction', self.whole_variance)
        self.mother_factor = QComboBox()
        self.mother_factor.addItems(('2', '4', '8'))
        self.mother_factor.setCurrentText('8')
        self.mother_factor.setToolTip('Mother width / initial rough width. A larger mother needs fresh local FEM; existing sources are not extended.')
        self.controls.form.addRow('Hybrid mother range factor', self.mother_factor)
        self.convergence_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Convergence')
        self.convergence_button.clicked.connect(self.start_convergence)
        self.ensemble_count = QSpinBox()
        self.ensemble_count.setRange(1, 128)
        self.ensemble_count.setValue(16)
        self.ensemble_count.setToolTip('Independent rough-surface samples. Statistical convergence needs '
                                      'more samples and a patch-width study.')
        self.ensemble_taper = QDoubleSpinBox()
        self.ensemble_taper.setRange(.001, 1.)
        self.ensemble_taper.setDecimals(3)
        self.ensemble_taper.setSingleStep(.005)
        self.ensemble_taper.setValue(.025)
        self.ensemble_taper.setSuffix(' um')
        self.ensemble_taper.setToolTip('Physical width tapered at each patch edge. '
                                      'Must cover at least one surface-grid interval.')
        self.ensemble_area = QDoubleSpinBox()
        self.ensemble_area.setRange(.0001, 10000)
        self.ensemble_area.setDecimals(4)
        self.ensemble_area.setValue(.04)
        self.ensemble_area.setSuffix(' um2')
        self.ensemble_area.setToolTip('Common area for FEM and NIST DSC comparison; '
                                     '0.04 um2 reproduces the earlier 0.2 x 0.2 um plot.')
        self.fit_button = QPushButton(self.style().standardIcon(QStyle.SP_TitleBarMaxButton), '')
        self.fit_button.setToolTip('Fit figure to view')
        self.fit_button.clicked.connect(self.results.fit_current)
        self.start.clicked.connect(self.start_job)
        self.cancel.clicked.connect(self.cancel_job)
        self.load_button.clicked.connect(self.choose_result)
        self.cancel.setEnabled(False)
        toolbar = QHBoxLayout()
        for button in (self.start, self.cancel, self.load_button, self.fit_button):
            toolbar.addWidget(button)
        toolbar.addStretch()
        imaging_toolbar = QHBoxLayout()
        imaging_toolbar.addWidget(self.image_current_button)
        imaging_toolbar.addWidget(self.reimage_button)
        imaging_toolbar.addWidget(self.compare_button)
        imaging_toolbar.addStretch()
        ensemble_toolbar = QHBoxLayout()
        ensemble_toolbar.addWidget(self.ensemble_button)
        ensemble_toolbar.addWidget(self.convergence_button)
        for label, field in (('Surfaces', self.ensemble_count),
                             ('Edge taper', self.ensemble_taper),
                             ('Compare area', self.ensemble_area)):
            ensemble_toolbar.addWidget(QLabel(label))
            ensemble_toolbar.addWidget(field)
        ensemble_toolbar.addStretch()
        whole_toolbar = QHBoxLayout()
        whole_toolbar.addWidget(self.whole_button)
        whole_toolbar.addWidget(self.hybrid_button)
        whole_toolbar.addWidget(self.hybrid_reimage)
        whole_toolbar.addWidget(self.validation_button)
        whole_toolbar.addStretch()
        right = QWidget()
        layout = QVBoxLayout(right)
        layout.addLayout(toolbar)
        layout.addLayout(imaging_toolbar)
        layout.addLayout(ensemble_toolbar)
        layout.addLayout(whole_toolbar)
        layout.addWidget(self.status)
        layout.addWidget(self.loaded)
        layout.addWidget(self.results, 1)
        layout.addWidget(self.progress)
        layout.addWidget(self.log)
        splitter = QSplitter()
        splitter.addWidget(scroll)
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([385, 1000])
        outer = QVBoxLayout(self)
        outer.addWidget(splitter)

    def choose_interpreter(self):
        path, _ = QFileDialog.getOpenFileName(self, 'FEM engine / Python', self.interpreter.text(), 'Executable (*.exe)')
        if path:
            self.interpreter.setText(path)

    def apply_runtime_preset(self, *_):
        if self.controls.preset.currentIndex() == 2:
            return
        self.memory.setValue(10)
        self.timeout.setValue(30)

    def set_running(self, running):
        self.controls.setEnabled(not running)
        self.runtime.setEnabled(not running)
        self.start.setEnabled(not running)
        self.load_button.setEnabled(not running)
        self.reimage_button.setEnabled(not running)
        self.image_current_button.setEnabled(not running and self.loaded_source is not None)
        self.compare_button.setEnabled(not running and self.loaded_source is not None and
                                       self.loaded_result is not None)
        for field in (self.ensemble_button, self.convergence_button, self.ensemble_count,
                      self.ensemble_taper, self.ensemble_area, self.whole_button, self.hybrid_button):
            field.setEnabled(not running)
        self.hybrid_reimage.setEnabled(not running and self.loaded_hybrid is not None)
        self.validation_button.setEnabled(not running and (self.loaded_hybrid is not None or self.loaded_whole is not None))
        self.cancel.setEnabled(running)
        self.progress.setRange(0, 0 if running else 1)
        self.progress.setValue(0)

    def choose_source(self):
        directory = QFileDialog.getExistingDirectory(self, 'Reuse scattering folder (imaging settings only)', str(results_root(ROOT)))
        if directory:
            source = Path(directory)
            if (source/'scattering').is_dir():
                source = source/'scattering'
            self.start_job(source=source)

    def image_current_result(self):
        if self.loaded_source is not None:
            self.start_job(source=self.loaded_source)

    def start_job(self, checked=False, source=None):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            source_report = None
            if source is not None:
                source_report = json.loads((Path(source)/'report.json').read_text(encoding='utf-8'))
                if source_report.get('source_format_version') != 1 or not (Path(source)/'paired_fields.npz').is_file():
                    raise ValueError('Select a complete reusable scattering source folder.')
            source_incidence = (source_report.get('parameters', {}).get('incidence_deg')
                                if source_report else None)
            if source_report and source_incidence is None:
                raise ValueError('Reusable scattering report has no illumination incidence angle.')
            arguments = self.controls.arguments(imaging_only=source is not None,
                                                incidence_deg=source_incidence,
                                                parameters=Parameters(**source_report['parameters']) if source_report else None)
            arguments += ['--pupil-cache', str(results_root(ROOT)/'pupil_cache')]
            if source is not None:
                arguments += ['--source', str(Path(source).resolve())]
            interpreter = Path(self.interpreter.text())
            if not interpreter.is_file():
                raise ValueError('FEM Python is missing. Select the isolated NGSolve environment.')
            self.output = results_root(ROOT)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            self.output.mkdir(parents=True, exist_ok=False)
            self.cancel_file = self.output/'cancel.request'
            self.active_report = self.output/'resources.json'
            self.job_kind = 'imaging' if source is not None else 'simulation'
            command = ['-m', 'fullwave_prototype.guarded_run', '--report',
                str(self.active_report), '--max-memory-gib', str(self.memory.value()),
                '--timeout-seconds', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file),
                'imaging_run', '--output', str(self.output)]+arguments
            (self.output/'request.json').write_text(json.dumps(dict(interpreter=str(interpreter),
                arguments=command), indent=2), encoding='utf-8')
            self.process.setWorkingDirectory(str(ROOT))
            env = QProcessEnvironment.systemEnvironment()
            env.insert('PYTHONIOENCODING', 'utf-8')
            env.insert('PYTHONUNBUFFERED', '1')
            env.insert('PYINSTALLER_RESET_ENVIRONMENT', '1')
            self.process.setProcessEnvironment(env)
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.imaging_only_job = source is not None
            self.log.clear()
            self.set_running(True)
            self.status.setText('Calculating new result | '+str(self.output))
            if source is not None:
                self.status.setText('Re-imaging saved scattering (physical controls not applied) | '+str(source))
            self.process.start(str(interpreter), command)
        except (ValueError, OSError, KeyError) as exc:
            QMessageBox.warning(self, 'Cannot start FEM', str(exc))

    def start_comparison(self):
        if self.process.state() != QProcess.NotRunning or self.loaded_result is None:
            return
        try:
            if self.loaded_source is None:
                raise ValueError('Load a complete reusable FEM result before comparing components.')
            interpreter = Path(self.interpreter.text())
            if not interpreter.is_file():
                raise ValueError('FEM Python is missing. Select the isolated NGSolve environment.')
            source_report = json.loads((self.loaded_source/'report.json').read_text(encoding='utf-8'))
            source_incidence = source_report.get('parameters', {}).get('incidence_deg')
            if source_incidence is None:
                raise ValueError('Reusable scattering report has no illumination incidence angle.')
            imaging_args = build_parser().parse_args(self.controls.arguments(
                imaging_only=True, incidence_deg=source_incidence,
                parameters=Parameters(**source_report['parameters'])))
            settings = imaging_settings(imaging_args, incidence_deg=source_incidence,
                                        parameters=Parameters(**source_report['parameters']))
            self.output = self.loaded_result
            stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            settings_path = self.output/('comparison_imaging_settings_'+stamp+'.json')
            settings_path.write_text(json.dumps(asdict(settings), indent=2), encoding='utf-8')
            self.cancel_file = self.output/('comparison_cancel_'+stamp+'.request')
            self.active_report = self.output/('comparison_resources_'+stamp+'.json')
            self.job_kind = 'comparison'
            command = ['-m', 'fullwave_prototype.guarded_run', '--report',
                str(self.active_report), '--max-memory-gib', str(self.memory.value()),
                '--timeout-seconds', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file),
                'component_compare', '--result', str(self.output), '--pupil-cache',
                str(results_root(ROOT)/'pupil_cache'), '--imaging-settings-json', str(settings_path)]
            (self.output/('comparison_request_'+stamp+'.json')).write_text(json.dumps(dict(
                interpreter=str(interpreter), arguments=command), indent=2), encoding='utf-8')
            self.process.setWorkingDirectory(str(ROOT))
            env = QProcessEnvironment.systemEnvironment()
            env.insert('PYTHONIOENCODING', 'utf-8')
            env.insert('PYTHONUNBUFFERED', '1')
            env.insert('PYINSTALLER_RESET_ENVIRONMENT', '1')
            self.process.setProcessEnvironment(env)
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.log.clear()
            self.set_running(True)
            self.status.setText('Computing particle-only reference | '+str(self.output))
            self.process.start(str(interpreter), command)
        except (ValueError, OSError) as exc:
            QMessageBox.warning(self, 'Cannot compare components', str(exc))

    def start_ensemble(self):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            from fullwave_prototype.far_field import angular_grid
            from fullwave_prototype.haze_crosscheck import statistical_reference, finite_patch_reference
            import numpy as np

            arguments = self.controls.arguments(particle_enabled=False)
            args = build_parser().parse_args(arguments)
            parameters, roughness = scattering_settings(args, particle_enabled=False)
            roughness = replace(roughness, taper_fraction=(self.ensemble_taper.value()/
                                                          roughness.width_um))
            roughness.validate(parameters, particle_enabled=False)
            if roughness.model not in ('gaussian', 'abc', 'measured_psd') or roughness.rms_nm <= 0 or parameters.beam_waist_um != 0:
                raise ValueError('Use Gaussian/ABC/measured PSD and plane/Gaussian footprint mode for haze averaging. '
                                 'Measured height is a single deterministic surface.')
            if roughness.seed+self.ensemble_count.value()-1 > 1000000:
                raise ValueError('Last roughness seed exceeds 1000000.')
            interpreter = Path(self.interpreter.text())
            if not interpreter.is_file():
                raise ValueError('FEM engine is missing. Select the isolated NGSolve environment.')
            theta, phi, _ = angular_grid(5)
            reference = finite_reference = None
            if parameters.beam_profile == 'plane':
                reference = statistical_reference(theta, phi, asdict(parameters), asdict(roughness),
                                               self.ensemble_area.value())
                finite_reference = finite_patch_reference(
                    theta, phi, asdict(parameters), asdict(roughness), self.ensemble_area.value())
            self.output = results_root(ROOT)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            self.output.mkdir(parents=True, exist_ok=False)
            reference_path = self.output/'nist_psd_reference.npz' if reference is not None else None
            measurement_sha = ''
            if roughness.model == 'measured_psd':
                from fullwave_prototype.experimental_surface import input_fingerprint
                measurement_sha = input_fingerprint(roughness.measurement_path)
            if reference_path is not None:
                np.savez_compressed(reference_path, theta=theta, phi=phi, dsc=reference,
                                finite_patch_dsc=finite_reference,
                                psd_model=roughness.model,
                                parameters_json=json.dumps(asdict(parameters)),
                                roughness_json=json.dumps(asdict(roughness)),
                                measurement_sha256=measurement_sha,
                                comparison_area_um2=self.ensemble_area.value())
            self.cancel_file = self.output/'cancel.request'
            self.active_report = self.output/'resources.json'
            self.job_kind = 'ensemble'
            command = ['-m', 'fullwave_prototype.guarded_run', '--report',
                str(self.active_report), '--max-memory-gib', str(self.memory.value()),
                '--timeout-seconds', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file),
                'haze_ensemble', '--output', str(self.output), '--samples',
                str(self.ensemble_count.value()), '--taper-fraction', str(roughness.taper_fraction),
                '--comparison-area-um2', str(self.ensemble_area.value()),
                ]+(['--reference', str(reference_path)] if reference_path is not None else [])+arguments
            (self.output/'request.json').write_text(json.dumps(dict(interpreter=str(interpreter),
                arguments=command), indent=2), encoding='utf-8')
            self.process.setWorkingDirectory(str(ROOT))
            env = QProcessEnvironment.systemEnvironment()
            env.insert('PYTHONIOENCODING', 'utf-8')
            env.insert('PYTHONUNBUFFERED', '1')
            env.insert('PYINSTALLER_RESET_ENVIRONMENT', '1')
            self.process.setProcessEnvironment(env)
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.log.clear()
            self.set_running(True)
            self.status.setText('Averaging {} independent FEM haze surfaces'.format(
                self.ensemble_count.value()))
            self.process.start(str(interpreter), command)
        except (ValueError, OSError, RuntimeError, KeyError) as exc:
            QMessageBox.warning(self, 'Cannot average haze', str(exc))

    def start_convergence(self):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            from .convergence_dialog import ConvergenceDialog
            from fullwave_prototype.haze_convergence import default_plan
            args = build_parser().parse_args(self.controls.arguments(particle_enabled=False))
            p, rough = scattering_settings(args, particle_enabled=False)
            rough = replace(rough, taper_fraction=self.ensemble_taper.value()/rough.width_um)
            samples = 1 if rough.model == 'measured_height' else min(8, self.ensemble_count.value())
            plan = default_plan(p, rough, samples, args.na, args.tilts)
            plan['comparison_area_um2'] = self.ensemble_area.value()
            dialog = ConvergenceDialog(plan, self)
            if not dialog.exec_():
                return
            interpreter = Path(self.interpreter.text())
            if not interpreter.is_file():
                raise ValueError('Select the isolated FEM Python environment')
            self.output = results_root(ROOT)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            self.output.mkdir(parents=True, exist_ok=False)
            plan_path = self.output/'convergence_plan.json'
            plan_path.write_text(json.dumps(dialog.plan, indent=2), encoding='utf-8')
            self.cancel_file = self.output/'cancel.request'
            self.active_report = self.output/'resources.json'
            self.job_kind = 'convergence'
            command = ['-m', 'fullwave_prototype.guarded_run', '--report', str(self.active_report),
                '--max-memory-gib', str(self.memory.value()), '--timeout-seconds', str(self.timeout.value()*60),
                '--cancel-file', str(self.cancel_file), 'haze_convergence', '--input', str(plan_path),
                '--output', str(self.output), '--max-memory-gib', str(self.memory.value()),
                '--sample-timeout-seconds', str(self.timeout.value()*60)]
            (self.output/'request.json').write_text(json.dumps(dict(interpreter=str(interpreter),
                arguments=command), indent=2), encoding='utf-8')
            self.process.setWorkingDirectory(str(ROOT))
            env = QProcessEnvironment.systemEnvironment()
            for key, value in (('PYTHONIOENCODING', 'utf-8'), ('PYTHONUNBUFFERED', '1'),
                               ('PYINSTALLER_RESET_ENVIRONMENT', '1')):
                env.insert(key, value)
            self.process.setProcessEnvironment(env)
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.log.clear()
            self.set_running(True)
            self.status.setText('Scanning haze convergence | '+str(self.output))
            self.process.start(str(interpreter), command)
        except (ValueError, OSError, KeyError) as exc:
            QMessageBox.warning(self, 'Cannot start convergence', str(exc))

    def start_whole_haze(self):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            from .whole_haze_job import start_whole_haze
            start_whole_haze(self)
        except (ValueError, OSError, RuntimeError, KeyError) as exc:
            QMessageBox.warning(self, 'Cannot calculate whole-beam haze', str(exc))

    def start_hybrid_job(self, checked=False, source=None):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            from .hybrid_job import start_hybrid
            start_hybrid(self, source)
        except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot calculate hybrid fields', str(exc))

    def reimage_hybrid(self):
        if self.loaded_hybrid is not None:
            self.start_hybrid_job(source=self.loaded_hybrid)

    def start_model_validation(self, mode):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            from .validation_job import start_validation
            start_validation(self, mode)
        except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot validate model', str(exc))

    def refresh_validation_actions(self):
        for mode, action in self.validation_actions.items():
            action.setEnabled(self.loaded_whole is not None if mode == 'ensemble128' else self.loaded_hybrid is not None)

    def load_model_validation(self, directory):
        try:
            report = self.results.load_validation(directory)
            self.results.setCurrentWidget(self.results.validation_tabs)
            conclusion = ('1% adjacent range target '+
                ('met' if report['last_adjacent_complex_target_met'] else 'NOT met')
                if 'last_adjacent_complex_target_met' in report else 'Selected-direction low-height check')
            if 'tested_spatial_taper_target_met' in report:
                conclusion = '1% grid/taper target '+('met' if report['tested_spatial_taper_target_met'] else 'NOT met')
            self.status.setText(conclusion+' | Not full-model certification')
            self.loaded.setText('Validation: '+str(Path(directory).resolve()))
            return True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load validation', str(exc))
            return False

    def read_output(self):
        text = self.decoder.decode(bytes(self.process.readAllStandardOutput()))
        self.log.moveCursor(self.log.textCursor().End)
        self.log.insertPlainText(text)
        self.log.ensureCursorVisible()

    def cancel_job(self):
        if self.process.state() != QProcess.NotRunning and self.cancel_file:
            try:
                self.cancel_file.touch()
                self.cancel.setEnabled(False)
                self.status.setText('Cancelling | Waiting for owned processes to exit')
            except OSError as exc:
                QMessageBox.warning(self, 'Cannot request cancellation', str(exc))

    def process_error(self, error):
        if error == QProcess.FailedToStart:
            self.set_running(False)
            self.status.setText('FEM interpreter failed to start: '+self.process.errorString())

    def finished(self, exit_code, exit_status):
        self.read_output()
        self.set_running(False)
        status = 'failed'
        report = self.active_report
        if report.exists():
            try:
                status = json.loads(report.read_text(encoding='utf-8'))['status']
            except (ValueError, KeyError, OSError):
                pass
        if exit_code == 0 and exit_status == QProcess.NormalExit and status == 'completed':
            if self.job_kind == 'model_validation':
                if self.load_model_validation(self.output):
                    self.progress.setValue(1)
            elif self.job_kind == 'hybrid':
                if self.load_hybrid_result(self.output):
                    self.progress.setValue(1)
            elif self.job_kind == 'whole_haze':
                self.load_whole_haze_result(self.output)
                self.progress.setValue(1)
            elif self.job_kind == 'convergence':
                self.results.load_convergence(self.output)
                self.results.setCurrentWidget(self.results.convergence_tabs)
                self.status.setText('Convergence study loaded')
                self.loaded.setText(str(self.output))
                self.progress.setValue(1)
            elif self.job_kind == 'ensemble':
                if self.load_ensemble_result(self.output):
                    self.results.setCurrentWidget(self.results.ensemble_tabs)
                    self.progress.setValue(1)
            elif self.load_result(self.output, restore_imaging_controls=self.job_kind != 'comparison'):
                if self.job_kind == 'comparison':
                    self.results.setCurrentWidget(self.results.comparison_tabs)
                    self.status.setText('Particle / haze comparison loaded')
                elif self.imaging_only_job:
                    self.results.setCurrentWidget(self.results.image_tabs)
                self.progress.setValue(1)
        else:
            reasons = {'stopped_memory_limit': 'Stopped: memory limit reached',
                       'stopped_system_memory_reserve': 'Stopped: low available system memory',
                       'stopped_timeout': 'Stopped: time limit reached',
                       'stopped_cancelled': 'Cancelled', 'failed': 'Calculation failed'}
            self.status.setText(reasons.get(status, status)+' | No new result | '+str(report))

    def choose_result(self):
        directory = QFileDialog.getExistingDirectory(self, 'Load FEM result', str(results_root(ROOT)))
        if directory:
            if (Path(directory)/'range_report.json').is_file() or (Path(directory)/'validation_report.json').is_file():
                self.load_model_validation(directory)
            elif (Path(directory)/'surface_mesh_report.json').is_file():
                self.load_model_validation(directory)
            elif (Path(directory)/'hybrid_report.json').is_file():
                self.load_hybrid_result(directory)
            elif (Path(directory)/'whole_haze_report.json').is_file():
                self.load_whole_haze_result(directory)
            elif (Path(directory)/'haze_convergence_report.json').is_file():
                try:
                    self.results.load_convergence(directory)
                    self.results.setCurrentWidget(self.results.convergence_tabs)
                    self.loaded.setText(str(Path(directory).resolve()))
                except (OSError, ValueError, KeyError) as exc:
                    QMessageBox.warning(self, 'Cannot load convergence', str(exc))
            elif (Path(directory)/'haze_ensemble_report.json').is_file():
                self.load_ensemble_result(directory)
            else:
                self.load_result(directory)

    def load_ensemble_result(self, directory):
        try:
            report = self.results.load_ensemble(directory)
            metrics = report['metrics']
            self.loaded.setText('Haze ensemble: {} surfaces | FEM {:.4g} um2 | '
                'NIST {:.4g} um2 | patch RMS {:.4g} nm'.format(
                    metrics['samples'], metrics['fem_mean_hemisphere_dsc_um2'],
                    metrics.get('nist_hemisphere_dsc_um2', float('nan')),
                    metrics['realized_patch_rms_mean_nm']))
            self.status.setText('Haze ensemble loaded | Physical convergence unverified')
            return True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load haze ensemble', str(exc))
            return False

    def load_whole_haze_result(self, directory):
        try:
            report = self.results.load_whole_haze(directory)
            self.controls.load_beam_parameters(report['parameters'])
            self.controls.load_imaging_settings(report['imaging_settings'], report['parameters']['incidence_deg'])
            self.loaded_whole = Path(directory).resolve()
            self.loaded_hybrid = self.loaded_source = self.loaded_result = None
            for button in (self.hybrid_reimage, self.image_current_button, self.compare_button):
                button.setEnabled(False)
            self.validation_button.setEnabled(self.process.state() == QProcess.NotRunning)
            self.results.setCurrentWidget(self.results.whole_haze_tabs)
            self.loaded.setText('Whole-beam SPM haze | {} surfaces | coverage {:.3%} | hemisphere {:.4g} W'.format(
                report['metrics']['samples'], report['plan']['checked_taper_weighted_power_coverage'],
                report['metrics']['mean_power_w']['hemisphere']))
            target = report.get('statistics', {}).get('statistical_target_met', False)
            self.status.setText('Whole-beam first-order surface result | NOT full-wafer FEM | statistical target '+
                                ('met' if target else 'NOT met'))
            return True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load whole-beam haze', str(exc))
            return False

    def load_hybrid_result(self, directory):
        try:
            report = json.loads((Path(directory)/'hybrid_report.json').read_text(encoding='utf-8'))
            if report.get('model') != 'whole_spm_local_fem_hybrid':
                raise ValueError('Select a phase-aligned hybrid result')
            powers = self.results.load(directory)
            self.loaded_whole = None
            self.controls.load_beam_parameters(report['scattering_parameters'])
            self.controls.load_imaging_settings(report['settings'],
                                                incidence_deg=report['scattering_parameters']['incidence_deg'])
            source = (Path(directory) if (Path(directory)/'hybrid_source_report.json').is_file()
                      else Path(report['source_directory']))
            self.loaded_hybrid = source.resolve() if (source/'hybrid_fields.npz').is_file() else None
            self.loaded_source = None
            self.loaded_result = None
            self.image_current_button.setEnabled(False)
            self.compare_button.setEnabled(False)
            self.hybrid_reimage.setEnabled(self.loaded_hybrid is not None and self.process.state() == QProcess.NotRunning)
            self.validation_button.setEnabled(self.loaded_hybrid is not None and self.process.state() == QProcess.NotRunning)
            self.loaded.setText('Hybrid: '+str(Path(directory).resolve())+'\nImage powers (pW): '+
                                ' / '.join('{:.4g}'.format(v) for v in powers))
            audit = report.get('surface_convergence', {})
            self.status.setText('Shared whole SPM + local FEM | Near fields LOCAL | outer numerical target '+
                ('met' if audit.get('outer_spm_numerical_target_met') else 'NOT met / unchecked')+
                ' | Local FEM/SPM complex difference {:.2%}'.format(report['local_fem_vs_spm_complex_relative_l2']))
            return True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load hybrid result', str(exc))
            return False

    def load_result(self, directory, restore_imaging_controls=True):
        if (Path(directory)/'hybrid_report.json').is_file():
            return self.load_hybrid_result(directory)
        try:
            powers = self.results.load(directory)
            self.loaded_hybrid = None
            self.loaded_whole = None
            self.hybrid_reimage.setEnabled(False)
            self.validation_button.setEnabled(False)
            if restore_imaging_controls:
                report = json.loads((Path(directory)/'imaging_report.json').read_text(encoding='utf-8'))
                if 'settings' in report:
                    self.controls.load_beam_parameters(report.get('scattering_parameters', {}))
                    incidence = report.get('scattering_parameters', {}).get('incidence_deg')
                    self.controls.load_imaging_settings(report['settings'], incidence_deg=incidence)
            self.loaded_result = Path(directory).resolve()
            self.loaded_source = reusable_source(directory)
            roughness_note = ''
            if self.loaded_source is not None:
                source_report = json.loads((self.loaded_source/'report.json').read_text(encoding='utf-8'))
                roughness = source_report.get('roughness')
                if roughness and all(key in roughness for key in ('model', 'grid_points', 'rms_nm')):
                    surface = source_report.get('surface_statistics', {})
                    roughness_note = '\nFEM roughness: {} | {} nodes/axis | input RMS {:.6g} nm | patch RMS {:.6g} nm'.format(
                        roughness['model'], roughness['grid_points'], roughness['rms_nm'],
                        surface.get('triangle_rms_nm', roughness['rms_nm']))
                    if roughness['model'] == 'legacy':
                        roughness_note += '\nLegacy surface has no equivalent Gaussian PSD length L.'
            self.image_current_button.setEnabled(self.loaded_source is not None and
                                                self.process.state() == QProcess.NotRunning)
            self.compare_button.setEnabled(self.loaded_source is not None and
                                           self.process.state() == QProcess.NotRunning)
            self.loaded.setText('Displayed result: '+str(Path(directory).resolve())+
                '\nImage powers (pW), channels 1 / 2 / 3: '+
                ' / '.join('{:.4g}'.format(power) for power in powers)+roughness_note)
            self.status.setText('Experimental FEM result loaded')
            self.progress.setValue(1)
            return True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load result', str(exc))
            return False

    def wait_for_shutdown(self):
        if self.process.state() == QProcess.NotRunning:
            return True
        self.cancel_job()
        return self.process.waitForFinished(25000)
