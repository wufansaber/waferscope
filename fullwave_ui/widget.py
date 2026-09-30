"""Asynchronous FEM subprocess UI with cooperative, process-tree cancellation."""
import codecs
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcess, QProcessEnvironment, Qt
from PyQt5.QtWidgets import (QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
    QScrollArea, QSplitter, QStyle, QVBoxLayout, QWidget)

from .controls import SimulationControls
from .results import ResultsTabs
from .runtime_paths import application_root, default_engine, results_root
from .saved_source import reusable_source
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.imaging_run import build_parser

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
        self.memory.setValue(12)
        self.memory.setSuffix(' GiB')
        self.timeout = QDoubleSpinBox()
        self.timeout.setRange(1, 240)
        self.timeout.setValue(10)
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
        self.reimage_button.setToolTip('Reuse saved scattering; imaging, camera and illumination power/area controls apply.')
        self.reimage_button.clicked.connect(self.choose_source)
        self.image_current_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Image current result')
        self.image_current_button.setToolTip('Imaging, camera and illumination power/area controls apply. '
                                            'Particle, wafer and illumination direction/shape stay from the loaded source.')
        self.image_current_button.setEnabled(False)
        self.image_current_button.clicked.connect(self.image_current_result)
        self.compare_button = QPushButton(self.style().standardIcon(QStyle.SP_MediaPlay), 'Compare particle / haze')
        self.compare_button.setToolTip('Compare using the current Ch1/2/3 NA and tilt settings; NA=0 disables a channel. '
                                       'Other imaging controls also apply. Particle and wafer settings come from the loaded source.')
        self.compare_button.setEnabled(False)
        self.compare_button.clicked.connect(self.start_comparison)
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
        right = QWidget()
        layout = QVBoxLayout(right)
        layout.addLayout(toolbar)
        layout.addLayout(imaging_toolbar)
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
        self.memory.setValue(12)
        self.timeout.setValue(10 if self.controls.preset.currentIndex() == 0 else 30)

    def set_running(self, running):
        self.controls.setEnabled(not running)
        self.runtime.setEnabled(not running)
        self.start.setEnabled(not running)
        self.load_button.setEnabled(not running)
        self.reimage_button.setEnabled(not running)
        self.image_current_button.setEnabled(not running and self.loaded_source is not None)
        self.compare_button.setEnabled(not running and self.loaded_source is not None and
                                       self.loaded_result is not None)
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
                                                incidence_deg=source_incidence)
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
        except (ValueError, OSError) as exc:
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
                imaging_only=True, incidence_deg=source_incidence))
            settings = imaging_settings(imaging_args, incidence_deg=source_incidence)
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
            if self.load_result(self.output, restore_imaging_controls=self.job_kind != 'comparison'):
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
            self.load_result(directory)

    def load_result(self, directory, restore_imaging_controls=True):
        try:
            powers = self.results.load(directory)
            if restore_imaging_controls:
                report = json.loads((Path(directory)/'imaging_report.json').read_text(encoding='utf-8'))
                if 'settings' in report:
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
