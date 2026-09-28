"""Asynchronous FEM subprocess UI with cooperative, process-tree cancellation."""
import codecs
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
        self.timeout.setValue(30)
        self.timeout.setSuffix(' min')
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
        self.reimage_button.setToolTip('Reuse saved scattering; only imaging and camera controls apply.')
        self.reimage_button.clicked.connect(self.choose_source)
        self.fit_button = QPushButton(self.style().standardIcon(QStyle.SP_TitleBarMaxButton), '')
        self.fit_button.setToolTip('Fit figure to view')
        self.fit_button.clicked.connect(self.results.fit_current)
        self.start.clicked.connect(self.start_job)
        self.cancel.clicked.connect(self.cancel_job)
        self.load_button.clicked.connect(self.choose_result)
        self.cancel.setEnabled(False)
        toolbar = QHBoxLayout()
        for button in (self.start, self.cancel, self.load_button, self.reimage_button, self.fit_button):
            toolbar.addWidget(button)
        toolbar.addStretch()
        right = QWidget()
        layout = QVBoxLayout(right)
        layout.addLayout(toolbar)
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

    def set_running(self, running):
        self.controls.setEnabled(not running)
        self.runtime.setEnabled(not running)
        self.start.setEnabled(not running)
        self.load_button.setEnabled(not running)
        self.reimage_button.setEnabled(not running)
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

    def start_job(self, checked=False, source=None):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            arguments = self.controls.arguments(imaging_only=source is not None)
            arguments += ['--pupil-cache', str(results_root(ROOT)/'pupil_cache')]
            if source is not None:
                report = json.loads((Path(source)/'report.json').read_text(encoding='utf-8'))
                if report.get('source_format_version') != 1 or not (Path(source)/'paired_fields.npz').is_file():
                    raise ValueError('Select a complete reusable scattering source folder.')
                arguments += ['--source', str(Path(source).resolve())]
            interpreter = Path(self.interpreter.text())
            if not interpreter.is_file():
                raise ValueError('FEM Python is missing. Select the isolated NGSolve environment.')
            self.output = results_root(ROOT)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            self.output.mkdir(parents=True, exist_ok=False)
            self.cancel_file = self.output/'cancel.request'
            command = ['-m', 'fullwave_prototype.guarded_run', '--report',
                str(self.output/'resources.json'), '--max-memory-gib', str(self.memory.value()),
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
            self.log.clear()
            self.set_running(True)
            self.status.setText('Calculating new result | '+str(self.output))
            if source is not None:
                self.status.setText('Re-imaging saved scattering (physical controls not applied) | '+str(source))
            self.process.start(str(interpreter), command)
        except (ValueError, OSError) as exc:
            QMessageBox.warning(self, 'Cannot start FEM', str(exc))

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
        report = self.output/'resources.json'
        if report.exists():
            try:
                status = json.loads(report.read_text(encoding='utf-8'))['status']
            except (ValueError, KeyError, OSError):
                pass
        if exit_code == 0 and exit_status == QProcess.NormalExit and status == 'completed':
            self.load_result(self.output)
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

    def load_result(self, directory):
        try:
            powers = self.results.load(directory)
            self.loaded.setText('Displayed result: '+str(Path(directory).resolve())+
                '\nImage powers (pW), channels 1 / 2 / 3: '+
                ' / '.join('{:.4g}'.format(power) for power in powers))
            self.status.setText('Experimental FEM result loaded')
            self.progress.setValue(1)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, 'Cannot load result', str(exc))

    def wait_for_shutdown(self):
        if self.process.state() == QProcess.NotRunning:
            return True
        self.cancel_job()
        return self.process.waitForFinished(25000)
