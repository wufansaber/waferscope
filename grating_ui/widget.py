"""Native S4 runs only in an independently guarded subprocess, never in Qt."""
import codecs
from datetime import datetime
import json
from pathlib import Path
import shutil

from PyQt5.QtCore import QProcess, QProcessEnvironment, Qt
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QProgressBar,
    QPushButton, QScrollArea, QSplitter, QStyle, QVBoxLayout, QWidget)

from fullwave_ui.runtime_paths import application_root, results_root
from grating_scene.pipeline import is_within, load_scene
from grating_scene.runtime import default_engine, default_s4
from grating_scene.examples import EXAMPLES, example_folder, example_inputs
from .controls import SceneControls
from .camera_controls import PeriodicCameraControls
from .results import SceneResults


class GratingWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.root = application_root()
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.read_output)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)
        self.loaded_result = self.output = self.cancel_file = None
        self.camera_result = None
        self.continuous_result = None
        self.job_kind = 'scene'
        self.controls, self.results = SceneControls(), SceneResults()
        self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
        self.pending_stage = ''
        self.engine = QLineEdit(str(default_engine(self.root)))
        self.executable = QLineEdit(str(default_s4(self.root)))
        self.memory_profile = QLineEdit()
        self.memory_profile.setToolTip('Optional JSON calibration for this exact S4 executable; guards stay unchanged')
        self.memory, self.timeout = QDoubleSpinBox(), QDoubleSpinBox()
        self.memory.setRange(1, 64)
        self.memory.setValue(8)
        self.memory.setSuffix(' GiB')
        self.timeout.setRange(1, 240)
        self.timeout.setValue(60)
        self.timeout.setSuffix(' min')
        left, runtime = QWidget(), QWidget()
        form = QFormLayout(runtime)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        for label, control in (('Worker Python / EXE', self.engine), ('S4 stripes executable', self.executable),
            ('Memory profile (optional)', self.memory_profile), ('Process-tree memory cap', self.memory), ('Time cap', self.timeout)):
            form.addRow(label, control)
        sidebar = QVBoxLayout(left)
        sidebar.addWidget(runtime)
        sidebar.addWidget(self.controls)
        self.camera_controls = PeriodicCameraControls()
        sidebar.addWidget(self.camera_controls)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(365)
        scroll.setWidget(left)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2500)
        self.log.setMaximumHeight(125)
        self.status = QLabel('S4 periodic model | No result loaded | isolation unchecked')
        self.status.setWordWrap(True)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.start = self.button(QStyle.SP_MediaPlay, 'Calculate', self.start_job)
        self.reimage = self.button(QStyle.SP_BrowserReload, 'Update optics', lambda: self.start_job(reuse=True))
        self.load_button = self.button(QStyle.SP_DirOpenIcon, 'Load', self.choose_result)
        self.export = self.button(QStyle.SP_DialogSaveButton, 'Export raw bundle', self.export_result)
        self.cancel = self.button(QStyle.SP_MediaStop, 'Cancel', self.cancel_job)
        self.check_button = self.button(QStyle.SP_FileDialogDetailedView, '', self.start_validation)
        self.check_button.setToolTip('Check current scene cutoff and slicing convergence')
        self.check_button.setEnabled(False)
        self.camera_button = self.button(QStyle.SP_ComputerIcon, '', self.start_camera)
        self.camera_button.setToolTip('Simulate area camera from saved periodic optical fields')
        self.camera_button.setEnabled(False)
        self.continuous_button = self.button(QStyle.SP_FileDialogDetailedView, 'Continuous angles', self.start_continuous)
        self.reimage.setEnabled(False)
        self.export.setEnabled(False)
        self.cancel.setEnabled(False)
        self.log_scale = QCheckBox('Log color scale')
        self.log_scale.setChecked(True)
        self.log_scale.toggled.connect(self.results.set_logarithmic)
        examples = QHBoxLayout()
        self.example_choice = QComboBox()
        self.example_choice.setMinimumContentsLength(20)
        self.example_choice.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        for row in EXAMPLES:
            self.example_choice.addItem(row['title'], row['key'])
        self.apply_example_button = self.button(QStyle.SP_DialogApplyButton, 'Apply parameters', self.apply_example)
        self.open_example_button = self.button(QStyle.SP_DirOpenIcon, 'Open example', self.open_example)
        examples.addWidget(self.example_choice, 1)
        examples.addWidget(self.apply_example_button)
        examples.addWidget(self.open_example_button)
        toolbar = QHBoxLayout()
        for button in (self.start, self.reimage, self.load_button, self.export, self.cancel):
            toolbar.addWidget(button)
        toolbar.addWidget(self.check_button)
        toolbar.addWidget(self.camera_button)
        toolbar.addWidget(self.continuous_button)
        toolbar.addWidget(self.log_scale)
        toolbar.addStretch()
        right, layout = QWidget(), QVBoxLayout()
        right.setLayout(layout)
        layout.addWidget(self.status)
        layout.addLayout(examples)
        layout.addLayout(toolbar)
        layout.addWidget(self.progress)
        layout.addWidget(self.results, 1)
        layout.addWidget(self.log)
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(scroll)
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([390, 1000])
        main = QVBoxLayout(self)
        main.addWidget(splitter)

    def button(self, icon, text, callback):
        button = QPushButton(self.style().standardIcon(icon), text)
        button.clicked.connect(callback)
        return button

    def set_busy(self, busy):
        self.controls.setEnabled(not busy)
        self.continuous_button.setEnabled(not busy)
        self.camera_controls.setEnabled(not busy)
        for control in (self.engine, self.executable, self.memory_profile, self.memory, self.timeout, self.start, self.load_button):
            control.setEnabled(not busy)
        for control in (self.example_choice, self.apply_example_button, self.open_example_button):
            control.setEnabled(not busy)
        self.cancel.setEnabled(busy)
        self.reimage.setEnabled(not busy and self.loaded_result is not None)
        self.export.setEnabled(not busy and (self.loaded_result is not None or self.continuous_result is not None))
        self.check_button.setEnabled(not busy and self.loaded_result is not None)
        self.camera_button.setEnabled(not busy and self.loaded_result is not None)
        self.progress.setRange(0, 0 if busy else 1)
        self.progress.setValue(0)

    def start_job(self, checked=False, reuse=False):
        if self.process.state() != QProcess.NotRunning:
            return
        try:
            config = self.controls.config()
            if not Path(self.engine.text()).is_file():
                raise ValueError('Worker Python executable not found')
            if not reuse and not Path(self.executable.text()).is_file():
                raise ValueError('S4 stripes executable not found')
            if reuse and self.loaded_result is None:
                raise ValueError('Load or calculate a scene before updating optics')
            parent = results_root(self.root)/'grating'
            parent.mkdir(parents=True, exist_ok=True)
            self.output = parent/datetime.now().strftime('scene_%Y%m%d_%H%M%S_%f')
            self.job_kind = 'scene'
            input_path = self.output.with_suffix('.input.json')
            input_path.write_text(json.dumps(config.to_dict(), indent=2), encoding='utf-8')
            self.cancel_file = self.output.with_suffix('.cancel')
            args = ['-m', 'grating_scene.run', '--config', str(input_path), '--output', str(self.output),
                '--executable', self.executable.text(), '--max-memory-gib', str(self.memory.value()),
                '--timeout', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file)]
            if reuse:
                args += ['--reuse-source', str(self.loaded_result)]
            if self.memory_profile.text().strip():
                args += ['--memory-profile', self.memory_profile.text().strip()]
            environment = QProcessEnvironment.systemEnvironment()
            environment.insert('PYTHONIOENCODING', 'utf-8')
            environment.insert('PYTHONUNBUFFERED', '1')
            self.process.setProcessEnvironment(environment)
            self.process.setWorkingDirectory(str(self.root))
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.pending_stage = ''
            self.log.clear()
            self.status.setText('S4 periodic calculation | '+str(self.output))
            self.set_busy(True)
            self.process.start(self.engine.text(), args)
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot start grating scene', str(error))

    def start_validation(self):
        if self.loaded_result is None or self.process.state() != QProcess.NotRunning:
            return
        try:
            load_scene(self.loaded_result)
            if not Path(self.engine.text()).is_file():
                raise ValueError('Worker Python executable not found')
            parent = results_root(self.root)/'grating'
            parent.mkdir(parents=True, exist_ok=True)
            self.output = parent/datetime.now().strftime('checks_%Y%m%d_%H%M%S_%f')
            self.cancel_file = self.output.with_suffix('.cancel')
            self.job_kind = 'validation'
            environment = QProcessEnvironment.systemEnvironment()
            environment.insert('PYTHONIOENCODING', 'utf-8')
            environment.insert('PYTHONUNBUFFERED', '1')
            self.process.setProcessEnvironment(environment)
            self.process.setWorkingDirectory(str(self.root))
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.pending_stage = ''
            self.log.clear()
            self.status.setText('Periodic cutoff/slicing checks | not an isolation certification')
            self.set_busy(True)
            args = ['-m', 'tools.grating_scene_validate',
                '--baseline', str(self.loaded_result), '--output', str(self.output),
                '--timeout', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file),
                '--max-memory-gib', str(self.memory.value()), '--executable', self.executable.text()]
            if self.memory_profile.text().strip():
                args += ['--memory-profile', self.memory_profile.text().strip()]
            self.process.start(self.engine.text(), args)
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot start convergence checks', str(error))

    def start_camera(self):
        if self.loaded_result is None or self.process.state() != QProcess.NotRunning:
            return
        try:
            load_scene(self.loaded_result)
            if not Path(self.engine.text()).is_file():
                raise ValueError('Worker Python executable not found')
            settings = self.camera_controls.settings()
            from detector_model import DetectorParameters
            from gray_acquisition import validate_acquisition
            validate_acquisition(DetectorParameters(**settings['detector']), settings['gray_pairs'], settings['gray_seed'])
            parent = results_root(self.root)/'grating'
            parent.mkdir(parents=True, exist_ok=True)
            self.output = parent/datetime.now().strftime('camera_%Y%m%d_%H%M%S_%f')
            settings_path = self.output.with_suffix('.settings.json')
            settings_path.write_text(json.dumps(settings, indent=2), encoding='utf-8')
            self.cancel_file = self.output.with_suffix('.cancel')
            self.job_kind = 'camera'
            environment = QProcessEnvironment.systemEnvironment()
            environment.insert('PYTHONIOENCODING', 'utf-8')
            environment.insert('PYTHONUNBUFFERED', '1')
            self.process.setProcessEnvironment(environment)
            self.process.setWorkingDirectory(str(self.root))
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.pending_stage = ''
            self.log.clear()
            self.status.setText('Area camera | periodic optical model | no TDI')
            self.set_busy(True)
            self.process.start(self.engine.text(), ['-m', 'grating_scene.camera',
                '--source', str(self.loaded_result), '--output', str(self.output), '--settings', str(settings_path),
                '--timeout', str(self.timeout.value()*60), '--cancel-file', str(self.cancel_file),
                '--max-memory-gib', str(self.memory.value())])
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot start periodic camera', str(error))

    def start_continuous(self):
        if self.process.state() != QProcess.NotRunning:
            return
        from .continuous_view import ContinuousInputs
        from rcwa_prototype.current_radiation import load_source, matching_scenes
        dialog = ContinuousInputs(self)
        if not dialog.exec_():
            return
        try:
            source, reference = Path(dialog.source.text()), Path(dialog.reference.text())
            data, scene, info, _ = load_source(source)
            _, other, ref_info, _ = load_source(reference)
            if (not matching_scenes(scene, other) or not info.get('particle') or ref_info.get('particle')
                    or 'source_electric_fourier' not in data):
                raise ValueError('Matching particle current and bare reference exports required')
            if not Path(self.engine.text()).is_file():
                raise ValueError('Worker Python / EXE not found')
            parent = results_root(self.root)/'grating'
            parent.mkdir(parents=True, exist_ok=True)
            self.output = parent/datetime.now().strftime('continuous_%Y%m%d_%H%M%S_%f')
            self.cancel_file = self.output.with_suffix('.cancel')
            self.job_kind = 'continuous'
            args = ['-m', 'rcwa_prototype.angular_distribution', '--source', str(source),
                '--reference', str(reference), '--output', str(self.output),
                '--polar-samples', str(dialog.polar.value()), '--azimuthal-samples', str(dialog.azimuthal.value()),
                '--reference-halfwidth', str(dialog.halfwidth.value()),
                '--max-memory-gib', str(self.memory.value()), '--timeout', str(self.timeout.value()*60),
                '--cancel-file', str(self.cancel_file)]
            environment = QProcessEnvironment.systemEnvironment()
            environment.insert('PYTHONIOENCODING', 'utf-8')
            environment.insert('PYTHONUNBUFFERED', '1')
            environment.insert('MPLBACKEND', 'Agg')
            self.process.setProcessEnvironment(environment)
            self.process.setWorkingDirectory(str(self.root))
            self.decoder = codecs.getincrementaldecoder('utf-8')('replace')
            self.pending_stage = ''
            self.log.clear()
            self.status.setText('Computing continuous current radiation | '+str(self.output))
            self.set_busy(True)
            self.process.start(self.engine.text(), args)
        except (OSError, ValueError, KeyError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot calculate continuous angles', str(error))

    def read_output(self):
        text = self.decoder.decode(bytes(self.process.readAllStandardOutput()))
        self.log.insertPlainText(text)
        lines = (self.pending_stage+text).split('\n')
        self.pending_stage = lines.pop()
        for line in lines:
            if line.startswith('GRATING_STAGE '):
                self.status.setText(line[len('GRATING_STAGE '):].strip()+' | periodic S4')
            elif line.startswith('CONTINUOUS_RING '):
                self.status.setText('Continuous angle rings '+line[len('CONTINUOUS_RING '):])
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    def finished(self, code, status):
        self.read_output()
        self.set_busy(False)
        if code == 0 and status == QProcess.NormalExit and self.output is not None:
            if self.job_kind == 'continuous':
                self.open_result(self.output)
            elif self.job_kind == 'camera':
                try:
                    from grating_scene.camera import load_camera
                    report, arrays = load_camera(self.output)
                    self.camera_controls.load_settings(report['settings'])
                    self.results.show_camera(report, arrays)
                    self.camera_result = self.output
                    self.continuous_result = None
                    self.status.setText('Camera DN completed | conditional periodic-model SNR | '+str(self.output))
                    self.progress.setValue(1)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    self.status.setText('Camera report could not be loaded: '+str(error))
            elif self.job_kind == 'validation':
                try:
                    report = json.loads((self.output/'validation_report.json').read_text(encoding='utf-8'))
                    self.results.show_validation(report)
                    self.status.setText('Convergence controls: '+report['status']+' | isolation NOT certified')
                    self.progress.setValue(1 if report.get('all_cases_solved', False) else 0)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    self.status.setText('Validation report could not be loaded: '+str(error))
            else:
                self.load_result(self.output)
        else:
            self.status.setText('Calculation cancelled or failed; incomplete data not loaded. Previous result preserved.')

    def process_error(self, error):
        if error == QProcess.FailedToStart:
            self.set_busy(False)
            self.status.setText('Worker failed to start: '+self.process.errorString())

    def cancel_job(self):
        if self.process.state() != QProcess.NotRunning and self.cancel_file is not None:
            self.cancel_file.touch()
            self.cancel.setEnabled(False)
            self.status.setText('Cancelling owned worker and native S4 process')

    def choose_result(self):
        folder = QFileDialog.getExistingDirectory(self, 'Open grating scene, camera or continuous radiation')
        if folder:
            self.open_result(folder)

    def open_result(self, folder):
        if (Path(folder)/'continuous_report.json').is_file():
            try:
                from rcwa_prototype.angular_distribution import load_distribution
                report, arrays = load_distribution(folder)
                self.results.continuous_view.show_distribution(report, arrays)
                self.results.setCurrentWidget(self.results.continuous_view)
                self.continuous_result = Path(folder).resolve()
                self.export.setEnabled(True)
                self.progress.setValue(1)
                self.status.setText('Continuous current radiation | isolation unchecked | '+str(folder))
                return True
            except (OSError, ValueError, KeyError, TypeError) as error:
                QMessageBox.warning(self, 'Cannot load continuous angles', str(error))
                return False
        if not (Path(folder)/'camera_report.json').is_file():
            return self.load_result(folder)
        try:
            from grating_scene.camera import load_camera
            report, arrays = load_camera(folder)
            if not self.load_result(Path(folder)/'optical_scene'):
                return False
            self.camera_controls.load_settings(report['settings'])
            self.camera_result = Path(folder).resolve()
            self.results.show_camera(report, arrays)
            self.status.setText('Periodic mixed-scene camera example | optical convergence NOT certified | '+str(folder))
            return True
        except (OSError, ValueError, KeyError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot load camera result', str(error))
            return False

    def apply_example(self):
        try:
            config, camera = example_inputs(self.root, self.example_choice.currentData())
            self.controls.load_config(config)
            if camera:
                self.camera_controls.load_settings(camera)
            self.status.setText('Mixed-scene example parameters applied | displayed result unchanged | convergence unchecked')
        except (OSError, ValueError, KeyError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot apply grating example', str(error))

    def open_example(self):
        self.open_result(example_folder(self.root, self.example_choice.currentData()))

    def load_result(self, folder):
        try:
            config, report, arrays = load_scene(folder)
            self.controls.load_config(config)
            self.results.set_scene(config, report, arrays)
            self.loaded_result = Path(folder).resolve()
            self.camera_result = None
            self.continuous_result = None
            self.reimage.setEnabled(True)
            self.export.setEnabled(True)
            self.check_button.setEnabled(True)
            self.camera_button.setEnabled(True)
            self.progress.setValue(1)
            self.status.setText('S4 periodic model | cutoff/isolation unchecked | ideal relay, not camera | '+str(self.loaded_result))
            return True
        except (OSError, ValueError, KeyError, TypeError) as error:
            QMessageBox.warning(self, 'Cannot load grating scene', str(error))
            return False

    def export_result(self):
        if self.loaded_result is None and self.continuous_result is None:
            return
        source = self.continuous_result or self.camera_result or self.loaded_result
        name, _ = QFileDialog.getSaveFileName(self, 'Export scene with original native data',
            str(source.parent/(source.name+'.zip')), 'ZIP (*.zip)')
        if not name:
            return
        try:
            target = Path(name).with_suffix('.zip').resolve()
            if is_within(target, source):
                raise ValueError('Export the ZIP outside its own source folder')
            if self.continuous_result:
                from rcwa_prototype.angular_distribution import load_distribution
                load_distribution(self.continuous_result)
            else:
                load_scene(self.loaded_result)
            if self.camera_result and not self.continuous_result:
                from grating_scene.camera import load_camera
                load_camera(self.camera_result)
            shutil.make_archive(str(target.with_suffix('')), 'zip', root_dir=source)
            self.status.setText('Exported original sample/reference fields, inputs and derived data: '+str(target))
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, 'Cannot export scene', str(error))

    def wait_for_shutdown(self):
        if self.process.state() == QProcess.NotRunning:
            return True
        self.cancel_job()
        return self.process.waitForFinished(25000)
