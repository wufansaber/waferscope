"""Qt worker object for running native SCATMECH calculations off the UI thread."""

import traceback

from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

from .backend import calculate_far_field
from .size_sweep import calculate_particle_size_sweep


class SimulationWorker(QObject):
    progress = pyqtSignal(int)
    completed = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, parameters):
        super().__init__()
        self.parameters = parameters

    @pyqtSlot()
    def run(self):
        try:
            result = calculate_far_field(self.parameters, self.progress.emit)
            self.completed.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())


class ParticleSizeSweepWorker(QObject):
    progress = pyqtSignal(int)
    completed = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, simulation_parameters, objective_parameters, sweep_parameters):
        super().__init__()
        self.simulation_parameters = simulation_parameters
        self.objective_parameters = objective_parameters
        self.sweep_parameters = sweep_parameters

    @pyqtSlot()
    def run(self):
        try:
            result = calculate_particle_size_sweep(
                self.simulation_parameters,
                self.objective_parameters,
                self.sweep_parameters,
                self.progress.emit,
            )
            self.completed.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())
