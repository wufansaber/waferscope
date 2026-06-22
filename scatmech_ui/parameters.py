"""Simulation parameter model and validation."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SimulationParameters:
    wafer_material: str = "Silicon (Si)"
    wafer_n: float = 5.57
    wafer_k: float = 0.387
    wafer_thickness_um: float = 725.0
    particle_material: str = "PSL"
    particle_n: float = 1.59
    particle_k: float = 0.0
    particle_diameter_nm: float = 100.0
    wavelength_nm: float = 355.0
    incidence_deg: float = 75.0
    polarization: str = "p"
    angular_step_deg: float = 4.0
    display_floor_db: float = -35.0

    def validate(self):
        if self.wafer_n <= 0 or self.particle_n <= 0:
            raise ValueError("Refractive-index real parts must be positive.")
        if self.wafer_k < 0 or self.particle_k < 0:
            raise ValueError("Extinction coefficients cannot be negative.")
        if self.wafer_thickness_um <= 0:
            raise ValueError("Wafer thickness must be positive.")
        if self.particle_diameter_nm <= 0 or self.wavelength_nm <= 0:
            raise ValueError("Particle diameter and wavelength must be positive.")
        if not 0 <= self.incidence_deg < 90:
            raise ValueError("Incidence angle must be in [0, 90) degrees.")
        if not 2.0 <= self.angular_step_deg <= 15.0:
            raise ValueError("Angular step must be between 2 and 15 degrees.")
