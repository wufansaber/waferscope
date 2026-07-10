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


@dataclass(frozen=True)
class ObjectiveCollectionParameters:
    numerical_aperture: float = 0.5
    focal_length_mm: float = 20.0
    tilt_degrees: tuple = (45.0, -45.0, 0.0)
    grid_size: int = 181
    display_floor_db: float = -35.0

    @property
    def tilt_deg(self):
        return self.tilt_degrees[0]

    def validate(self):
        if not 0.0 < self.numerical_aperture < 1.0:
            raise ValueError("Objective NA must be in the open interval (0, 1).")
        if self.focal_length_mm <= 0.0:
            raise ValueError("Objective focal length must be positive.")
        if len(self.tilt_degrees) != 3:
            raise ValueError("Exactly three objective tilt angles are required.")
        for tilt_deg in self.tilt_degrees:
            if not -89.9 <= tilt_deg <= 89.9:
                raise ValueError(
                    "Objective tilt angles must be between -89.9 and 89.9 degrees."
                )
        if self.grid_size < 51:
            raise ValueError("Objective example grid size must be at least 51.")
