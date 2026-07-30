"""Simulation parameter models and validation."""

from dataclasses import dataclass
from math import floor


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
    channel_numerical_apertures: tuple = None

    def __post_init__(self):
        channel_na = self.channel_numerical_apertures
        if channel_na is None:
            channel_na = (self.numerical_aperture,) * 3
        else:
            channel_na = tuple(channel_na)
        object.__setattr__(self, "channel_numerical_apertures", channel_na)
        object.__setattr__(self, "numerical_aperture", channel_na[0])

    @property
    def tilt_deg(self):
        return self.tilt_degrees[0]

    def validate(self):
        if len(self.channel_numerical_apertures) != 3:
            raise ValueError("Exactly three objective NA values are required.")
        for numerical_aperture in self.channel_numerical_apertures:
            if not 0.0 <= numerical_aperture < 1.0:
                raise ValueError("Objective NA values must be in [0, 1).")
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


@dataclass(frozen=True)
class ParticleSizeSweepParameters:
    start_diameter_nm: float = 30.0
    end_diameter_nm: float = 200.0
    step_diameter_nm: float = 10.0
    maximum_points: int = 200

    @property
    def point_count(self):
        if self.step_diameter_nm <= 0.0 or self.end_diameter_nm < self.start_diameter_nm:
            return 0
        span = self.end_diameter_nm - self.start_diameter_nm
        return int(floor(span / self.step_diameter_nm + 1e-9)) + 1

    @property
    def diameters_nm(self):
        return tuple(
            self.start_diameter_nm + index * self.step_diameter_nm
            for index in range(self.point_count)
        )

    def validate(self):
        if self.start_diameter_nm <= 0.0 or self.end_diameter_nm <= 0.0:
            raise ValueError("Sweep particle diameters must be positive.")
        if self.end_diameter_nm < self.start_diameter_nm:
            raise ValueError("Sweep end diameter must not be smaller than start diameter.")
        if self.step_diameter_nm <= 0.0:
            raise ValueError("Sweep diameter step must be positive.")
        if self.point_count < 2:
            raise ValueError("Particle-size sweep requires at least two points.")
        if self.point_count > self.maximum_points:
            raise ValueError(
                "Particle-size sweep is limited to {} points.".format(
                    self.maximum_points
                )
            )
