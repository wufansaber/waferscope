"""Explicit scope and numerical controls for the first FEM benchmark."""
from dataclasses import dataclass, asdict
import math


@dataclass(frozen=True)
class Parameters:
    wavelength_nm: float = 355.0
    diameter_nm: float = 100.0
    particle_n: float = 1.59
    particle_k: float = 0.0
    substrate_n: float = 5.57
    substrate_k: float = 0.387
    incidence_deg: float = 75.0
    polarization: str = 'p'
    gap_nm: float = 0.0
    domain_radius_um: float = 0.22
    pml_thickness_um: float = 0.16
    pml_strength: float = 2.0
    pml_power: int = 1
    pml_mesh_scale: float = 1.0
    air_mesh_um: float = 0.08
    substrate_mesh_um: float = 0.04
    particle_mesh_um: float = 0.015
    order: int = 2
    threads: int = 2
    max_dofs: int = 350000
    static_condensation: bool = False
    direct_solver: str = 'umfpack'
    integration_bonus: int = 0
    pml_order: int = 0
    pml_inner_fraction: float = 0.0
    beam_waist_um: float = 0.0
    beam_order: int = 5
    beam_profile: str = 'plane'
    beam_length_um: float = 0.0
    beam_width_um: float = 0.0
    beam_definition: str = 'wafer_gaussian'
    beam_measurement_path: str = ''
    beam_measurement_sha256: str = ''
    beam_phase_assumption: str = 'flat'

    def validate(self):
        from .beam_inputs import DEFINITIONS, fingerprint, measured_envelope
        if self.beam_definition not in DEFINITIONS or self.beam_phase_assumption not in ('flat', 'reject'):
            raise ValueError('Invalid finite illumination definition or phase assumption')
        if self.beam_definition == 'measured_wafer' and self.beam_profile == 'gaussian':
            if not self.beam_measurement_path or fingerprint(self.beam_measurement_path) != self.beam_measurement_sha256:
                raise ValueError('Measured illumination requires its unchanged file and SHA256')
            measured_envelope(self.beam_measurement_path, self.beam_measurement_sha256, self.beam_phase_assumption)
        if self.beam_waist_um < 0 or not isinstance(self.beam_order, int) or self.beam_order not in (3, 5, 7, 9, 15, 21, 31, 41, 61, 81, 101, 151):
            raise ValueError('Beam waist must be nonnegative; supported quadrature orders: 3/5/7/9/15/21/31/41/61.')
        if self.beam_profile not in ('plane', 'gaussian'):
            raise ValueError('Use plane or gaussian illumination')
        if min(self.beam_length_um, self.beam_width_um) < 0:
            raise ValueError('Beam diameters must be nonnegative')
        if self.beam_profile == 'gaussian' and (
                min(self.beam_length_um, self.beam_width_um) <= 0 or self.beam_waist_um != 0):
            raise ValueError('Gaussian footprint requires positive diameters and legacy beam waist = 0')
        for key, value in asdict(self).items():
            if key not in ('polarization', 'direct_solver', 'beam_profile', 'beam_definition',
                           'beam_measurement_path', 'beam_measurement_sha256', 'beam_phase_assumption') and not math.isfinite(value):
                raise ValueError(key + ' must be finite')
        for key in ('wavelength_nm', 'diameter_nm', 'particle_n', 'substrate_n',
                    'domain_radius_um', 'pml_thickness_um', 'pml_strength',
                    'air_mesh_um', 'substrate_mesh_um', 'particle_mesh_um', 'pml_mesh_scale'):
            if getattr(self, key) <= 0:
                raise ValueError(key + ' must be positive')
        if min(self.particle_k, self.substrate_k, self.gap_nm) < 0:
            raise ValueError('Loss and gap must be nonnegative')
        if self.domain_radius_um <= (self.diameter_nm + self.gap_nm)/1000:
            raise ValueError('The particle must be entirely inside the physical domain')
        if not 0 <= self.incidence_deg < 90 or self.polarization not in ('s', 'p'):
            raise ValueError('Use s/p polarization and incidence in [0,90)')
        if self.order not in (1, 2, 3) or not isinstance(self.threads, int) or self.threads < 1:
            raise ValueError('Use order 1/2/3 and a positive integer thread count')
        if not isinstance(self.max_dofs, int) or self.max_dofs < 1:
            raise ValueError('max_dofs must be a positive integer')
        if not isinstance(self.pml_power, int) or self.pml_power not in (1, 2, 3):
            raise ValueError('pml_power must be 1, 2 or 3')
        if not isinstance(self.static_condensation, bool):
            raise ValueError('static_condensation must be boolean')
        if self.direct_solver not in ('umfpack', 'sparsecholesky'):
            raise ValueError('Use umfpack or sparsecholesky direct solver')
        if not isinstance(self.integration_bonus, int) or not 0 <= self.integration_bonus <= 12:
            raise ValueError('integration_bonus must be an integer in 0..12')
        if not isinstance(self.pml_order, int) or not 0 <= self.pml_order <= self.order:
            raise ValueError('pml_order must be 0 (same as physical region) or between 1 and order')
        if not 0 <= self.pml_inner_fraction < 1:
            raise ValueError('pml_inner_fraction must be in [0,1)')

    @property
    def k0(self):
        return 2*math.pi/(self.wavelength_nm/1000)

    @property
    def eps_particle(self):
        return complex(self.particle_n, self.particle_k)**2

    @property
    def eps_substrate(self):
        return complex(self.substrate_n, self.substrate_k)**2
