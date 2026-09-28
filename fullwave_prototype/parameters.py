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

    def validate(self):
        for key, value in asdict(self).items():
            if key not in ('polarization', 'direct_solver') and not math.isfinite(value):
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
