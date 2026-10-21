"""One sphere on an infinite line grating, repeated on an explicit supercell."""
from dataclasses import asdict, dataclass
import math

import numpy as np


@dataclass(frozen=True)
class Scene:
    wavelength_um: float = .266
    incidence_deg: float = 75.
    pitch_um: float = .5
    ridge_width_um: float = .25
    depth_um: float = .1
    substrate_n: float = 5.57
    substrate_k: float = .387
    particle_n: float = 1.59
    particle_k: float = 0.
    diameter_um: float = .1
    gap_um: float = .01
    lx_um: float = .5
    periods_y: int = 2
    slices: int = 8
    grid_x: int = 128
    grid_y: int = 256
    modes: int = 121
    factorization: str = 'laurent'
    kind: str = 'grating'

    def validate(self):
        numeric = [v for v in asdict(self).values() if isinstance(v, (float, int))]
        if not all(math.isfinite(v) for v in numeric):
            raise ValueError('Finite inputs required')
        if min(self.wavelength_um, self.pitch_um, self.depth_um, self.substrate_n,
               self.particle_n, self.diameter_um, self.lx_um) <= 0:
            raise ValueError('Positive dimensions and indices required')
        if not 0 < self.ridge_width_um <= self.pitch_um or not 0 <= self.incidence_deg < 90:
            raise ValueError('Invalid ridge width or incidence')
        if min(self.substrate_k, self.particle_k, self.gap_um) < 0:
            raise ValueError('Loss and particle gap must be nonnegative')
        for key in ('periods_y', 'slices', 'grid_x', 'grid_y', 'modes'):
            v = getattr(self, key)
            if isinstance(v, bool) or not isinstance(v, int) or v < 1:
                raise ValueError(key+' must be a positive integer')
        if self.kind not in ('flat', 'grating') or self.factorization not in ('pol', 'laurent'):
            raise ValueError('Unsupported scene or Fourier factorization')

    @property
    def ly_um(self):
        return self.periods_y*self.pitch_um

    @property
    def cap_top_um(self):
        return self.gap_um+self.diameter_um

    @property
    def eps_substrate(self):
        return complex(self.substrate_n, self.substrate_k)**2

    @property
    def eps_particle(self):
        return complex(self.particle_n, self.particle_k)**2


def layers(scene, particle):
    """Layers descend from the sphere top; z=0 is the unetched ridge top."""
    scene.validate()
    if particle and min(scene.lx_um, scene.ly_um) <= scene.diameter_um:
        raise ValueError('Particle must not touch its periodic copies')
    x = np.arange(scene.grid_x)*scene.lx_um/scene.grid_x-scene.lx_um/2
    y = np.arange(scene.grid_y)*scene.ly_um/scene.grid_y-scene.ly_um/2
    xx, yy = np.meshgrid(x, y, indexing='ij')
    result = []
    dz = scene.diameter_um/scene.slices
    radius = scene.diameter_um/2
    volume = 0.
    for i in range(scene.slices):
        z = scene.cap_top_um-(i+.5)*dz
        epsilon = np.ones(xx.shape, complex)
        if particle:
            cross_r2 = max(0., radius**2-(z-scene.gap_um-radius)**2)
            mask = xx**2+yy**2 < cross_r2
            epsilon[mask] = scene.eps_particle
            volume += mask.sum()*scene.lx_um*scene.ly_um/xx.size*dz
        result.append((dz, epsilon))
    if scene.gap_um:
        result.append((scene.gap_um, None))
    offset = (yy+scene.pitch_um/2) % scene.pitch_um-scene.pitch_um/2
    ridge = abs(offset) < scene.ridge_width_um/2
    epsilon = np.full(xx.shape, scene.eps_substrate, complex)
    if scene.kind == 'grating':
        epsilon[~ridge] = 1.
    result.append((scene.depth_um, epsilon))
    return result, dict(particle_present=particle, voxel_particle_volume_um3=float(volume),
        exact_particle_volume_um3=4*np.pi*radius**3/3 if particle else 0.,
        scope=('Infinite straight grooves with one sphere per supercell, not a finite two-trench patch.'
               if scene.kind == 'grating' else 'Flat wafer with one sphere per supercell.'))


def matrix_memory_estimate(modes, patterned_layers):
    """Planning allowance, not a peak-memory guarantee or convergence criterion."""
    one = 16*(2*modes)**2
    return dict(modes=modes, complex_matrix_dimension=2*modes,
        one_complex_matrix_mib=one/2**20,
        planning_gib=(one*(5*patterned_layers+30))/2**30+.4,
        note='Complex128 dense matrices; includes a workspace allowance. Measure actual peak before extrapolating.')
