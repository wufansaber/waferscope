"""Finite, piecewise-planar roughness patch; flat surface outside the square."""
from dataclasses import dataclass, asdict, field
import math

import numpy as np


@dataclass(frozen=True)
class Roughness:
    width_um: float = .20
    rms_nm: float = 2.0
    grid_points: int = 5
    smoothing_cells: float = .6
    seed: int = 7
    mesh_um: float = .025
    model: str = 'legacy'
    correlation_length_um: float = .05
    psd_exponent: float = 3.0
    taper_fraction: float = .5
    measurement_path: str = ''
    detrend: str = 'plane'
    missing_policy: str = 'reject'
    sampling_grid_levels: list = field(default_factory=list)

    def __post_init__(self):
        object.__setattr__(self, 'sampling_grid_levels', list(self.sampling_grid_levels))

    def validate(self, p, particle_enabled=True):
        for name, value in asdict(self).items():
            if name != 'sampling_grid_levels' and not isinstance(value, str) and not math.isfinite(value):
                raise ValueError(name + ' must be finite')
        if self.sampling_grid_levels:
            from .coupled_surface import joint_nodes
            _, levels = joint_nodes(tuple(self.sampling_grid_levels))
            if self.model not in ('gaussian', 'abc') or self.grid_points not in levels:
                raise ValueError('Coupled refinement requires Gaussian/ABC and a listed height grid')
        if self.model not in ('legacy', 'gaussian', 'abc', 'measured_height', 'measured_psd'):
            raise ValueError('Unknown roughness model')
        if self.model.startswith('measured_'):
            from .experimental_surface import measurement_report
            if not self.measurement_path:
                raise ValueError('Select an experimental measurement file')
            report = measurement_report(self)
            if not np.isclose(self.rms_nm, report['rms_nm'], rtol=1e-8, atol=1e-10):
                raise ValueError('Experimental RMS must come from the measurement, without amplitude rescaling')
        if self.width_um <= 0 or self.rms_nm < 0 or self.mesh_um <= 0 or self.smoothing_cells < 0:
            raise ValueError('Invalid roughness length, RMS or smoothing')
        if self.correlation_length_um <= 0:
            raise ValueError('PSD length must be positive')
        if self.model == 'abc' and not 2.01 <= self.psd_exponent <= 10:
            raise ValueError('ABC PSD exponent C must be in [2.01, 10] for the shared UI')
        if not 0 < self.taper_fraction <= .5:
            raise ValueError('Edge taper fraction must be in (0, 0.5]')
        if not isinstance(self.grid_points, int) or not 4 <= self.grid_points <= 41:
            raise ValueError('Use 4..41 surface grid points per axis')
        if (self.model == 'gaussian' and self.rms_nm > 0 and
                self.width_um/(self.grid_points-1) > self.correlation_length_um/2):
            raise ValueError('Gaussian surface needs at least two grid intervals per correlation length')
        if self.taper_fraction*(self.grid_points-1) < 1-1e-12:
            raise ValueError('Edge taper must span at least one surface grid interval')
        if not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError('seed must be a nonnegative integer')
        _, triangles = height_mesh(self)
        if np.max(np.linalg.norm(triangles, axis=-1)) >= .9*p.domain_radius_um:
            raise ValueError('Rough patch must stay clear of the PML')
        if (particle_enabled and self.rms_nm > 0 and
                p.gap_nm/1000 <= triangles[..., 2].max()+1e-6):
            raise ValueError('Sphere bottom must exceed the highest roughness peak; increase gap_nm')


def raw_height_grid(r):
    from scipy.ndimage import gaussian_filter
    if r.sampling_grid_levels:
        from .coupled_surface import height_grid
        return height_grid(r)
    rng = np.random.default_rng(r.seed)
    raw = rng.standard_normal((r.grid_points, r.grid_points))
    if r.model == 'legacy':
        return gaussian_filter(raw, r.smoothing_cells)
    if r.model == 'abc':
        from .rough_statistics import abc_node_factor
        return (abc_node_factor(r.grid_points, r.width_um, r.correlation_length_um,
                                r.psd_exponent) @ raw.ravel()).reshape(raw.shape)
    if r.model != 'gaussian':
        raise ValueError('Unknown roughness model')
    axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
    difference = (axis[:, None]-axis[None, :])/r.correlation_length_um
    covariance = np.exp(-difference**2)
    factor = np.linalg.cholesky(covariance+1e-12*np.eye(r.grid_points))
    return factor @ raw @ factor.T


def _taper_grid(r):
    distance_from_edge = np.minimum(np.arange(r.grid_points),
                                    np.arange(r.grid_points)[::-1])/(r.grid_points-1)
    edge_phase = np.minimum(distance_from_edge/r.taper_fraction, 1)
    taper_1d = np.sin(np.pi*edge_phase/2)**2
    return np.outer(taper_1d, taper_1d)


def _triangle_indices(n):
    indices = []
    for i in range(n-1):
        for j in range(n-1):
            a, b, c, d = i*n+j, (i+1)*n+j, (i+1)*n+j+1, i*n+j+1
            indices.extend(((a, b, c), (a, c, d)))
    return np.asarray(indices)


def taper_effective_area(r):
    """Integral of the squared piecewise-linear height window over the patch."""
    tapered = _taper_grid(r).ravel()[_triangle_indices(r.grid_points)]
    triangle_area = (r.width_um/(r.grid_points-1))**2/2
    return float(triangle_area*np.sum((np.sum(tapered*tapered, axis=1)+
                                      np.sum(tapered, axis=1)**2)/12))


def height_mesh(r):
    """Build a tapered triangular patch; legacy alone has exact per-sample mean/RMS."""
    axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    taper = _taper_grid(r)
    taper[[0, -1], :] = 0
    taper[:, [0, -1]] = 0
    indices = _triangle_indices(r.grid_points)
    if r.model.startswith('measured_'):
        from .experimental_surface import measured_nodes
        z = measured_nodes(r)*taper
    elif r.model in ('gaussian', 'abc'):
        raw = raw_height_grid(r)
        z = raw*taper*(r.rms_nm/1000)
    else:
        raw = raw_height_grid(r)
        mean_raw = (raw*taper).ravel()[indices].mean()
        mean_taper = taper.ravel()[indices].mean()
        z = (raw-mean_raw/mean_taper)*taper
        ht = z.ravel()[indices]
        second = np.mean((np.sum(ht*ht, axis=1)+np.sum(ht, axis=1)**2)/12)
        if second <= 0:
            raise ValueError('Degenerate surface realization')
        z *= (r.rms_nm/1000)/np.sqrt(second)
    vertices = np.stack((x, y, z), axis=-1)
    return vertices, vertices.reshape(-1, 3)[indices]


def surface_statistics(r, vertices, triangles):
    heights = triangles[..., 2]
    second = np.mean((np.sum(heights*heights, axis=1)+np.sum(heights, axis=1)**2)/12)
    spacing = r.width_um/(r.grid_points-1)
    summary = dict(model=r.model, seed=r.seed, width_um=r.width_um,
                   grid_points=r.grid_points, grid_spacing_um=spacing,
                   taper_fraction=r.taper_fraction,
                   taper_effective_area_um2=taper_effective_area(r),
                   triangle_mean_nm=float(heights.mean()*1000),
                   triangle_rms_nm=float(np.sqrt(second)*1000),
                   input_rms_nm=r.rms_nm,
                   minimum_nm=float(vertices[..., 2].min()*1000),
                   maximum_nm=float(vertices[..., 2].max()*1000),
                   sampling_grid_levels=list(r.sampling_grid_levels),
                   scope='One tapered finite realization; not an ensemble-averaged wafer PSD.')
    if r.model in ('gaussian', 'abc'):
        from .rough_statistics import nyquist_variance_fraction
        length = r.correlation_length_um
        summary.update(correlation_length_um=length,
                       normalization='Input sigma is stationary-parent RMS; taper and one realization change patch mean/RMS.',
                       nyquist_variance_fraction=nyquist_variance_fraction(
                           r.model, spacing, length, r.psd_exponent))
        if r.model == 'gaussian':
            summary.update(correlation_lengths_across_patch=r.width_um/length,
                           target_covariance='C(r)=sigma^2 exp(-(r/L)^2) before finite-patch taper',
                           target_psd='S(f)=pi sigma^2 L^2 exp(-(pi L |f|)^2) before finite-patch taper',
                           target_neighbor_correlation=float(np.exp(-(spacing/length)**2)))
        else:
            from .rough_statistics import abc_amplitude, abc_correlation
            summary.update(psd_lengths_across_patch=r.width_um/length,
                           psd_exponent=r.psd_exponent,
                           psd_length_definition='B is NIST ABC roll-off length, not Gaussian L',
                           target_psd='S(f)=A/(1+(B |f|)^2)^(C/2) before finite-patch taper',
                           abc_amplitude_um4=abc_amplitude(r.rms_nm, length, r.psd_exponent),
                           target_neighbor_correlation=float(
                               abc_correlation(spacing, length, r.psd_exponent)))
    elif r.model.startswith('measured_'):
        from .experimental_surface import measurement_report
        summary['measurement'] = measurement_report(r)
        summary['normalization'] = 'Measured absolute height/PSD; no RMS rescaling, followed by edge taper.'
    else:
        summary['normalization'] = 'Exact zero mean and projected triangle RMS for each legacy sample.'
    return summary


def surface_solid(r, depth, scale=1.):
    from .rough_geometry import surface_solid as construct
    return construct(r, depth, scale)


def build_rough_mesh(p, roughness, particle_enabled=True, cad_scale=1000.):
    from .rough_geometry import build_rough_mesh as construct
    return construct(p, roughness, particle_enabled, cad_scale)
