"""Finite, piecewise-planar roughness patch; flat surface outside the square."""
from dataclasses import dataclass, asdict
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

    def validate(self, p):
        for name, value in asdict(self).items():
            if name != 'model' and not math.isfinite(value):
                raise ValueError(name + ' must be finite')
        if self.model not in ('legacy', 'gaussian'):
            raise ValueError('Roughness model must be legacy or gaussian')
        if self.width_um <= 0 or self.rms_nm < 0 or self.mesh_um <= 0 or self.smoothing_cells < 0:
            raise ValueError('Invalid roughness length, RMS or smoothing')
        if self.correlation_length_um <= 0:
            raise ValueError('Correlation length must be positive')
        if not isinstance(self.grid_points, int) or not 4 <= self.grid_points <= 17:
            raise ValueError('Use 4..17 surface grid points per axis')
        if (self.model == 'gaussian' and self.rms_nm > 0 and
                self.width_um/(self.grid_points-1) > self.correlation_length_um/2):
            raise ValueError('Gaussian surface needs at least two grid intervals per correlation length')
        if not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError('seed must be a nonnegative integer')
        _, triangles = height_mesh(self)
        if np.max(np.linalg.norm(triangles, axis=-1)) >= .9*p.domain_radius_um:
            raise ValueError('Rough patch must stay clear of the PML')
        if self.rms_nm > 0 and p.gap_nm/1000 <= triangles[..., 2].max()+1e-6:
            raise ValueError('Sphere bottom must exceed the highest roughness peak; increase gap_nm')


def raw_height_grid(r):
    from scipy.ndimage import gaussian_filter
    rng = np.random.default_rng(r.seed)
    raw = rng.standard_normal((r.grid_points, r.grid_points))
    if r.model == 'legacy':
        return gaussian_filter(raw, r.smoothing_cells)
    if r.model != 'gaussian':
        raise ValueError('Unknown roughness model')
    axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
    difference = (axis[:, None]-axis[None, :])/r.correlation_length_um
    covariance = np.exp(-difference**2)
    factor = np.linalg.cholesky(covariance+1e-12*np.eye(r.grid_points))
    return factor @ raw @ factor.T


def height_mesh(r):
    """Build a tapered triangular patch; legacy alone has exact per-sample mean/RMS."""
    axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    taper = np.outer(np.sin(np.linspace(0, np.pi, r.grid_points))**2,
                     np.sin(np.linspace(0, np.pi, r.grid_points))**2)
    taper[[0, -1], :] = 0
    taper[:, [0, -1]] = 0
    raw = raw_height_grid(r)
    indices = []
    n = r.grid_points
    for i in range(n-1):
        for j in range(n-1):
            a, b, c, d = i*n+j, (i+1)*n+j, (i+1)*n+j+1, i*n+j+1
            indices.extend(((a, b, c), (a, c, d)))
    indices = np.asarray(indices)
    if r.model == 'gaussian':
        z = raw*taper*(r.rms_nm/1000)
    else:
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
                   triangle_mean_nm=float(heights.mean()*1000),
                   triangle_rms_nm=float(np.sqrt(second)*1000),
                   input_rms_nm=r.rms_nm,
                   minimum_nm=float(vertices[..., 2].min()*1000),
                   maximum_nm=float(vertices[..., 2].max()*1000),
                   scope='One tapered finite realization; not an ensemble-averaged wafer PSD.')
    if r.model == 'gaussian':
        length = r.correlation_length_um
        summary.update(correlation_length_um=length,
                       target_covariance='C(r)=sigma^2 exp(-(r/L)^2) before finite-patch taper',
                       target_psd='S(f)=pi sigma^2 L^2 exp(-(pi L |f|)^2) before finite-patch taper',
                       normalization='Input sigma is stationary-parent RMS; taper and one realization change patch mean/RMS.',
                       target_neighbor_correlation=float(np.exp(-(spacing/length)**2)),
                       correlation_lengths_across_patch=r.width_um/length)
    else:
        summary['normalization'] = 'Exact zero mean and projected triangle RMS for each legacy sample.'
    return summary


def surface_solid(r, depth):
    """Sew planar triangular faces into a closed OCC solid above z=-depth."""
    from netgen.occ import Vertex, Face, MakePolygon, Sew, Solid
    vertices, triangles = height_mesh(r)

    def face(points):
        corners = [Vertex(tuple(point)) for point in points]
        return Face(MakePolygon(corners+[corners[0]]))

    faces = [face(t) for t in triangles]
    edge = np.concatenate((vertices[:-1, 0], vertices[-1, :-1],
                           vertices[:0:-1, -1], vertices[0, :0:-1]))
    bottom = edge.copy()
    bottom[:, 2] = -depth
    for i in range(len(edge)):
        j = (i+1) % len(edge)
        faces.append(face((edge[j], edge[i], bottom[i], bottom[j])))
    half = r.width_um/2
    faces.append(face(((-half, -half, -depth), (-half, half, -depth),
                       (half, half, -depth), (half, -half, -depth))))
    return Solid(Sew(faces, tolerance=1e-9))


def build_rough_mesh(p, roughness):
    import ngsolve as ng
    from netgen.occ import Sphere, Box, Glue, OCCGeometry
    from .solver import build_mesh
    from .absorbing_boundary import set_boundary
    roughness.validate(p)
    if roughness.rms_nm == 0:
        return build_mesh(p)
    radius = p.domain_radius_um
    outer_r = radius+p.pml_thickness_um
    outer = Sphere((0, 0, 0), outer_r)
    outer.faces.name = 'outer'
    inner = Sphere((0, 0, 0), radius)
    lower = Box((-2*outer_r,)*3, (2*outer_r, 2*outer_r, 0))
    half = roughness.width_um/2
    vertices, _ = height_mesh(roughness)
    depth = max(.03, 2*abs(vertices[..., 2]).max())
    slab = surface_solid(roughness, depth)
    column = Box((-half, -half, -depth), (half, half, 0))
    bump, pit = slab-lower, column-slab
    ball = Sphere((0, 0, (p.diameter_nm/2+p.gap_nm)/1000), p.diameter_nm/2000)
    regions = [(inner-lower-bump-ball, 'air', p.air_mesh_um),
               ((inner*lower)-pit, 'substrate', p.substrate_mesh_um),
               (bump, 'bump', roughness.mesh_um), (pit, 'pit', roughness.mesh_um),
               (ball, 'particle', p.particle_mesh_um),
               ((outer-inner)-lower, 'pml_air', p.air_mesh_um*p.pml_mesh_scale),
               ((outer-inner)*lower, 'pml_substrate', p.substrate_mesh_um*p.pml_mesh_scale)]
    for shape, name, size in regions:
        shape.mat(name)
        shape.maxh = size
        if name in ('bump', 'pit', 'particle'):
            shape.faces.maxh = size
    mesh = ng.Mesh(OCCGeometry(Glue([s for s, _, _ in regions])).GenerateMesh(maxh=p.air_mesh_um))
    mesh.Curve(max(3, p.order+1))
    set_boundary(mesh, p)
    return mesh
