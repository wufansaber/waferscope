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

    def validate(self, p):
        for name, value in asdict(self).items():
            if not math.isfinite(value):
                raise ValueError(name + ' must be finite')
        if self.width_um <= 0 or self.rms_nm < 0 or self.mesh_um <= 0 or self.smoothing_cells < 0:
            raise ValueError('Invalid roughness length, RMS or smoothing')
        if not isinstance(self.grid_points, int) or not 4 <= self.grid_points <= 17:
            raise ValueError('Use 4..17 surface grid points per axis')
        if not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError('seed must be a nonnegative integer')
        _, triangles = height_mesh(self)
        if np.max(np.linalg.norm(triangles, axis=-1)) >= .9*p.domain_radius_um:
            raise ValueError('Rough patch must stay clear of the PML')
        if self.rms_nm > 0 and p.gap_nm/1000 <= triangles[..., 2].max()+1e-6:
            raise ValueError('Sphere bottom must exceed the highest roughness peak; increase gap_nm')


def height_mesh(r):
    """Normalize projected-area mean/RMS of the actual linear triangles, not nodes."""
    from scipy.ndimage import gaussian_filter
    axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    taper = np.outer(np.sin(np.linspace(0, np.pi, r.grid_points))**2,
                     np.sin(np.linspace(0, np.pi, r.grid_points))**2)
    taper[[0, -1], :] = 0
    taper[:, [0, -1]] = 0
    raw = gaussian_filter(np.random.default_rng(r.seed).standard_normal(x.shape), r.smoothing_cells)
    indices = []
    n = r.grid_points
    for i in range(n-1):
        for j in range(n-1):
            a, b, c, d = i*n+j, (i+1)*n+j, (i+1)*n+j+1, i*n+j+1
            indices.extend(((a, b, c), (a, c, d)))
    indices = np.asarray(indices)
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
