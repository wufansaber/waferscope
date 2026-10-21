"""Finite etched-line verification geometry, not an infinite periodic wafer."""
from dataclasses import dataclass, asdict

import numpy as np


@dataclass(frozen=True)
class LinePatch:
    pitch_um: float = 1.
    trench_width_um: float = .5
    trench_depth_um: float = .1
    line_length_um: float = 1.4
    half_width_um: float = .9
    air_height_um: float = .30
    substrate_depth_um: float = .18

    def validate(self, p):
        if not all(np.isfinite(v) and v > 0 for v in asdict(self).values()):
            raise ValueError('Positive finite pattern dimensions required')
        if self.trench_width_um >= self.pitch_um:
            raise ValueError('Trenches must leave a central ridge')
        if max(self.line_length_um/2, (self.pitch_um+self.trench_width_um)/2) >= self.half_width_um:
            raise ValueError('Pattern must remain inside the physical box')
        if self.trench_depth_um >= self.substrate_depth_um:
            raise ValueError('Trench must remain above the lower PML')
        if (p.diameter_nm+p.gap_nm)/1000 >= self.air_height_um:
            raise ValueError('Particle must remain below the upper PML')
        if p.diameter_nm/1000 >= self.pitch_um-self.trench_width_um:
            raise ValueError('Particle must fit on the central ridge')

    @property
    def centers(self):
        return (-self.pitch_um/2, self.pitch_um/2)

    @property
    def trench_volume_um3(self):
        return 2*self.line_length_um*self.trench_width_um*self.trench_depth_um


def build_mesh(p, pattern, half_domain=False, patterned=True):
    import ngsolve as ng
    from netgen.occ import Box, Sphere, Glue, OCCGeometry
    pattern.validate(p)
    scale = 1000.
    h, z, d, t = (pattern.half_width_um, pattern.air_height_um,
                  pattern.substrate_depth_um, p.pml_thickness_um)

    def box(lo, hi):
        return Box(tuple(scale*np.asarray(lo)), tuple(scale*np.asarray(hi)))

    if half_domain and (p.polarization != 'p' or p.beam_profile != 'plane'):
        raise ValueError('Half-domain PMC example requires x-invariant p-polarized plane-wave illumination')
    outer = box((0 if half_domain else -h-t, -h-t, -d-t), (h+t, h+t, z+t))
    outer.faces.name = 'outer'
    inner = box((0 if half_domain else -h, -h, -d), (h, h, z))
    lower = box((-2*(h+t),)*3, (2*(h+t), 2*(h+t), 0))
    trenches = [box((-pattern.line_length_um/2, c-pattern.trench_width_um/2, -pattern.trench_depth_um),
                    (pattern.line_length_um/2, c+pattern.trench_width_um/2, 0)) for c in pattern.centers]
    pit = (trenches[0]+trenches[1])*inner
    ball = Sphere((0, 0, scale*(p.gap_nm+p.diameter_nm/2)/1000), scale*p.diameter_nm/2000)
    if half_domain:
        ball = ball*inner
    substrate = (inner*lower)-pit if patterned else inner*lower
    regions = [(inner-lower-ball, 'air', p.air_mesh_um),
               (substrate, 'substrate', p.substrate_mesh_um),
               (ball, 'particle', p.particle_mesh_um),
               ((outer-inner)-lower, 'pml_air', p.air_mesh_um*p.pml_mesh_scale),
               ((outer-inner)*lower, 'pml_substrate', p.substrate_mesh_um*p.pml_mesh_scale)]
    if patterned:
        regions.append((pit, 'pit', min(.04, p.air_mesh_um)))
    for shape, name, size in regions:
        shape.mat(name)
        shape.maxh = size*scale
        if name in ('pit', 'particle'):
            shape.faces.maxh = size*scale
    shape = Glue([shape for shape, _, _ in regions])
    if half_domain:
        # Natural H(curl) boundary n x curl(E)=0 is the p-wave PMC symmetry plane.
        for face in shape.faces:
            if abs(face.center.x) < 1e-8:
                face.name = 'symmetry'
    native = OCCGeometry(shape).GenerateMesh(maxh=p.air_mesh_um*scale)
    native.Scale(1/scale)
    native.SetGeometry(OCCGeometry(shape.Scale((0, 0, 0), 1/scale)))
    mesh = ng.Mesh(native)
    mesh.Curve(max(3, p.order+1))
    transform = ng.pml.Cartesian(mins=(0 if half_domain else -h, -h, -d), maxs=(h, h, z), alpha=1j*p.pml_strength)
    for material in ('pml_air', 'pml_substrate'):
        mesh.SetPML(transform, definedon=material)
    volume = float(ng.Integrate(1, mesh, definedon=mesh.Materials('pit'))) if patterned else 0.
    expected = pattern.trench_volume_um3/(2 if half_domain else 1) if patterned else 0.
    if not np.isclose(volume, expected, rtol=1e-9, atol=1e-14):
        raise RuntimeError('Etched trench material volume mismatch')
    return mesh


def mirror_sources(points, weights):
    """Reconstruct the full vector polarization from an exact PMC half-domain solve."""
    reflection = np.array((-1, 1, 1))
    return np.concatenate((points, points*reflection)), np.concatenate((weights, weights*reflection))


def sample_half_planes(mesh, field, p, planes, state, data):
    from .solver import evaluate_fields
    for plane in planes:
        points = plane.points.copy()
        negative = points[:, 0] < 0
        points[:, 0] = abs(points[:, 0])
        electric, magnetic, _ = evaluate_fields(mesh, field, points, p)
        electric[negative] *= (-1, 1, 1)
        magnetic[negative] *= (1, -1, -1)
        prefix = 'near_'+plane.name+'_'
        data[prefix+'u_um'], data[prefix+'v_um'] = plane.u, plane.v
        data[prefix+'position_um'] = np.asarray(plane.offset)
        shape = (len(plane.v), len(plane.u), 3)
        data[prefix+'e_'+state] = electric.reshape(shape)
        data[prefix+'hz0_'+state] = magnetic.reshape(shape)


def surface_triangles(pattern):
    triangles = []

    def quad(points):
        triangles.extend((np.asarray(points)[[0, 1, 2]], np.asarray(points)[[0, 2, 3]]))

    x, h, depth = pattern.line_length_um/2, pattern.half_width_um, pattern.trench_depth_um
    edges = [-h]+[v for c in pattern.centers for v in
                   (c-pattern.trench_width_um/2, c+pattern.trench_width_um/2)]+[h]
    for a, b in zip(edges[:-1], edges[1:]):
        trench = any(a < c < b for c in pattern.centers)
        z = -depth if trench else 0
        quad(((-x, a, z), (x, a, z), (x, b, z), (-x, b, z)))
        for lo, hi in ((-h, -x), (x, h)):
            quad(((lo, a, 0), (hi, a, 0), (hi, b, 0), (lo, b, 0)))
        if trench:
            for y in (a, b):
                quad(((-x, y, 0), (x, y, 0), (x, y, -depth), (-x, y, -depth)))
            for xx in (-x, x):
                quad(((xx, a, 0), (xx, b, 0), (xx, b, -depth), (xx, a, -depth)))
    return np.asarray(triangles)
