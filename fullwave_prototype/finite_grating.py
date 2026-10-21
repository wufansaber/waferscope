"""Finite 3D etched gratings with real line ends; no periodic boundary."""
from dataclasses import dataclass, asdict

import numpy as np


@dataclass(frozen=True)
class FiniteGrating:
    trench_count: int = 2
    pitch_um: float = .5
    top_width_um: float = .25
    bottom_width_um: float = .25
    depth_um: float = .1
    length_um: float = .8
    margin_um: float = .15
    air_height_um: float = .35
    substrate_depth_um: float = .18
    corner_mesh_um: float = 0.

    def validate(self):
        if isinstance(self.trench_count, bool) or not isinstance(self.trench_count, int) or not 1 <= self.trench_count <= 32:
            raise ValueError('One to 32 finite trenches required')
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite grating dimensions required')
        dimensions = (self.pitch_um, self.top_width_um, self.bottom_width_um, self.depth_um,
            self.length_um, self.margin_um, self.air_height_um, self.substrate_depth_um)
        if min(dimensions) <= 0 or self.corner_mesh_um < 0:
            raise ValueError('Positive dimensions and nonnegative corner size required')
        if max(self.top_width_um, self.bottom_width_um) >= self.pitch_um:
            raise ValueError('Trenches must leave silicon between adjacent grooves')
        if self.depth_um >= self.substrate_depth_um:
            raise ValueError('Trench must remain above the lower PML')

    @property
    def centers(self):
        return (np.arange(self.trench_count)-(self.trench_count-1)/2)*self.pitch_um

    @property
    def extent_y_um(self):
        return (self.trench_count-1)*self.pitch_um/2+max(self.top_width_um, self.bottom_width_um)/2

    @property
    def half_width_um(self):
        return max(self.length_um/2, self.extent_y_um)+self.margin_um

    @property
    def trench_volume_um3(self):
        return self.trench_count*self.length_um*self.depth_um*(self.top_width_um+self.bottom_width_um)/2

    def surface_bounds(self, fraction=.5):
        if not 0 < fraction < 1:
            raise ValueError('Surface clearance fraction must lie strictly between zero and one')
        x = self.length_um/2+(self.half_width_um-self.length_um/2)*fraction
        y = self.extent_y_um+(self.half_width_um-self.extent_y_um)*fraction
        return [[-x, x], [-y, y],
            [-self.depth_um-(self.substrate_depth_um-self.depth_um)*fraction,
             self.air_height_um*fraction]]


def trench_vertices(g, center):
    x = g.length_um/2
    return np.array([[-x, center-g.top_width_um/2, 0], [x, center-g.top_width_um/2, 0],
        [x, center+g.top_width_um/2, 0], [-x, center+g.top_width_um/2, 0],
        [-x, center-g.bottom_width_um/2, -g.depth_um], [x, center-g.bottom_width_um/2, -g.depth_um],
        [x, center+g.bottom_width_um/2, -g.depth_um], [-x, center+g.bottom_width_um/2, -g.depth_um]])


def build_mesh(p, g):
    import ngsolve as ng
    from netgen.occ import Box, Glue, OCCGeometry, Face, Wire, Segment
    g.validate()
    if p.polarization != 'p' or p.beam_profile != 'plane':
        raise ValueError('PMC half-domain currently requires p-polarized plane illumination')
    scale = 1000.
    h, top, bottom, t = g.half_width_um, g.air_height_um, g.substrate_depth_um, p.pml_thickness_um
    def box(lo, hi):
        return Box(tuple(scale*np.asarray(lo)), tuple(scale*np.asarray(hi)))
    outer = box((0, -h-t, -bottom-t), (h+t, h+t, top+t))
    outer.faces.name = 'outer'
    inner = box((0, -h, -bottom), (h, h, top))
    lower = box((-2*(h+t),)*3, (2*(h+t), 2*(h+t), 0))
    grooves = []
    for center in g.centers:
        # A yz trapezoid extruded along the finite x length retains both line ends.
        vertices = scale*trench_vertices(g, center)[[0, 3, 7, 4]]
        wire = Wire([Segment(tuple(a), tuple(b)) for a, b in zip(vertices, np.roll(vertices, -1, axis=0))])
        grooves.append(Face(wire).Extrude((scale*g.length_um, 0, 0)))
    pit = grooves[0]
    for groove in grooves[1:]:
        pit = pit+groove
    pit = pit*inner
    regions = [(inner-lower, 'air', p.air_mesh_um), ((inner*lower)-pit, 'substrate', p.substrate_mesh_um),
        (pit, 'pit', min(.04, p.air_mesh_um)),
        ((outer-inner)-lower, 'pml_air', p.air_mesh_um*p.pml_mesh_scale),
        ((outer-inner)*lower, 'pml_substrate', p.substrate_mesh_um*p.pml_mesh_scale)]
    for shape, name, size in regions:
        shape.mat(name)
        shape.maxh = size*scale
        if name == 'pit':
            shape.faces.maxh = size*scale
    shape = Glue([r[0] for r in regions])
    for face in shape.faces:
        if abs(face.center.x) < 1e-8:
            face.name = 'symmetry'
    if g.corner_mesh_um:
        for edge in shape.edges:
            c = np.array(tuple(edge.center))/scale
            if abs(c[0]) <= g.length_um/2+1e-9 and -g.depth_um-1e-9 <= c[2] <= 1e-9 and abs(c[1]) <= g.extent_y_um+1e-9:
                edge.maxh = g.corner_mesh_um*scale
    native = OCCGeometry(shape).GenerateMesh(maxh=p.air_mesh_um*scale)
    native.Scale(1/scale)
    native.SetGeometry(OCCGeometry(shape.Scale((0, 0, 0), 1/scale)))
    mesh = ng.Mesh(native)
    mesh.Curve(max(3, p.order+1))
    stretch = ng.pml.Cartesian(mins=(0, -h, -bottom), maxs=(h, h, top), alpha=1j*p.pml_strength)
    for name in ('pml_air', 'pml_substrate'):
        mesh.SetPML(stretch, definedon=name)
    volume = float(ng.Integrate(1, mesh, definedon=mesh.Materials('pit')))
    if not np.isclose(2*volume, g.trench_volume_um3, rtol=1e-8):
        raise RuntimeError('3D finite trench volume mismatch')
    return mesh
