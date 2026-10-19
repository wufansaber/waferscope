"""OCC roughness geometry built at nanometre CAD scale, solved in micrometres."""
import math


def surface_solid(r, depth, scale=1.):
    """Sew the unchanged height triangles into a closed, scaled CAD solid."""
    import numpy as np
    from netgen.occ import Vertex, Face, MakePolygon, Sew, Solid
    from .rough_surface import height_mesh
    vertices, triangles = height_mesh(r)
    vertices, triangles, depth = vertices*scale, triangles*scale, depth*scale

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
    half = r.width_um*scale/2
    faces.append(face(((-half, -half, -depth), (-half, half, -depth),
                       (half, half, -depth), (half, -half, -depth))))
    return Solid(Sew(faces, tolerance=1e-9))


def build_rough_mesh(p, roughness, particle_enabled=True, cad_scale=1000.):
    import ngsolve as ng
    from netgen.occ import Sphere, Box, Glue, OCCGeometry
    from .rough_surface import height_mesh
    from .solver import build_mesh
    from .defects import cad_solid
    from .absorbing_boundary import set_boundary
    if not math.isfinite(cad_scale) or cad_scale <= 0:
        raise ValueError('CAD construction scale must be finite and positive')
    roughness.validate(p, particle_enabled=particle_enabled)
    if roughness.rms_nm == 0:
        return build_mesh(p)
    radius = p.domain_radius_um*cad_scale
    outer_r = radius+p.pml_thickness_um*cad_scale
    outer = Sphere((0, 0, 0), outer_r)
    outer.faces.name = 'outer'
    inner = Sphere((0, 0, 0), radius)
    lower = Box((-2*outer_r,)*3, (2*outer_r, 2*outer_r, 0))
    half = roughness.width_um*cad_scale/2
    vertices, _ = height_mesh(roughness)
    depth = max(.03, 2*abs(vertices[..., 2]).max())*cad_scale
    slab = surface_solid(roughness, depth/cad_scale, cad_scale)
    column = Box((-half, -half, -depth), (half, half, 0))
    bump, pit = slab-lower, column-slab
    air = inner-lower-bump
    if particle_enabled:
        ball = cad_solid(p, cad_scale)
        air = air-ball
    regions = [(air, 'air', p.air_mesh_um),
               ((inner*lower)-pit, 'substrate', p.substrate_mesh_um),
               (bump, 'bump', roughness.mesh_um), (pit, 'pit', roughness.mesh_um),
               ((outer-inner)-lower, 'pml_air', p.air_mesh_um*p.pml_mesh_scale),
               ((outer-inner)*lower, 'pml_substrate', p.substrate_mesh_um*p.pml_mesh_scale)]
    if particle_enabled:
        regions.append((ball, 'particle', p.particle_mesh_um))
    for shape, name, size in regions:
        shape.mat(name)
        shape.maxh = size*cad_scale
        if name in ('bump', 'pit', 'particle'):
            shape.faces.maxh = size*cad_scale
    shape = Glue([s for s, _, _ in regions])
    native = OCCGeometry(shape).GenerateMesh(maxh=p.air_mesh_um*cad_scale)
    # Scale linear vertices and replace CAD geometry before constructing curved maps.
    if cad_scale != 1:
        native.Scale(1/cad_scale)
        native.SetGeometry(OCCGeometry(shape.Scale((0, 0, 0), 1/cad_scale)))
    mesh = ng.Mesh(native)
    mesh.Curve(max(3, p.order+1))
    set_boundary(mesh, p)
    return mesh
