"""H(curl) space with optional lower order strictly inside the PML."""
import numpy as np


def make_space(mesh, p):
    import ngsolve as ng
    options = dict(order=p.order, complex=True, dirichlet='outer')
    if p.pml_order in (0, p.order):
        return ng.HCurl(mesh, **options)
    space = ng.HCurl(mesh, order_policy=ng.ORDER_POLICY.VARIABLE, **options)
    physical_edges, physical_faces, pml_cells = set(), set(), []
    for element in mesh.Elements(ng.VOL):
        lower = element.mat.startswith('pml_')
        if lower and p.pml_inner_fraction > 0:
            center = np.mean([mesh[vertex].point for vertex in element.vertices], axis=0)
            lower = np.linalg.norm(center) > p.domain_radius_um+p.pml_inner_fraction*p.pml_thickness_um
        if lower:
            pml_cells.append(element.nr)
        else:
            physical_edges.update(edge.nr for edge in element.edges)
            physical_faces.update(face.nr for face in element.faces)
    for edge in mesh.edges:
        if edge.nr not in physical_edges:
            space.SetOrder(ng.NodeId(ng.EDGE, edge.nr), p.pml_order)
    for face in mesh.faces:
        if face.nr not in physical_faces:
            space.SetOrder(ng.NodeId(ng.FACE, face.nr), p.pml_order)
    for cell in pml_cells:
        space.SetOrder(ng.NodeId(ng.CELL, cell), p.pml_order)
    space.Update()
    return space
