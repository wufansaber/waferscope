"""Extract polarization quadrature directly from curved FEM volume elements."""
import numpy as np

from .solver import background_coefficient, contrast_materials


def polarization_sources(mesh, field, p, particle_enabled=True, order=5, vectorized_background=False):
    import ngsolve as ng
    total = field if vectorized_background else field+background_coefficient(p)
    coordinates = ng.CF((ng.x, ng.y, ng.z))
    jacobian = ng.specialcf.JacobianMatrix(3)
    contrasts = contrast_materials(p, particle_enabled)
    points, sources, contrast_weights = [], [], []
    volumes = {name: 0. for name in contrasts if name in mesh.GetMaterials()}
    rules = {}
    for element in mesh.Elements(ng.VOL):
        contrast = contrasts.get(element.mat, 0)
        if contrast == 0:
            continue
        if element.type not in rules:
            rules[element.type] = ng.IntegrationRule(element.type, order)
        rule = rules[element.type]
        mapped = mesh.GetTrafo(element)(rule)
        determinants = np.linalg.det(np.asarray(jacobian(mapped)).reshape(-1, 3, 3))
        if not np.isfinite(determinants).all() or np.any(determinants == 0):
            raise ValueError('Invalid element Jacobian in source quadrature')
        weights = np.abs(determinants)*np.array([ip.weight for ip in rule])
        volumes[element.mat] += float(weights.sum())
        points.append(np.asarray(coordinates(mapped)).reshape(-1, 3))
        sources.append(contrast*np.asarray(total(mapped)).reshape(-1, 3)*weights[:, None])
        if vectorized_background:
            contrast_weights.append(contrast*weights)
    if not points:
        return np.empty((0, 3)), np.empty((0, 3), complex), volumes
    points, sources = np.concatenate(points), np.concatenate(sources)
    if vectorized_background:
        from .background_sampling import electric_background
        sources += electric_background(points, p)*np.concatenate(contrast_weights)[:, None]
    return points, sources, volumes
