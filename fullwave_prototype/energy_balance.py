"""Independent finite-volume Poynting balance; never equate DSC with total energy."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .background_sampling import sampled_background, electric_background


def sphere_points(radius, order):
    nodes, weights = leggauss(order)
    phi = np.arange(2*order)*np.pi/order
    directions, areas = [], []
    for sign in (1, -1):
        mu = sign*(nodes+1)/2
        mm, pp = np.meshgrid(mu, phi, indexing='ij')
        radial = np.sqrt(1-mm*mm)
        directions.append(np.stack((radial*np.cos(pp), radial*np.sin(pp), mm), axis=-1).reshape(-1, 3))
        areas.append(np.broadcast_to(weights[:, None]/2*np.pi/order*radius**2, mm.shape).ravel())
    return np.concatenate(directions), np.concatenate(areas)


def surface_flux(mesh, field, p, radius, order):
    import ngsolve as ng
    normals, areas = sphere_points(radius, order)
    points = normals*radius
    locations = mesh(points[:, 0], points[:, 1], points[:, 2])
    e, h = sampled_background(points, p, magnetic=True)
    e += np.asarray(field(locations)).reshape(-1, 3)
    h += np.asarray(ng.curl(field)(locations)).reshape(-1, 3)/(1j*p.k0)
    flux = np.sum(np.real(np.cross(e, h.conj()))*normals, axis=-1)*areas
    return dict(outward=float(flux.sum()), incoming=float(-flux[flux < 0].sum()),
                outgoing=float(flux[flux > 0].sum()))


def absorption(mesh, field, p, radius, order):
    import ngsolve as ng
    coordinates = ng.CF((ng.x, ng.y, ng.z))
    jacobian = ng.specialcf.JacobianMatrix(3)
    rules, points, values, weighted_loss = {}, [], [], []
    losses = {'substrate': p.eps_substrate.imag, 'bump': p.eps_substrate.imag,
              'particle': p.eps_particle.imag}
    for element in mesh.Elements(ng.VOL):
        loss = losses.get(element.mat, 0)
        if loss <= 0:
            continue
        if element.type not in rules:
            rules[element.type] = ng.IntegrationRule(element.type, order)
        rule = rules[element.type]
        mapped = mesh.GetTrafo(element)(rule)
        xyz = np.asarray(coordinates(mapped)).reshape(-1, 3)
        selected = np.sum(xyz*xyz, axis=1) < radius**2
        if not selected.any():
            continue
        det = np.linalg.det(np.asarray(jacobian(mapped)).reshape(-1, 3, 3))
        weights = abs(det)*np.array([ip.weight for ip in rule])
        points.append(xyz[selected])
        values.append(np.asarray(field(mapped)).reshape(-1, 3)[selected])
        weighted_loss.append(p.k0*loss*weights[selected])
    if not points:
        return 0.
    xyz = np.concatenate(points)
    e = np.concatenate(values)+electric_background(xyz, p)
    return float(np.sum(np.concatenate(weighted_loss)*np.sum(abs(e)**2, axis=1)))


def check_balance(mesh, field, p, irradiance, tolerance=.01):
    radius = .95*p.domain_radius_um
    coarse = surface_flux(mesh, field, p, radius, 24)
    fine = surface_flux(mesh, field, p, radius, 40)
    a5 = absorption(mesh, field, p, radius, 5)
    a9 = absorption(mesh, field, p, radius, 9)
    scale = max(fine['incoming'], a9, 1e-300)
    closure = abs(fine['outward']+a9)/scale
    angular = abs(fine['outward']-coarse['outward'])/scale
    volume = abs(a9-a5)/scale
    return dict(control_radius_um=radius, incoming_w=fine['incoming']*irradiance,
        outgoing_w=fine['outgoing']*irradiance, absorbed_w=a9*irradiance,
        signed_balance_w=(fine['outward']+a9)*irradiance, relative_closure_error=closure,
        surface_quadrature_relative_change=angular, clipped_volume_quadrature_relative_change=volume,
        target=tolerance, target_met=bool(max(closure, angular, volume) <= tolerance),
        scope='Total incident+reflected+scattered field and material absorption inside a sphere, excluding PML. '
              'Clipped element quadrature is approximate; volume and surface checks must both pass. '
              'This is not hemisphere DSC=laser power or a certified FEM discretization error bound.')
