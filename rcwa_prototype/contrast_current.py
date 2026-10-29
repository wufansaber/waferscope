"""Frozen local contrast-current data and reciprocal radiation postprocessing."""
import numpy as np
from scipy.special import j1


def sample_current(model, scene, analytic_spheres, z_order):
    nodes, weights = np.polynomial.legendre.leggauss(z_order)
    k0 = 2*np.pi/scene.wavelength_um
    phase = np.exp(1j*(model.kx*scene.lx_um/2+model.ky*scene.ly_um/2)-
                   1j*k0*np.cos(np.deg2rad(scene.incidence_deg))*scene.cap_top_um-
                   1j*k0*np.sin(np.deg2rad(scene.incidence_deg))*scene.ly_um/2)
    fields, z, wz, radii, indices = [], [], [], [], []
    for layer, radius in analytic_spheres.items():
        thickness = model.thickness_list[layer]
        offsets = (nodes+1)*thickness/2
        data = np.asarray(model.Solve_FieldFourier(layer, offsets))[:, 0].transpose(0, 2, 1).copy()
        data[:, :, 2] *= -1
        fields.extend(data*phase[None, :, None])
        z.extend(scene.cap_top_um-sum(model.thickness_list[:layer])-offsets)
        wz.extend(weights*thickness/2)
        radii.extend([radius]*z_order)
        indices.extend([layer]*z_order)
    return dict(source_electric_fourier=np.asarray(fields, complex).reshape(-1, model.nG, 3),
        source_kxy_per_um=np.column_stack((model.kx.real, model.ky.real)),
        source_z_um=np.asarray(z), source_weight_um=np.asarray(wz),
        source_radius_um=np.asarray(radii), source_layer=np.asarray(indices, dtype=int),
        source_z_order=np.asarray(z_order))


def disk_transform(qx, qy, radius):
    u = radius*np.hypot(qx, qy)
    value = np.ones_like(u)
    keep = u != 0
    value[keep] = 2*j1(u[keep])/u[keep]
    return np.pi*radius**2*value


def project_ray(source, scene, direction, reference_orders=None, reference_halfwidth=20):
    from dataclasses import replace
    from .reciprocal_background import spectrum
    if isinstance(reference_halfwidth, bool) or not isinstance(reference_halfwidth, int) or reference_halfwidth < 1:
        raise ValueError('Positive integer reciprocal-background halfwidth required')
    if reference_orders is None:
        ref_scene = replace(scene, periods_y=1)
        orders = np.column_stack((np.zeros(2*reference_halfwidth+1, int),
                                  np.arange(-reference_halfwidth, reference_halfwidth+1)))
        bloch = None
    else:
        ref_scene, orders = scene, reference_orders
        k0 = 2*np.pi/scene.wavelength_um
        bloch = np.array([0., -k0*np.sin(np.deg2rad(scene.incidence_deg))])
    waves = spectrum(ref_scene, direction, orders, bloch)
    kxy = source['source_kxy_per_um']
    qx = kxy[:, 0, None]+waves['kxy'][None, :, 0]
    qy = kxy[:, 1, None]+waves['kxy'][None, :, 1]
    overlap = np.zeros(2, complex)
    for layer in np.unique(source['source_layer']):
        take = source['source_layer'] == layer
        phase = np.exp(1j*source['source_z_um'][take, None]*waves['kz'][None, :])
        mixed = np.einsum('zgk,zr,z->grk', source['source_electric_fourier'][take],
                           phase, source['source_weight_um'][take], optimize=True)
        disk = disk_transform(qx, qy, source['source_radius_um'][take][0])
        # Lorentz reciprocity is bilinear, not a complex-conjugated overlap.
        overlap += np.einsum('gr,grk,rkp->p', disk, mixed, waves['electric'], optimize=True)
    factor = (2*np.pi/scene.wavelength_um)**2/(4*np.pi)*(scene.eps_particle-1)
    return waves['polarizations']@(factor*overlap)


def modal_reconstruction(source, scene, reference):
    from .comparison import field_change, modal_power
    keep = np.flatnonzero(source['propagating'])
    reconstructed = np.zeros_like(source['electric_reflected_z0'])
    k0 = 2*np.pi/scene.wavelength_um
    area = scene.lx_um*scene.ly_um
    for index in keep:
        g = source['orders'][index]
        sector = (source['orders'][:, 0] == g[0]) & (
            source['orders'][:, 1] % scene.periods_y == g[1] % scene.periods_y)
        far = project_ray(source, scene, source['directions'][index], -source['orders'][sector])
        reconstructed[index] = 2j*np.pi/(area*k0*source['directions'][index, 2])*far
    direct = dict(source, electric_reflected_z0=source['electric_reflected_z0']-
                  reference['electric_reflected_z0'])
    rebuilt = dict(source, electric_reflected_z0=reconstructed)
    result = field_change(rebuilt, direct, scene.incidence_deg)
    a = modal_power(reconstructed, source['directions'], source['propagating'], scene.incidence_deg)
    b = modal_power(direct['electric_reflected_z0'], source['directions'], source['propagating'], scene.incidence_deg)
    result.update(reconstructed_difference_field_efficiency=float(a.sum()),
        direct_difference_field_efficiency=float(b.sum()),
        relative_total_power_difference=float(a.sum()/b.sum()-1) if b.sum() else None,
        caveat='Same finite-basis implementation identity, not truncation convergence or an isolated-defect validation.')
    return reconstructed, result
