"""Native S4 total E/H observations on three explicitly positioned planes."""
import numpy as np


def observation_planes(config):
    if not config.near_enabled:
        return {}
    s, n = config.scene, config.near_samples
    x = np.linspace(-s.lx_um/2, s.lx_um/2, n)
    y = np.linspace(-s.ly_um/2, s.ly_um/2, n)
    z = np.linspace(-.9*s.depth_um, s.cap_top_um+config.near_air_height_um, n)
    a, b = np.meshgrid(x, y)
    xy = np.stack((a, b, np.full_like(a, config.near_xy_z_um)), axis=-1)
    a, b = np.meshgrid(x, z)
    xz = np.stack((a, np.zeros_like(a), b), axis=-1)
    a, b = np.meshgrid(y, z)
    yz = np.stack((np.zeros_like(a), a, b), axis=-1)
    return dict(xy=xy, xz=xz, yz=yz)


def unpack_planes(reference, sample, planes):
    if not planes:
        return {}
    points = np.concatenate([p.reshape(-1, 3) for p in planes.values()])
    for data in (reference, sample):
        np.testing.assert_allclose(data['near_points_um'], points, rtol=0., atol=1e-13)
    output, begin = {}, 0
    for name, grid in planes.items():
        count = grid.shape[0]*grid.shape[1]
        output[name+'_points_um'] = grid
        for label, data in (('reference', reference), ('sample', sample)):
            for key, raw in (('electric', 'near_electric'), ('magnetic', 'near_magnetic_impedance_scaled')):
                output[name+'_'+label+'_'+key] = data[raw][begin:begin+count].reshape(grid.shape)
        for key in ('electric', 'magnetic'):
            output[name+'_delta_'+key] = output[name+'_sample_'+key]-output[name+'_reference_'+key]
        begin += count
    return output


def air_modal_closure(data, config, tolerance=1e-7):
    """Reconstruct exterior E/H from all reflected orders plus the incident wave."""
    if 'near_points_um' not in data:
        return dict(sampled=False)
    s = config.scene
    points = data['near_points_um']
    keep = points[:, 2] > s.cap_top_um+1e-10
    points = points[keep]
    if not len(points):
        return dict(sampled=False)
    k0 = 2*np.pi/s.wavelength_um
    theta = np.deg2rad(s.incidence_deg)
    kx = 2*np.pi*data['orders'][:, 0]/s.lx_um
    ky = k0*np.sin(theta)+2*np.pi*data['orders'][:, 1]/s.ly_um
    q = np.sqrt(k0*k0-kx*kx-ky*ky+0j)
    wavevectors = np.column_stack((kx, ky, q))
    electric = data['electric_reflected_z0']
    magnetic = np.cross(wavevectors/k0, electric)
    incoming_k = np.array([0., np.sin(theta), -np.cos(theta)])
    incoming_e = np.array([-1., 0., 0.]) if config.polarization == 's' else np.array([0., np.cos(theta), np.sin(theta)])
    incoming_h = np.cross(incoming_k, incoming_e)
    predicted = np.empty((len(points), 6), complex)
    for begin in range(0, len(points), 64):
        xyz = points[begin:begin+64]
        phase = np.exp(1j*xyz@wavevectors.T)
        incident = np.exp(1j*k0*(xyz@incoming_k))[:, None]
        predicted[begin:begin+64, :3] = phase@electric+incident*incoming_e
        predicted[begin:begin+64, 3:] = phase@magnetic+incident*incoming_h
    actual = np.column_stack((data['near_electric'][keep], data['near_magnetic_impedance_scaled'][keep]))
    errors = {}
    for label, subset in (('electric', slice(0, 3)), ('magnetic', slice(3, 6))):
        norm = np.sum(abs(actual[:, subset])**2)
        errors[label+'_relative_l2'] = float(np.sqrt(np.sum(abs(predicted[:, subset]-actual[:, subset])**2)/max(norm, 1e-300)))
    passed = all(np.isfinite(value) and value < tolerance for value in errors.values())
    return dict(sampled=True, samples=len(points), passed=passed, tolerance=tolerance, **errors,
        scope='Same finite-basis identity above the entire structure; not truncation or isolated-defect certification.')
