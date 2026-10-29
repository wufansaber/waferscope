"""Validate tagged S4 output and reconstruct global air-side electric fields."""
import numpy as np


def parse(stdout):
    records = {}
    for line in stdout.splitlines():
        if not line.startswith('WS4\t'):
            continue
        _, tag, *values = line.split('\t')
        row = np.asarray([float(v) for v in values])
        if not np.isfinite(row).all():
            raise ValueError('Nonfinite S4 raw output')
        records.setdefault(tag, []).append(row)
    for tag, width in (('COUNT', 1), ('TOP', 4), ('BOTTOM', 4), ('DONE', 1)):
        if len(records.get(tag, [])) != 1 or records[tag][0].size != width:
            raise ValueError('Missing or malformed S4 record: '+tag)
    n = int(records['COUNT'][0][0])
    modes = np.asarray(records.get('MODE', []))
    if n < 1 or n != records['COUNT'][0][0] or records['DONE'][0][0] != 1 or \
            modes.shape != (n, 15) or not np.array_equal(modes[:, 0], np.arange(1, n+1)):
        raise ValueError('Incomplete S4 modal output')
    if not np.array_equal(modes[:, 1:3], modes[:, 1:3].astype(int)):
        raise ValueError('Noninteger reciprocal vectors')
    if len(np.unique(modes[:, 1:3], axis=0)) != n:
        raise ValueError('Duplicated S4 reciprocal vectors')
    raw = dict(raw_orders=modes[:, 1:3].astype(int),
        raw_forward_hxy=modes[:, [3, 5]]+1j*modes[:, [4, 6]],
        raw_backward_hxy=modes[:, [7, 9]]+1j*modes[:, [8, 10]],
        raw_flux_by_order=modes[:, 11:15], raw_top_flux=records['TOP'][0],
        raw_bottom_flux=records['BOTTOM'][0])
    if 'FIELD' in records:
        fields = np.asarray(records['FIELD'])
        if fields.ndim != 2 or fields.shape[1] != 15:
            raise ValueError('Malformed near-field samples')
        raw['near_points_um'] = fields[:, :3]
        raw['raw_near_eh_s4'] = fields[:, 3:9]+1j*fields[:, 9:15]
    return raw


def reconstruct(raw, scene, one_dimensional=False):
    orders = raw['raw_orders'].copy()
    if one_dimensional:
        orders = np.column_stack((np.zeros(len(orders), int), orders[:, 0]*scene.periods_y))
    k0 = 2*np.pi/scene.wavelength_um
    incidence = np.deg2rad(scene.incidence_deg)
    kx = 2*np.pi*orders[:, 0]/scene.lx_um
    ky = k0*np.sin(incidence)+2*np.pi*orders[:, 1]/scene.ly_um
    q = np.sqrt(k0*k0-kx*kx-ky*ky+0j)
    if np.any(abs(q) < 1e-12*k0):
        raise ValueError('Exactly grazing diffraction order; shift test wavelength away from anomaly')
    hxy = raw['raw_backward_hxy']
    if one_dimensional:
        # Native 1D S4 varies along x': x'=global y, y'=-global x.
        hxy = np.column_stack((-hxy[:, 1], hxy[:, 0]))
    h = np.column_stack((hxy, (kx*hxy[:, 0]+ky*hxy[:, 1])/q))
    electric = -np.cross(np.column_stack((kx, ky, -q)), h)/k0
    electric[:, 2] *= -1
    electric *= np.exp(-1j*(k0*np.cos(incidence)+q)*scene.cap_top_um)[:, None]
    valid = (abs(q.imag) < 1e-8*k0) & (q.real > 1e-9*k0)
    power = np.zeros(len(orders))
    power[valid] = q[valid].real/k0/np.cos(incidence)*np.sum(abs(electric[valid])**2, axis=1)
    directions = np.full((len(orders), 3), np.nan)
    directions[valid] = np.column_stack((kx[valid], ky[valid], q[valid].real))/k0
    incident = raw['raw_top_flux'][0]
    if incident <= 0:
        raise ValueError('S4 incident flux must be positive')
    library_power = -raw['raw_flux_by_order'][:, 1]/incident
    error = float(np.max(abs(library_power-power)))
    if error > 1e-8 or not np.isfinite(electric).all():
        raise ValueError('S4 modal E/flux mismatch: '+str(error))
    result = dict(raw, orders=orders, electric_reflected_z0=electric,
        reflectance_by_order=power, directions=directions, propagating=valid)
    if 'raw_near_eh_s4' in raw:
        eh = raw['raw_near_eh_s4'].copy()
        if one_dimensional:
            eh[:, [0, 1]] = np.column_stack((-eh[:, 1], eh[:, 0]))
            eh[:, [3, 4]] = np.column_stack((-eh[:, 4], eh[:, 3]))
        eh[:, 2] *= -1
        eh[:, 3:5] *= -1
        eh *= np.exp(-1j*k0*np.cos(incidence)*scene.cap_top_um)
        result.update(near_electric=eh[:, :3], near_magnetic_impedance_scaled=eh[:, 3:])
    return result, dict(reflectance=float(power.sum()),
        flux_into_substrate=float(raw['raw_bottom_flux'][0]/incident),
        field_vs_library_power_absolute_difference=error)
