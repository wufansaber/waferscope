"""Replayable finite-grating fields and independent near-to-far diagnostics."""
import numpy as np

from .patterned_wafer import mirror_sources
from .radiation import radiation_amplitude
from .reciprocity_surface import box_quadrature, amplitude, evaluate_half_field


def relative(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Matching finite arrays required')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-300))


def cap_checks(mesh, field, p, g, source_order=9, surface_orders=(48, 96)):
    from .direct_imaging import polar_cap
    from .volume_sources import polarization_sources
    _, area, directions = polar_cap(.6, 0, 4, 12)
    data = dict(cap_directions=directions, cap_solid_angle_weights=area/directions[:, 2])
    checks = dict(surface_checks=[], numerical_status='mesh_and_PML_convergence_pending')
    previous = None
    for order in (source_order, source_order+2):
        print('Volume radiation quadrature {}'.format(order), flush=True)
        points, weights, volumes = polarization_sources(mesh, field, p, False, order)
        points, weights = mirror_sources(points, weights)
        result = radiation_amplitude(points, weights, directions, p)
        data['volume_q{}'.format(order)] = result
        if previous is not None:
            checks['volume_quadrature_relative_l2'] = relative(previous, result)
        previous = result
    checks['trench_volume_um3'] = 2*volumes['pit']
    checks['CH3_NA_0p6_DSC_um2'] = float(np.sum(data['cap_solid_angle_weights']*
        np.sum(abs(result)**2, axis=1)))
    checks['cap_quadrature'] = dict(radial_order=4, azimuth_count=12,
        scope='Selected-ray extraction diagnostic; collection-power angular convergence not certified')
    for fraction in (.35, .65):
        previous = None
        for order in surface_orders:
            print('Reciprocity box {}, quadrature {}'.format(fraction, order), flush=True)
            surface, normals, areas = box_quadrature(g.surface_bounds(fraction), order)
            electric, magnetic = evaluate_half_field(mesh, field, surface, p)
            surface_result = amplitude(surface, normals, areas, electric, magnetic, directions, p)
            checks['surface_checks'].append(dict(clearance_fraction=fraction, order=order,
                volume_vs_surface_relative_l2=relative(surface_result, result),
                quadrature_change_relative_l2=relative(previous, surface_result) if previous is not None else None))
            data['surface_f{}_q{}'.format(fraction, order)] = surface_result
            previous = surface_result
    a, b = [data['surface_f{}_q{}'.format(f, surface_orders[-1])] for f in (.35, .65)]
    checks['two_surface_relative_l2'] = relative(a, b)
    return data, checks, points, weights


def near_planes(mesh, field, p, g, samples=100):
    from .planar import background_fields
    data = {}
    lateral = np.linspace(-g.half_width_um*.95, g.half_width_um*.95, samples)
    vertical = np.linspace(-g.substrate_depth_um*.9, g.air_height_um*.9, samples)
    configurations = [('xy', lateral, lateral, .02),
        ('xz', lateral, vertical, float(g.centers[0])), ('yz', lateral, vertical, 0.)]
    for name, u, v, offset in configurations:
        uu, vv = np.meshgrid(u, v)
        positions = np.zeros((uu.size, 3))
        positions[:, 'xyz'.index(name[0])] = uu.ravel()
        positions[:, 'xyz'.index(name[1])] = vv.ravel()
        fixed = next(i for i, axis in enumerate('xyz') if axis not in name)
        positions[:, fixed] = offset
        # Pick the upper-side trace at the planar interface, never average E_n across it.
        positions[np.abs(positions[:, 2]) < 1e-12, 2] = 1e-10
        esc, hsc = evaluate_half_field(mesh, field, positions, p)
        ebg, hbg = background_fields(positions, p)
        prefix = 'near_'+name+'_'
        data[prefix+'u'], data[prefix+'v'] = u, v
        data[prefix+'offset_um'] = np.asarray(offset)
        data[prefix+'electric'] = (esc+ebg).reshape(samples, samples, 3)
        data[prefix+'magnetic_z0'] = (hsc+hbg).reshape(samples, samples, 3)
        data[prefix+'intensity'] = np.sum(abs(esc+ebg)**2, axis=1).reshape(samples, samples)
    return data


def far_disk(points, weights, p, samples=48):
    axis = np.linspace(-1, 1, samples)
    xx, yy = np.meshgrid(axis, axis)
    valid = xx**2+yy**2 < 1-1e-8
    directions = np.column_stack((xx[valid], yy[valid], np.sqrt(1-xx[valid]**2-yy[valid]**2)))
    far = radiation_amplitude(points, weights, directions, p)
    disk = np.full_like(xx, np.nan)
    disk[valid] = np.sum(abs(far)**2, axis=1)
    return dict(disk_axis=axis, disk_dsc=disk, disk_directions=directions, disk_amplitude=far)
