"""Compare periodic modal powers without treating them as continuous DSC."""
import numpy as np


def modal_power(electric, directions, propagating, incidence_deg):
    electric = np.asarray(electric)
    result = np.zeros(len(electric))
    valid = np.asarray(propagating, dtype=bool)
    result[valid] = (directions[valid, 2] / np.cos(np.deg2rad(incidence_deg))
                     * np.sum(abs(electric[valid])**2, axis=1))
    return result


def compare_pair(reference, sample, scene, na=.6, tilts=(45., -45., 0.)):
    if not 0 <= na <= 1:
        raise ValueError('NA must lie in [0, 1] in air')
    for key in ('orders', 'propagating'):
        if not np.array_equal(reference[key], sample[key]):
            raise ValueError('Reference and sample must use identical modal grids')
    np.testing.assert_allclose(reference['directions'], sample['directions'], equal_nan=True)
    directions = sample['directions']
    delta = sample['electric_reflected_z0']-reference['electric_reflected_z0']
    power = modal_power(delta, directions, sample['propagating'], scene.incidence_deg)
    factor = scene.lx_um*scene.ly_um*np.cos(np.deg2rad(scene.incidence_deg))
    channels = []
    for tilt in tilts:
        angle = np.deg2rad(tilt)
        axis = np.array((np.sin(angle), 0., np.cos(angle)))
        keep = sample['propagating'] & (directions@axis >= np.sqrt(1-na**2))
        if na == 0:
            keep[:] = False
        s = float(np.sum(sample['reflectance_by_order'][keep]))
        r = float(np.sum(reference['reflectance_by_order'][keep]))
        d = float(np.sum(power[keep]))
        channels.append(dict(tilt_deg=tilt, na=na, collected_orders=int(keep.sum()),
            sample_efficiency=s, reference_efficiency=r, actual_power_difference=s-r,
            difference_field_efficiency=d, interference_efficiency=s-r-d,
            periodic_difference_field_cross_section_um2=factor*d))
    return dict(channels=channels, hemisphere_difference_field_efficiency=float(power.sum()),
        scope='Periodic array, discrete diffraction orders; not isolated-particle DSC.',
        normalization='Efficiencies relative to unit plane-wave power incident on this supercell.',
        cross_section_definition='Lx*Ly*cos(theta_i)*difference_field_efficiency; periodic, not isolated.')


def field_change(coarse, fine, incidence_deg):
    def table(data):
        return {tuple(g): i for i, g in enumerate(data['orders']) if data['propagating'][i]}
    a, b = table(coarse), table(fine)
    if a.keys() != b.keys():
        return dict(comparable=False, reason='Propagating orders differ; enlarge truncation first.')
    keys = sorted(a)
    ei = coarse['electric_reflected_z0'][[a[k] for k in keys]]
    ef = fine['electric_reflected_z0'][[b[k] for k in keys]]
    direction = fine['directions'][[b[k] for k in keys]]
    weight = direction[:, 2]/np.cos(np.deg2rad(incidence_deg))
    denominator = np.sum(weight[:, None]*abs(ef)**2)
    numerator = np.sum(weight[:, None]*abs(ei-ef)**2)
    relative = float(np.sqrt(numerator/denominator)) if denominator else None
    return dict(comparable=True, common_propagating_orders=len(keys),
        relative_power_weighted_complex_field_change=relative,
        caveat='Adjacent truncation change, not an error bound or proof of convergence.')
