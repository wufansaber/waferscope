"""Exact sums of periodic propagating orders, never relabelled as continuous DSC."""
import numpy as np

from rcwa_prototype.comparison import compare_pair, modal_power


def collect(reference, sample, config):
    s = config.scene
    # Reuse the existing input-grid and modal-power checks.
    compare_pair(reference, sample, s, na=0.)
    masks = []
    for na, tilt in zip(config.numerical_apertures, config.tilts_deg):
        t = np.deg2rad(tilt)
        axis = np.array([np.sin(t), 0., np.cos(t)])
        masks.append(sample['propagating'] & (na > 0) &
                     (sample['directions']@axis >= np.sqrt(1-na*na)))
    masks = np.array(masks)
    delta = sample['electric_reflected_z0']-reference['electric_reflected_z0']
    efficiency = {
        'reference': reference['reflectance_by_order'], 'sample': sample['reflectance_by_order'],
        'delta': modal_power(delta, sample['directions'], sample['propagating'], s.incidence_deg)}
    efficiency['difference'] = efficiency['sample']-efficiency['reference']
    efficiency['interference'] = efficiency['difference']-efficiency['delta']
    arrays = dict(orders=sample['orders'], directions=sample['directions'],
        propagating=sample['propagating'], channel_masks=masks, union_mask=np.any(masks, axis=0),
        reflected_electric_reference=reference['electric_reflected_z0'],
        reflected_electric_sample=sample['electric_reflected_z0'], reflected_electric_delta=delta)
    metrics = []
    for i, (mask, na, tilt) in enumerate(zip(masks, config.numerical_apertures, config.tilts_deg)):
        row = dict(channel=i+1, na=na, tilt_deg=tilt, orders=int(mask.sum()), enabled=na > 0)
        for name, values in efficiency.items():
            row[name+'_power_w'] = float(values[mask].sum()*config.incident_cell_power_w)
        metrics.append(row)
    totals = {}
    for name, values in efficiency.items():
        arrays['efficiency_'+name] = values
        totals[name+'_reflected_power_w'] = float(values.sum()*config.incident_cell_power_w)
        totals[name+'_union_power_w'] = float(values[arrays['union_mask']].sum()*config.incident_cell_power_w)
    return arrays, dict(channels=metrics, totals=totals,
        incident_cell_power_w=config.incident_cell_power_w,
        meaning='Periodic diffraction orders; union counts each order once. Not isolated DSC.')
