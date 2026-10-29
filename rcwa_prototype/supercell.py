"""Periodic-cell sensitivity; never subtract fields at unlike wave vectors."""
import numpy as np

from .geometry import Scene


def compare_cells(coarse, fine):
    for key in ('backend', 'backend_version', 'factorization', 'executable_sha256',
                'epsilon_method', 'stripe_rule', 'particle_method'):
        if coarse['info'].get(key) != fine['info'].get(key):
            raise ValueError('Same solver and material rules required: '+key)
    a, b = coarse['info']['scene'], fine['info']['scene']
    sa, sb = Scene(**a), Scene(**b)
    if any(a[k] != b[k] for k in a if k not in
           ('lx_um', 'periods_y', 'harmonics_x', 'harmonics_y')):
        raise ValueError('Cell comparison must preserve materials, illumination, particle and slicing')
    if (sa.lx_um, sa.ly_um) == (sb.lx_um, sb.ly_um):
        raise ValueError('Cell dimensions must differ')
    if sb.lx_um < sa.lx_um or sb.ly_um < sa.ly_um:
        raise ValueError('Cell dimensions must not decrease')
    if None in (sa.harmonics_x, sa.harmonics_y, sb.harmonics_x, sb.harmonics_y):
        raise ValueError('Explicit rectangular cutoffs required')
    fa = np.array([sa.harmonics_x/sa.lx_um, sa.harmonics_y/sa.ly_um])
    fb = np.array([sb.harmonics_x/sb.lx_um, sb.harmonics_y/sb.ly_um])
    if not np.allclose(fa, fb, rtol=1e-12, atol=1e-12):
        raise ValueError('Hold physical spatial-frequency cutoffs, not order counts, constant')
    factor_a = sa.lx_um*sa.ly_um*np.cos(np.deg2rad(sa.incidence_deg))
    factor_b = sb.lx_um*sb.ly_um*np.cos(np.deg2rad(sb.incidence_deg))
    ha = coarse['summary']['hemisphere_difference_field_efficiency']*factor_a
    hb = fine['summary']['hemisphere_difference_field_efficiency']*factor_b
    rows = []
    if len(coarse['summary']['channels']) != len(fine['summary']['channels']):
        raise ValueError('Identical channel count required')
    for x, y in zip(coarse['summary']['channels'], fine['summary']['channels']):
        if (x['na'], x['tilt_deg']) != (y['na'], y['tilt_deg']):
            raise ValueError('Identical collection apertures required')
        ca, cb = x['difference_field_efficiency']*factor_a, y['difference_field_efficiency']*factor_b
        rows.append(dict(na=x['na'], tilt_deg=x['tilt_deg'],
            coarse_collected_orders=x['collected_orders'], fine_collected_orders=y['collected_orders'],
            coarse_periodic_cross_section_um2=ca, fine_periodic_cross_section_um2=cb,
            relative_cross_section_change=cb/ca-1 if ca else None,
            coarse_actual_power_difference_cross_section_um2=x['actual_power_difference']*factor_a,
            fine_actual_power_difference_cross_section_um2=y['actual_power_difference']*factor_b))
    return dict(coarse_cell_um=[sa.lx_um, sa.ly_um], fine_cell_um=[sb.lx_um, sb.ly_um],
        physical_frequency_halfwidth_cycles_per_um=fa.tolist(), channels=rows,
        coarse_hemisphere_cross_section_um2=ha, fine_hemisphere_cross_section_um2=hb,
        relative_hemisphere_cross_section_change=hb/ha-1 if ha else None,
        status='PERIODIC_CELL_SENSITIVITY_NOT_ISOLATED_DEFECT_CONVERGENCE',
        field_difference_computed=False,
        caveat='Lattice interactions AND discrete angular/aperture sampling can change; two cells do not separate them.')
