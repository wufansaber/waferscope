"""Unfitted flux, transverse-field and mirror checks on saved periodic orders."""
import numpy as np


def scene_diagnostics(far, native, config):
    keep = far['propagating']
    rays = far['directions'][keep]
    weights = rays[:, 2]/np.cos(np.deg2rad(config.scene.incidence_deg))
    orders = far['orders'][keep]
    lookup = {tuple(g): i for i, g in enumerate(orders)}
    mirrors = np.array([lookup[(-gx, gy)] for gx, gy in orders])
    parity = np.array([-1., 1., 1.])*(1 if config.polarization == 'p' else -1)
    checks = {}
    for label in ('reference', 'sample', 'delta'):
        field = far['reflected_electric_'+label][keep]
        norm = float(np.sum(weights[:, None]*abs(field)**2))
        transverse = float(np.sum(weights*abs(np.sum(rays*field, axis=1))**2))
        mirror = float(np.sum(weights[:, None]*abs(field-field[mirrors]*parity)**2))
        checks[label] = dict(
            transverse_relative_l2=float(np.sqrt(transverse/norm)) if norm > 1e-28 else None,
            mirror_relative_l2=float(np.sqrt(mirror/norm)) if norm > 1e-28 else None,
            normalization='Flux-weighted norm of this field; no phase or amplitude fitting.')
        if label != 'delta':
            row = native[label]
            reflection, transmission = row['reflectance'], row['flux_into_substrate']
            balance = 1-reflection-transmission
            checks[label].update(reflectance=reflection, flux_into_substrate=transmission,
                finite_stack_absorption_by_balance=balance,
                passive_flux_check_passed=bool(np.isfinite([reflection, transmission, balance]).all()
                    and min(reflection, transmission, balance) >= -1e-8),
                modal_e_vs_native_flux_max_absolute_error=row['field_vs_library_power_absolute_difference'])
    symmetric = (config.numerical_apertures[0] == config.numerical_apertures[1]
        and config.tilts_deg[0] == -config.tilts_deg[1])
    asymmetry = {}
    for name in ('reference', 'sample', 'delta'):
        value = far['efficiency_'+name]
        a, b = [float(value[m].sum()) for m in far['channel_masks'][:2]]
        asymmetry[name] = 2*abs(a-b)/(a+b) if symmetric and a+b > 1e-28 else None
    mirror_passed = all(row['mirror_relative_l2'] is None or row['mirror_relative_l2'] < .02
                        for row in checks.values())
    return dict(fields=checks, mirror_check_passed=mirror_passed, mirror_relative_tolerance=.02,
        side_channels_mirror_matched=symmetric,
        side_channel_power_asymmetry=asymmetry,
        asymmetry_definition='2*abs(P1-P2)/(P1+P2); not a peak-normalized comparison.',
        scope='Centered sphere and X-invariant line grating with kx incident=0; no forced symmetry.',
        energy_scope='1-R-T is flux bookkeeping, not independently integrated material absorption.',
        certification='These are consistency diagnostics, not cutoff, isolation or camera certification.')
