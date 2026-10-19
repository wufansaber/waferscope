"""Same-realization outer-SPM range, spatial quadrature and FFT diagnostics."""
from dataclasses import replace

import numpy as np
from scipy.signal import resample

from .finite_illumination import power_area
from .shared_surface import sample_mother, load_mother
from .surface_scattering import SurfaceRadiation, surface_fields, direct_amplitude
from .whole_surface import window
from .range_observables import sample_caps, cap_powers, cap_field_changes


def refined_mother(axis, height):
    """Fourier-refine the stored band-limited realization, adding no random modes."""
    n = len(axis)
    result = resample(resample(height, 2*n, axis=0), 2*n, axis=1).real
    spacing = (axis[1]-axis[0])/2
    return (np.arange(2*n)-n)*spacing, result


def _observables(radiation, apertures, tilts, order=64):
    return cap_powers(sample_caps(radiation, apertures, tilts, order))


def _relative_change(current, reference):
    return (abs(current-reference)/np.maximum(abs(reference), 1e-300)).tolist()


def checked_audit_beam(p, width):
    axis = (np.arange(64)-32)*width/64
    previous = surface_fields(axis, axis, p)[0]
    checks = []
    for order in (41, 61, 81, 101, 151):
        if order <= p.beam_order:
            continue
        candidate = replace(p, beam_order=order)
        fields = surface_fields(axis, axis, candidate)[0]
        change = float(np.linalg.norm(fields-previous)/max(np.linalg.norm(fields), 1e-300))
        checks.append(dict(orders=[p.beam_order, order], relative_l2_change=change))
        p, previous = candidate, fields
        if change <= .001:
            return p, checks, True
    return p, checks, False


def surface_variant(mother_axis, mother_height, n, spacing, taper_um, p, padding, fields_provider=None):
    grid = (np.arange(n)-n//2)*spacing
    height = sample_mother(mother_axis, mother_height, grid, grid)
    taper = window(grid, grid, n*spacing, taper_um)
    fields_provider = fields_provider or surface_fields
    weighted = np.empty((n, n, 3), complex)
    incident_sum = 0.
    # Only weighted E is retained; chunk incident E/H to bound wide-field memory.
    for start in range(0, n, 256):
        stop = min(start+256, n)
        fields, ei, hi = fields_provider(grid, grid[start:stop], p)
        weighted[start:stop] = height[start:stop, :, None]*taper[start:stop, :, None]*fields
        flux = -np.real(np.cross(ei, hi.conj()))[..., 2]
        incident_sum += float(np.sum(flux*taper[start:stop]**2))
    del fields, ei, hi, flux
    coverage = incident_sum*spacing**2/power_area(p)
    return SurfaceRadiation(grid, grid, weighted, p, padding=padding), coverage, grid, weighted


def audit_padding(n, requested=4, maximum_fft_axis=4096):
    if n*requested <= maximum_fft_axis:
        return requested
    if 2*n <= maximum_fft_axis:
        return 2
    raise ValueError('Surface audit exceeds the bounded FFT workspace; reduce the field or use a larger-memory audit')


def audit_surface(source, data, metadata, settings):
    """Not a FEM mesh, PSD-tail, ensemble or optical pupil convergence certificate."""
    from .parameters import Parameters
    p = Parameters(**metadata['parameters'])
    axis, height, mother = load_mother(source)
    fine_axis, fine_height = refined_mother(axis, height)
    plan = metadata['plan']
    n0, step = plan['grid_points'], plan['spacing_um']
    provider = None
    if p.beam_definition in ('wafer_rectangle', 'measured_wafer'):
        from .surface_beam import checked_surface_order, surface_fields_order
        audit_order, beam_checks, beam_met = checked_surface_order(p, 2*round(1.875*n0/2)*step)
        provider = lambda x, y, parameters: surface_fields_order(x, y, parameters, audit_order)
    else:
        p, beam_checks, beam_met = checked_audit_beam(p, 2*round(1.875*n0/2)*step)
        audit_order = p.beam_order
    probes = np.array([[0, 0, 1], [.3, .4, np.sqrt(.75)],
        [-.4, .1, np.sqrt(.83)], [0, -.8, .6], [.65, -.3, np.sqrt(1-.65**2-.3**2)]])
    base_padding = audit_padding(n0)
    base, _, grid, weighted = surface_variant(axis, height, n0, step, plan['taper_width_um'], p, base_padding,
                                             fields_provider=provider)
    source_beam_change = float(np.linalg.norm(weighted-data['whole_weighted_fields'])/
                               max(np.linalg.norm(weighted), 1e-300))
    baseline_fields = base.amplitude(probes)
    baseline_caps = sample_caps(base, settings.numerical_apertures, settings.tilt_degrees)
    baseline_power = cap_powers(baseline_caps)
    angular_check = _observables(base, settings.numerical_apertures, settings.tilt_degrees, order=96)
    direct = direct_amplitude(grid, grid, weighted, probes, p)
    direct_error = float(np.linalg.norm(baseline_fields-direct)/max(np.linalg.norm(direct), 1e-300))
    del base, grid, weighted
    variants = []
    range_samples = []
    for label, n, spacing, padding in (
            ('smaller_range', 2*round(.8*n0/2), step, 4),
            ('larger_range', 2*round(1.25*n0/2), step, 4),
            ('extended_range', 2*round(1.5625*n0/2), step, 4),
            ('largest_shared_range', 2*round(1.875*n0/2), step, 4),
            ('finer_spatial_quadrature', 2*n0, step/2, 4),
            ('less_fft_padding', n0, step, 2)):
        chosen_axis, chosen_height = (fine_axis, fine_height) if spacing < step else (axis, height)
        requested_padding = padding
        padding = audit_padding(n, padding)
        radiation, coverage, variant_grid, variant_weighted = surface_variant(chosen_axis, chosen_height, n, spacing,
                                                    plan['taper_width_um'], p, padding, fields_provider=provider)
        variant_direct = direct_amplitude(variant_grid, variant_grid, variant_weighted, probes, p)
        variant_error = float(np.linalg.norm(radiation.amplitude(probes)-variant_direct)/
                              max(np.linalg.norm(variant_direct), 1e-300))
        del variant_grid, variant_weighted
        fields_change = float(np.linalg.norm(radiation.amplitude(probes)-baseline_fields)/
                              max(np.linalg.norm(baseline_fields), 1e-300))
        caps = sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees)
        power = cap_powers(caps)
        if label in ('larger_range', 'extended_range', 'largest_shared_range'):
            range_samples.append((n*spacing, radiation.amplitude(probes), power, caps))
        variants.append(dict(test=label, width_um=n*spacing, grid_points=n, spacing_um=spacing,
            fft_padding=padding, requested_fft_padding=requested_padding,
            selected_fft_vs_direct_complex_relative_l2=variant_error, taper_width_um=plan['taper_width_um'],
            input_power_coverage=coverage,
            selected_complex_field_relative_l2_change=fields_change,
            collected_spm_power_w=(power*settings.irradiance_w_um2).tolist(),
            collected_spm_complex_relative_l2_changes=cap_field_changes(caps, baseline_caps),
            collected_spm_relative_changes=_relative_change(power, baseline_power)))
        del radiation
    large = next(v for v in variants if v['test'] == 'larger_range')
    fine = next(v for v in variants if v['test'] == 'finer_spatial_quadrature')
    numerical_met = all(v <= .01 for item in (large, fine)
                        for v in item['collected_spm_relative_changes'])
    numerical_met &= all(item['selected_complex_field_relative_l2_change'] <= .01 for item in (large, fine))
    numerical_met &= all(v <= .01 for item in (large, fine)
                         for v in item['collected_spm_complex_relative_l2_changes'])
    angular_changes = _relative_change(angular_check, baseline_power)
    adjacent = []
    previous = (plan['width_um'], baseline_fields, baseline_power, baseline_caps)
    for current in range_samples:
        adjacent.append(dict(widths_um=[previous[0], current[0]],
            complex_field_relative_l2_change=float(np.linalg.norm(current[1]-previous[1])/
                                                   max(np.linalg.norm(current[1]), 1e-300)),
            cap_complex_relative_l2_changes=cap_field_changes(current[3], previous[3]),
            collected_spm_relative_changes=_relative_change(current[2], previous[2])))
        previous = current
    return dict(model='same_realization_outer_spm_audit', mother_sha256=mother['sha256'],
        baseline=dict(width_um=plan['width_um'], grid_points=n0, spacing_um=step,
            fft_padding=base_padding,
            collected_spm_power_w=(baseline_power*settings.irradiance_w_um2).tolist()),
        variants=variants, angular_orders=[64, 96], angular_relative_changes=angular_changes,
        adjacent_range_checks=adjacent,
        audit_beam_order=audit_order, beam_quadrature_checks=beam_checks,
        beam_quadrature_target_met=beam_met, source_vs_audit_beam_weighted_field_relative_l2=source_beam_change,
        selected_fft_vs_direct_complex_relative_l2=direct_error, numerical_target=.01,
        outer_spm_numerical_target_met=bool(numerical_met and direct_error <= .01 and
                                           all(v['selected_fft_vs_direct_complex_relative_l2'] <= .01 for v in variants) and
                                           all(v <= .01 for v in angular_changes) and
                                           beam_met and source_beam_change <= .001),
        full_model_converged=False,
        memory_policy='Sequential variants; chunked E/H; audit FFT axis <=4096. Padding can be 2 instead '
                      'of requested 4, explicitly recorded and checked against direct Fourier quadrature. '
                      'Physical surface extent, height samples and PSD bandwidth are not reduced.',
        scope='Same stored band-limited surface, same physical taper width and coherent beam. '
              'Finer quadrature uses Fourier interpolation with NO new random modes or PSD tail. '
              'The same physical beam is integrated at an independently checked higher ray order '
              'over the largest audit window; source illumination is also compared at its own window. '
              'Only outer SPM is checked; local FEM correction is fixed. '
              'Not FEM mesh/PML, PSD bandwidth, ensemble, pupil or experimental validation.')
