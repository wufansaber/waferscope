"""Bracket/bisect incident coverage; this does not certify scattering convergence."""
from dataclasses import replace
import numpy as np

from .finite_illumination import power_area
from .surface_scattering import surface_fields


def coverage_at_width(p, width, spacing, taper_um):
    from .whole_surface import window
    n = max(64, 2**int(np.ceil(np.log2(width/spacing))))
    taper = max(taper_um, 2*width/n)
    axis = (np.arange(96)-48)*width/96
    previous, stable = None, 0
    general = p.beam_definition in ('wafer_rectangle', 'measured_wafer')
    for order in ((31, 41, 61, 81, 101, 151, 251, 401, 601, 1001, 1601) if general else (31, 41, 61, 81, 101, 151)):
        candidate = replace(p, beam_order=max(p.beam_order, 61) if general else order)
        if general:
            from .surface_beam import surface_fields_order
            total, ei, hi = surface_fields_order(axis, axis, candidate, order)
        else:
            total, ei, hi = surface_fields(axis, axis, candidate)
        if previous is not None:
            change = float(np.linalg.norm(total-previous)/max(np.linalg.norm(total), 1e-300))
            stable = stable+1 if change < .001 else 0
            if stable >= (2 if general else 1):
                break
        previous = total
    else:
        raise ValueError('Illumination quadrature did not converge over the whole surface')
    if general:
        # Pixel envelopes can contain counter-propagating beat fringes; 96 points
        # can alias their power even when the angular integral is stable.
        resolved = (np.arange(n)-n//2)*width/n
        incident_sum = 0.
        for start in range(0, n, 256):
            rows = resolved[start:start+256]
            _, incident, magnetic = surface_fields_order(resolved, rows, candidate, order)
            flux = -np.real(np.cross(incident, magnetic.conj()))[..., 2]
            incident_sum += float(np.sum(flux*window(resolved, rows, width, taper)**2))
        fraction = incident_sum*(width/n)**2/power_area(candidate)
    else:
        flux = -np.real(np.cross(ei, hi.conj()))[..., 2]
        fraction = float(np.sum(flux*window(axis, axis, width, taper)**2)*(width/96)**2/power_area(candidate))
    if not np.isfinite(fraction) or not 0 <= fraction <= 1.01:
        raise ValueError('Coverage quadrature is unresolved; cannot certify a whole-beam surface')
    return candidate, dict(width_um=width, grid_points=n, spacing_um=width/n,
        taper_width_um=taper, planning_power_coverage=fraction,
        field_quadrature_relative_l2_change=change, beam_order=candidate.beam_order,
        surface_beam_order=order, power_coverage_samples_per_axis=n if general else 96)


def coverage_plan(p, spacing, target, taper_um, max_points):
    from .beam_inputs import wafer_dimensions
    low = max(wafer_dimensions(p))
    high = 3*low
    trials = []

    def evaluate(width):
        n = max(64, 2**int(np.ceil(np.log2(width/spacing))))
        if n > max_points:
            raise ValueError('Whole-beam surface needs {} nodes/axis (limit {}). '
                'Do not coarsen the PSD silently; increase resources or revise the physical bandwidth.'.format(n, max_points))
        candidate, result = coverage_at_width(p, width, spacing, taper_um)
        trials.append(dict(width_um=width, coverage=result['planning_power_coverage'], beam_order=candidate.beam_order))
        return candidate, result

    _, bottom = evaluate(low)
    if bottom['planning_power_coverage'] >= target:
        low = 0
    for _ in range(8):
        candidate, result = evaluate(high)
        if result['planning_power_coverage'] >= target:
            break
        low, high = high, 2*high
    else:
        raise ValueError('Cannot bracket requested illumination coverage')
    for _ in range(12):
        if high-low <= .01*high:
            break
        mid = (low+high)/2
        new_candidate, new_result = evaluate(mid)
        if new_result['planning_power_coverage'] >= target:
            high, candidate, result = mid, new_candidate, new_result
        else:
            low = mid
    result.update(selection_method='Bracket then bisect to 1% width tolerance; illumination coverage only',
                  coverage_bracket_um=[low, high], planning_trials=trials, scattering_convergence_verified=False)
    return candidate, result
