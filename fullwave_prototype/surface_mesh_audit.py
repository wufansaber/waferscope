"""Same-parent spatial quadrature and edge-taper checks at the largest verified range."""
import argparse
import json
from pathlib import Path
import numpy as np
from .hybrid_run import load_sources, optical_settings
from .hybrid_surface_audit import surface_variant, refined_mother
from .parameters import Parameters
from .shared_surface import recover_synthesis_parent
from .surface_beam import checked_surface_order, surface_fields_order
from .range_observables import sample_caps, cap_field_changes, cap_powers


def run(source, output):
    _, metadata = load_sources(source)
    p = Parameters(**metadata['parameters'])
    inputs = json.loads((Path(source)/'hybrid_input.json').read_text(encoding='utf-8'))
    settings = optical_settings(inputs['imaging_settings'], p)
    axis, height, parent = recover_synthesis_parent(source)
    fine_axis, fine_height = refined_mother(axis, height)
    plan = metadata['plan']
    n, step = 2*round(3.75*plan['grid_points']/2), plan['spacing_um']
    taper = plan['taper_width_um']
    order, checks, met = checked_surface_order(p, n*step)
    provider = lambda x, y, parameters: surface_fields_order(x, y, parameters, order)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = dict(model='same_parent_spatial_taper_audit', status='running', parent=parent,
        width_um=n*step, beam_order=order, beam_checks=checks, beam_target_met=met, cases=[],
        full_model_converged=False, scope='Spatial integration of the same band-limited realization; '
            'no new high-frequency PSD modes, no new realization, outer SPM only. Taper changes affect physical truncation.')
    baseline = None
    for name, count, spacing, edge, ax, h in (
        ('baseline', n, step, taper, axis, height),
        ('coarser_grid', n//2, 2*step, taper, axis, height),
        ('narrower_taper', n, step, taper/2, axis, height),
        ('wider_taper', n, step, 2*taper, axis, height),
        ('finer_same_band_grid', 2*n, step/2, taper, fine_axis, fine_height)):
        print('Spatial/taper audit: '+name, flush=True)
        radiation, coverage, _, _ = surface_variant(ax, h, count, spacing, edge, p, 4, provider)
        samples = sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees, 96)
        powers = cap_powers(samples)
        case = dict(name=name, spacing_um=spacing, taper_width_um=edge, coverage=coverage,
                    cap_power_w=(powers*settings.irradiance_w_um2).tolist())
        if baseline is None:
            baseline = samples
        else:
            case['complex_relative_l2'] = cap_field_changes(samples, baseline)
            case['power_relative_changes'] = (abs(powers-cap_powers(baseline))/np.maximum(powers, 1e-300)).tolist()
        report['cases'].append(case)
        (output/'surface_mesh_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        del radiation
    relevant = [c for c in report['cases'] if c['name'] not in ('baseline', 'coarser_grid')]
    report.update(status='complete_with_physical_limits', tested_spatial_taper_target_met=bool(met and all(
        max(c['complex_relative_l2']+c['power_relative_changes']) <= .01 for c in relevant)))
    (output/'surface_mesh_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    from .validation_plotting import plot_mesh
    plot_mesh(output, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), indent=2), flush=True)


if __name__ == '__main__':
    main()
