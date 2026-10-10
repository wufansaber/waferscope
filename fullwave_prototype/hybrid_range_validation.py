"""Validate larger whole-beam windows on the exactly recovered original FFT parent."""
import argparse
import json
from pathlib import Path

import numpy as np

from .hybrid_run import load_sources, optical_settings
from .hybrid_surface_audit import surface_variant
from .surface_beam import checked_surface_order, surface_fields_order
from .parameters import Parameters
from .range_observables import sample_caps, cap_powers, cap_field_changes
from .shared_surface import recover_synthesis_parent
from .surface_scattering import direct_amplitude
from .local_spm_validation import probe_directions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--imaging-settings-json')
    args = parser.parse_args()
    data, metadata = load_sources(args.source)
    settings_values = (json.loads(Path(args.imaging_settings_json).read_text(encoding='utf-8'))
                      if args.imaging_settings_json else
                      json.loads((Path(args.source)/'hybrid_input.json').read_text(encoding='utf-8'))['imaging_settings'])
    p = Parameters(**metadata['parameters'])
    settings = optical_settings(settings_values, p)
    axis, height, parent = recover_synthesis_parent(args.source)
    plan = metadata['plan']
    n0, spacing = plan['grid_points'], plan['spacing_um']
    multiples = ((1., 2., 4., 8., 12., 15.) if len(axis) >= 16*n0 else
                 (1., 1.25, 1.875, 2.5, 3.125, 3.75))
    order, beam_checks, beam_met = checked_surface_order(p, 2*round(multiples[-1]*n0/2)*spacing)
    def fields_provider(x, y, parameters):
        return surface_fields_order(x, y, parameters, order)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    report = dict(model='recovered_parent_outer_spm_range_validation', status='running',
        parent=parent, source=str(Path(args.source).resolve()), audit_beam_order=order,
        scattering_parameters=metadata['parameters'], imaging_settings=settings_values,
        beam_quadrature_checks=beam_checks, beam_quadrature_target_met=beam_met, cases=[],
        full_model_converged=False,
        scope='Only outer first-order SPM on one realization. Same physical edge taper width. '
              'Caps use solid-angle-weighted complex field errors. '
              'No local FEM correction, new random draw or periodic repetition is introduced.')
    previous = baseline = previous_finest = None
    for multiple in multiples:
        n = 2*round(multiple*n0/2)
        print('Recovered-parent window {:.6g} um, {} nodes/axis...'.format(n*spacing, n), flush=True)
        padding = 2 if n > 1536 else 4
        radiation, coverage, grid, weighted = surface_variant(axis, height, n, spacing,
            plan['taper_width_um'], p, padding, fields_provider=fields_provider)
        samples = sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees, 96)
        powers = cap_powers(samples)
        fine = sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees, 144)
        finest = sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees, 192)
        case = dict(width_um=n*spacing, grid_points=n, fft_padding=padding, input_power_coverage=coverage,
                    collected_spm_power_w=(powers*settings.irradiance_w_um2).tolist())
        if baseline is None:
            baseline = samples
            case['source_vs_refined_beam_weighted_field_relative_l2'] = float(
                np.linalg.norm(weighted-data['whole_weighted_fields'])/max(np.linalg.norm(weighted), 1e-300))
        else:
            case['from_baseline_cap_complex_relative_l2'] = cap_field_changes(samples, baseline)
            case['from_previous_cap_complex_relative_l2'] = cap_field_changes(samples, previous)
            prev_power = cap_powers(previous)
            case['from_previous_collected_power_relative_changes'] = (
                abs(powers-prev_power)/np.maximum(abs(powers), 1e-300)).tolist()
        # Independently check the highest window's angular quadrature as its spectrum is sharper.
        if multiple == multiples[-1]:
            fine_power = cap_powers(finest)
            case['angular_orders'] = [96, 144, 192]
            case['angular_power_relative_changes'] = (
                abs(fine_power-cap_powers(fine))/np.maximum(abs(fine_power), 1e-300)).tolist()
            fine_changes = cap_field_changes(finest, previous_finest)
            case['fine_from_previous_cap_complex_relative_l2'] = fine_changes
            case['angular_complex_error_absolute_changes'] = (
                abs(np.asarray(fine_changes)-case['from_previous_cap_complex_relative_l2'])).tolist()
            probes = probe_directions()[::7]
            direct = direct_amplitude(grid, grid, weighted, probes, p)
            case['fft_vs_direct_probe_complex_relative_l2'] = float(
                np.linalg.norm(radiation.amplitude(probes)-direct)/max(np.linalg.norm(direct), 1e-300))
        report['cases'].append(case)
        previous = samples
        previous_finest = finest
        del radiation, weighted, grid
        (output/'range_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    last = report['cases'][-1]
    report.update(status='complete_with_physical_limits', numerical_target=.01,
        last_adjacent_complex_target_met=bool(beam_met and
            report['cases'][0]['source_vs_refined_beam_weighted_field_relative_l2'] <= .001 and
            last['fft_vs_direct_probe_complex_relative_l2'] <= .001 and all(
            v <= .01 for v in last['fine_from_previous_cap_complex_relative_l2']) and all(
            v <= .001 for v in last['angular_complex_error_absolute_changes']) and all(
            v <= .01 for v in last['angular_power_relative_changes'])))
    (output/'range_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    from .validation_plotting import plot_range
    plot_range(output, report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
