"""Independent beam/angular refinement of the final two saved range-audit windows."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from .hybrid_run import load_sources, optical_settings
from .parameters import Parameters
from .shared_surface import recover_synthesis_parent
from .surface_beam import checked_surface_order, surface_fields_order
from .hybrid_surface_audit import surface_variant
from .range_observables import sample_caps, cap_powers, cap_field_changes


def run(previous, output):
    previous, output = Path(previous), Path(output)
    old = json.loads((previous/'range_report.json').read_text(encoding='utf-8'))
    if old['status'] != 'complete_with_physical_limits' or len(old['cases']) < 2:
        raise ValueError('A completed multi-window range audit is required')
    source = Path(old['source'])
    source_data, metadata = load_sources(source)
    del source_data
    axis, height, parent = recover_synthesis_parent(source)
    if parent['parent_sha256'] != old['parent']['parent_sha256']:
        raise ValueError('Range audit parent changed')
    p = Parameters(**old['scattering_parameters'])
    settings = optical_settings(old['imaging_settings'], p)
    order, checks, beam_met = checked_surface_order(p, old['cases'][-1]['width_um'])
    output.mkdir(parents=True, exist_ok=True)
    samples = []
    for case in old['cases'][-2:]:
        checkpoint = output/('window_{}_caps.npz'.format(case['grid_points']))
        signature = hashlib.sha256(json.dumps(dict(parent=parent['parent_sha256'],
            case=case, settings=old['imaging_settings'], parameters=old['scattering_parameters'],
            beam_order=order), sort_keys=True).encode('utf-8')).hexdigest()
        if checkpoint.is_file():
            with np.load(checkpoint, allow_pickle=False) as saved:
                if str(saved['signature']) != signature:
                    raise ValueError('Range refinement checkpoint inputs changed')
                samples.append({n: [(saved['a{}_{}'.format(n, c)], saved['w{}_{}'.format(n, c)])
                                     for c in range(3)] for n in (96, 144, 192)})
            if any(not np.isfinite(a).all() or not np.isfinite(w).all()
                   for group in samples[-1].values() for a, w in group):
                raise ValueError('Nonfinite range refinement checkpoint')
            continue
        print('Refining final range window {} um...'.format(case['width_um']), flush=True)
        radiation, coverage, grid, weighted = surface_variant(axis, height, case['grid_points'], parent['spacing_um'],
            metadata['plan']['taper_width_um'], p, case['fft_padding'],
            fields_provider=lambda x, y, parameters: surface_fields_order(x, y, parameters, order))
        samples.append({n: sample_caps(radiation, settings.numerical_apertures, settings.tilt_degrees, n)
                        for n in (96, 144, 192)})
        arrays = {'signature': signature}
        for n, group in samples[-1].items():
            for c, (a, w) in enumerate(group):
                arrays['a{}_{}'.format(n, c)], arrays['w{}_{}'.format(n, c)] = a, w
        np.savez_compressed(checkpoint, **arrays)
        del radiation, weighted, grid
    changes = {n: cap_field_changes(samples[1][n], samples[0][n]) for n in (96, 144, 192)}
    power_changes = abs(cap_powers(samples[1][192])-cap_powers(samples[1][144]))/np.maximum(
        abs(cap_powers(samples[1][192])), 1e-300)
    complex_changes = abs(np.asarray(changes[192])-changes[144])
    passed = bool(beam_met and max(changes[192]) <= .01 and max(power_changes) <= .01 and
                  max(complex_changes) <= .001 and old['cases'][0][
                  'source_vs_refined_beam_weighted_field_relative_l2'] <= .001 and old['cases'][-1][
                  'fft_vs_direct_probe_complex_relative_l2'] <= .001)
    old.update(previous_audit=str(previous.resolve()), audit_beam_order=order,
        beam_quadrature_checks=checks, beam_quadrature_target_met=beam_met,
        last_adjacent_complex_target_met=passed,
        independent_refinement=dict(angular_orders=[96, 144, 192],
            last_adjacent_complex_changes=changes,
            angular_power_relative_changes_144_to_192=power_changes.tolist(),
            angular_complex_error_absolute_changes_144_to_192=complex_changes.tolist(),
            collected_spm_power_w_192=(cap_powers(samples[1][192])*settings.irradiance_w_um2).tolist()),
        production_window_promoted=False)
    (output/'range_report.json').write_text(json.dumps(old, indent=2), encoding='utf-8')
    from .validation_plotting import plot_range
    plot_range(output, old)
    return old['independent_refinement'], passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.previous, args.output), indent=2), flush=True)


if __name__ == '__main__':
    main()
