"""Collect completed image checks, sampling sensitivity and measured resource use."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .parameters import Parameters


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def summarize(output, preview, fine, previous):
    result = dict(scope='Selected benchmark cases and order-2 rough pipeline only; not universal accuracy.',
                  flat=[], sampling=[], solver_comparison={}, resources={})
    for pol in ('p', 's'):
        report = read(output/('flat_'+pol)/'report.json')
        result['flat'].append(dict(polarization=pol, channels=report['channels']))
    base, refined = read(preview/'imaging_report.json'), read(fine/'imaging_report.json')
    result['scattering_parameters'] = base['scattering_parameters']
    result['imaging_settings'] = base['settings']
    result['roughness'] = base['roughness']
    for before, after in zip(base['channels'], refined['channels']):
        low, high = before['with']['image_power_w'], after['with']['image_power_w']
        result['sampling'].append(dict(tilt_deg=before['tilt_deg'], na=before['na'],
            power64_pw=low*1e12, power128_pw=high*1e12,
            relative_change_to_128=abs(low-high)/max(abs(high), 1e-300),
            edge_fraction64=before['fft_edge_abs_flux_fraction'],
            edge_fraction128=after['fft_edge_abs_flux_fraction']))
    old_report = read(previous/'scattering/report.json')
    new_report = read(preview/'scattering/report.json')
    old_parameters, new_parameters = (asdict(Parameters(**r['parameters'])) for r in (old_report, new_report))
    old_parameters.pop('direct_solver')
    new_parameters.pop('direct_solver')
    if old_parameters != new_parameters or old_report['roughness'] != new_report['roughness']:
        raise ValueError('Solver comparison requires identical physical and discretization parameters')
    with np.load(previous/'scattering/paired_fields.npz') as old, np.load(preview/'scattering/paired_fields.npz') as new:
        for key in ('near_y_um', 'near_z_um', 'theta', 'phi', 'surface_vertices_um'):
            if not np.array_equal(old[key], new[key]):
                raise ValueError('Solver comparison grids differ: '+key)
        for key in ('near_e_with', 'near_e_without', 'amplitude_with', 'amplitude_without'):
            result['solver_comparison'][key+'_relative_l2'] = float(
                np.linalg.norm(new[key]-old[key])/max(np.linalg.norm(old[key]), 1e-300))
    result['solver_comparison']['previous_report'] = str((previous/'scattering/report.json').resolve())
    result['solver_comparison']['current_report'] = str((preview/'scattering/report.json').resolve())
    for name in ('flat', 'pipeline', 'pipeline_umfpack', 'pipeline_preview', 'pupil128'):
        result['resources'][name] = read(output/('resources_'+name+'.json'))
    (output/'followup_summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/followup_validation')
    parser.add_argument('--preview', default='build/fullwave_imaging/complete_preview')
    parser.add_argument('--fine', default='build/fullwave_imaging/complete_pupil128')
    parser.add_argument('--previous', default='build/fullwave_imaging/default')
    args = parser.parse_args()
    summarize(Path(args.output), Path(args.preview), Path(args.fine), Path(args.previous))


if __name__ == '__main__':
    main()
