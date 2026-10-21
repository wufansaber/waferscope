"""Compare two unchanged-geometry Maxwell meshes without fitted scaling."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse', type=Path, required=True)
    parser.add_argument('--fine', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    coarse, fine = [read(folder/'imaging_report.json') for folder in (args.coarse, args.fine)]
    assert coarse['pattern'] == fine['pattern']
    assert coarse['settings'] == fine['settings']
    with np.load(args.coarse/'image_fields.npz') as saved:
        old = {k: saved[k] for k in saved.files if '_detector_' not in k}
    with np.load(args.fine/'image_fields.npz') as saved:
        new = {k: saved[k] for k in saved.files if '_detector_' not in k}
    rows = []
    for i, (a, b) in enumerate(zip(coarse['channels'], fine['channels']), 1):
        prefix = 'ch{}_'.format(i)
        np.testing.assert_array_equal(old[prefix+'image_axis_um'], new[prefix+'image_axis_um'])
        errors = {state: float(np.linalg.norm(new[prefix+'image_'+state+'_w_um2']-
            old[prefix+'image_'+state+'_w_um2'])/np.linalg.norm(new[prefix+'image_'+state+'_w_um2']))
            for state in ('with', 'without', 'difference')}
        power = {state: (a[state]['image_power_w']-b[state]['image_power_w'])/b[state]['image_power_w']
                 for state in ('with', 'without')}
        rows.append(dict(channel=i, flux_relative_l2=errors, power_relative_error=power))
    scattering = read(args.fine/'scattering_report.json')
    quadrature = {state: scattering['metrics'][state]['radiation_quadrature_relative_l2']
                  for state in ('with', 'without')}
    passed = (max(quadrature.values()) < .01 and
              all(max(abs(v) for v in row['power_relative_error'].values()) < .05 and
                  max(row['flux_relative_l2'].values()) < .1 for row in rows))
    result = dict(status='passed_selected_checks' if passed else 'not_converged',
        channels=rows, quadrature_relative_l2=quadrature,
        mesh_dofs=[read(folder/'scattering_report.json')['metrics']['with']['dofs']
                   for folder in (args.coarse, args.fine)],
        symmetry_scope='Left/right optical equality imposed by physically valid PMC mirror boundary for this synthetic symmetric input; not an independent convergence test.',
        limitations='Two finite trenches, P1 tetrahedra, unvalidated domain/PML truncation and substrate optical constants. '
                    'Even passing selected checks would not establish agreement with a real SEM specimen.')
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'mesh_comparison.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))

    from fullwave_prototype.detection_replay import run
    from tools.plot_patterned_wafer_case import plot_case
    replay_folder = args.output/'camera_1us'
    camera = dict(camera_profile='gemi_area', gray_snr_pairs=1000, gray_snr_seed=20261016)
    replay = run(args.fine, replay_folder, exposure_ms=.001, camera_settings=camera)
    replay['mesh_validation'] = result
    replay['status'] = result['status']+'; synthetic finite patch, not a calibrated SEM prediction'
    (replay_folder/'imaging_report.json').write_text(json.dumps(replay, indent=2), encoding='utf-8')
    with np.load(args.fine/'scattering_fields.npz') as saved:
        data = {k: saved[k] for k in saved.files}
    with np.load(replay_folder/'image_fields.npz') as saved:
        arrays = {k: saved[k] for k in saved.files}
    for i, channel in enumerate(replay['channels'], 1):
        series = arrays['ch{}_detector_gray_roi_series_dn'.format(i)][:, 2]
        empirical = abs(series.mean())/series.std(ddof=1)
        if channel['detector']['gray_valid']:
            assert np.isclose(empirical, channel['detector']['gray_snr'], rtol=1e-13)
    plot_case(replay_folder, data, arrays, replay)
    destination = Path('output/patterned_wafer_validation_20261019')
    destination.mkdir(parents=True, exist_ok=True)
    for name in ('optical_comparison.png', 'camera_comparison.png', 'geometry_far_field.png',
                 'case_summary.json', 'imaging_report.json'):
        shutil.copy2(replay_folder/name, destination/name)
    shutil.copy2(args.fine/'near_field.png', destination/'near_field.png')
    shutil.copy2(args.fine/'scattering_report.json', destination/'scattering_report.json')
    shutil.copy2(args.output/'mesh_comparison.json', destination/'mesh_comparison.json')
    print('Camera 1 us replay complete; physical field and geometry unchanged.')


if __name__ == '__main__':
    main()
