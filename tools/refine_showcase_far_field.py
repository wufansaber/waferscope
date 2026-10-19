"""Re-evaluate the saved hybrid source on a dense Cartesian direction disk."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'build/channel_symmetry_20261016/low_haze_02nm'
OUTPUT = ROOT/'build/far_field_figure_20261019'
ASSETS = ROOT/'docs/showcase_20261017/assets'


def source_fingerprint():
    return hashlib.sha256((SOURCE/'hybrid_source_report.json').read_bytes()).hexdigest()


def evaluate_source(data, p, directions, surface, progress=False):
    from fullwave_prototype.hybrid_radiation import triangle_amplitude
    from fullwave_prototype.radiation import radiation_amplitude
    result = surface.amplitude(directions)
    for start in range(0, len(directions), 512):
        q = directions[start:start+512]
        result[start:start+len(q)] += radiation_amplitude(data['source_points_with'],
            data['source_polarization_weights_with'], q, p, block_size=32)
        result[start:start+len(q)] -= triangle_amplitude(data['local_spm_positions_um'],
            data['local_spm_weights'], q, p)
        if progress and (start//512) % 12 == 0:
            print('Evaluated %d / %d physical directions'%(start+len(q), len(directions)), flush=True)
    return result


def worker():
    from fullwave_prototype.hybrid_run import load_sources
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.surface_scattering import SurfaceRadiation, direct_amplitude
    start = time.perf_counter()
    data, metadata = load_sources(SOURCE)
    p = Parameters(**metadata['parameters'])
    surface = SurfaceRadiation(data['whole_axis_um'], data['whole_axis_um'],
                               data['whole_weighted_fields'], p, padding=4)
    # A Cartesian disk removes polar-sector wedges and samples the finite-window spectrum.
    axis = np.linspace(-1., 1., 401)
    x, y = np.meshgrid(axis, axis)
    r2 = x*x+y*y
    mask = r2 < 1.-1e-12
    directions = np.column_stack((x[mask], y[mask], np.sqrt(1-r2[mask])))
    amplitude = evaluate_source(data, p, directions, surface, progress=True)
    dsc = np.full(x.shape, np.nan)
    dsc[mask] = np.sum(abs(amplitude)**2, axis=-1)
    with np.load(SOURCE/'hybrid_far_fields.npz') as saved:
        phi, theta = np.meshgrid(saved['phi'], saved['theta'])
        q = np.stack((np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)), axis=-1)
        selection = np.arange(q.size//3)[::11]
        coarse_q = q.reshape(-1,3)[selection]
        reference = saved['with_particle'].reshape(-1,2)[selection]
        new = evaluate_source(data, p, coarse_q, surface)
    relative = float(np.linalg.norm(new-reference)/np.linalg.norm(reference))
    assert relative < 1e-10, relative
    # Check the FFT surface transform against direct spatial sums at fixed directions.
    chosen = directions[np.linspace(0, len(directions)-1, 25).astype(int)]
    direct = direct_amplitude(data['whole_axis_um'], data['whole_axis_um'], data['whole_weighted_fields'], chosen, p)
    fft = surface.amplitude(chosen)
    fft_error = float(np.linalg.norm(fft-direct)/np.linalg.norm(direct))
    assert fft_error < .01, fft_error
    max_step = (p.wavelength_nm/1000)/(2*metadata['plan']['width_um'])
    assert axis[1]-axis[0] < max_step
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(ASSETS/'case01_dense_far.npz', axis=axis, dsc_total=dsc, mask=mask)
    report = dict(date='2026-10-19', source=str(SOURCE.relative_to(ROOT)),
        source_metadata_sha256=source_fingerprint(), coarse_grid_shape=[19,72], coarse_step_deg=5,
        new_grid_shape=list(dsc.shape), evaluated_directions=len(directions), cartesian_step=float(axis[1]-axis[0]),
        finite_window_sampling_reference=float(max_step), coarse_direction_complex_relative_l2=relative,
        coarse_direction_checks=len(selection), surface_fft_vs_direct_relative_l2=fft_error,
        surface_direct_direction_checks=len(chosen), elapsed_seconds=time.perf_counter()-start,
        field_combination='A_whole_SPM + A_local_FEM_with - A_local_SPM',
        scope='Same saved physical fields. Local volume and triangle radiation evaluated at every plotted direction. '
              'Whole-surface Fourier transform uses the existing padded FFT with selected direct-sum checks. '
              'No intensity-image smoothing, new FEM solve, forced symmetry or ensemble averaging. '
              'This is an angular-visualization refinement, not a FEM/rough-surface convergence certificate.')
    (ASSETS/'case01_dense_far_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker()
    else:
        from fullwave_prototype.guarded_run import run_guarded
        OUTPUT.mkdir(parents=True, exist_ok=True)
        report = run_guarded(['D:/anaconda/envs/NIST2/python.exe', '-m',
                           'tools.refine_showcase_far_field', '--worker'], OUTPUT/'resources.json',
                           maximum=6., reserve=3., timeout=900.)
        if report['status'] != 'completed' or report['returncode'] != 0:
            raise SystemExit(report.get('returncode') or 1)


if __name__ == '__main__':
    main()
