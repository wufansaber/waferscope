"""Diagnose channel symmetry and single-realization haze, without model edits."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from optical_propagation import pupil_geometry
from fullwave_prototype.hybrid_run import load_sources, optical_settings
from fullwave_prototype.hybrid_radiation import HybridRadiation
from fullwave_prototype.imaging import image_from_amplitude
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness
from fullwave_prototype.surface_beam import surface_fields_order
from fullwave_prototype.surface_scattering import SurfaceRadiation
from fullwave_prototype.whole_surface import random_surface


SOURCE = Path('build/model_campaign_20261010/large_hybrid')
OUTPUT = Path('build/particle100_validation_20261016/symmetry_haze')


def image(amplitude, geometry, p, settings):
    local, directions, rotation, mask, du = geometry
    field = np.zeros(mask.shape+(2,), complex)
    field[mask] = amplitude
    result = image_from_amplitude(field, local, directions, rotation, du,
                                  p.wavelength_nm/1000, settings)
    return result['flux']*settings.irradiance_w_um2, result['axis_um']


def power(amplitude, geometry):
    local, _, _, mask, du = geometry
    return float(np.sum(np.sum(abs(amplitude)**2, axis=-1)/local[mask, 2])*du**2)


def flat_reference(p, settings):
    from pySCATMECH.local import Local_BRDF_Model
    model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {
        'lambda': p.wavelength_nm/1000, 'radius': p.diameter_nm/2000,
        'sphere': complex(p.particle_n, p.particle_k),
        'substrate': complex(p.substrate_n, p.substrate_k), 'delta': p.gap_nm/1000,
        'type': 0, 'lmax': 14, 'order': -1, 'Norm_Inc_Approx': 0})
    result = []
    for tilt in (45, -45):
        geometry = pupil_geometry(.5, tilt, 64)
        q = geometry[1][geometry[3]]
        amplitude = np.array([model.JonesDSC(np.deg2rad(p.incidence_deg),
            np.arccos(k[2]), np.arctan2(k[1], k[0])-np.pi/2, coords='psps')[:, 0 if p.polarization == 's' else 1]
            for k in q])
        result.append(power(amplitude, geometry))
    return dict(model='Flat wafer, plane wave Bobbert-Vlieger; not the Gaussian rough FEM case',
                dsc_um2=result, relative_difference=abs(result[0]-result[1])/np.mean(result))


def mirror_sources(data):
    result = dict(data)
    reflection = np.array([-1, 1, 1])
    for name in ('with', 'without'):
        for key in ('source_points_', 'source_polarization_weights_'):
            result[key+name] = data[key+name]*reflection
    result['local_spm_positions_um'] = data['local_spm_positions_um']*[-1, 1]
    result['local_spm_weights'] = data['local_spm_weights']*reflection
    # Origin-centered even grid: -x uses roll(flip(x), 1), not just flip(x).
    result['whole_weighted_fields'] = np.roll(data['whole_weighted_fields'][:, ::-1], 1, axis=1)*reflection
    return result


def components_check(data, p):
    original, mirrored = HybridRadiation(data, p), HybridRadiation(mirror_sources(data), p)
    powers = []
    first = None
    for tilt in (45, -45, 0):
        geometry = pupil_geometry(.5, tilt, 64)
        q = geometry[1][geometry[3]]
        terms = original.components(q)
        terms['particle_change'] = terms['local_fem_with']-terms['local_fem_without']
        terms['haze_local_correction'] = terms['local_fem_without']-terms['local_spm']
        values = {k: power(v, geometry) for k, v in terms.items()}
        values['integrated_interference'] = values['with_particle']-values['without']-values['particle_change']
        powers.append(values)
        if first is None:
            first = (q, terms)
    q, terms = first
    reflected = mirrored.components(q*[-1, 1, 1])
    errors = {}
    for name in ('with_particle', 'without', 'whole_spm', 'local_fem_with', 'local_fem_without'):
        a = np.sum(abs(terms[name])**2, axis=-1)
        b = np.sum(abs(reflected[name])**2, axis=-1)
        errors[name] = float(np.linalg.norm(a-b)/np.linalg.norm(a))
    return dict(component_dsc_um2=powers, mirror_radiation_intensity_relative_l2=errors,
        mirror_scope='Reflects archived quadrature sources, not an independent remesh or FEM solve')


def ensemble(data, metadata, p, settings, count=64):
    axis = data['whole_axis_um']
    background = surface_fields_order(axis, axis, p, p.beam_order)[0]
    beam_error = float(np.linalg.norm(background[:, 1:]-background[:, :0:-1]*[-1, 1, 1])/
                       np.linalg.norm(background[:, 1:]))
    taper = data['whole_height_window']
    roughness = Roughness(**metadata['stationary_roughness'])
    plan = dict(grid_points=len(axis), spacing_um=float(axis[1]-axis[0]))
    geometries = [pupil_geometry(.5, tilt, 64) for tilt in (45, -45, 0)]
    first, sums, square_sums, powers = [], [], [], []
    for n in range(count):
        height = data['whole_height_um'] if n == 0 else random_surface(roughness, plan, seed=10000+n)
        radiation = SurfaceRadiation(axis, axis, height[..., None]*taper[..., None]*background, p)
        row = []
        for i, geometry in enumerate(geometries):
            a = radiation.amplitude(geometry[1][geometry[3]])
            flux, image_axis = image(a, geometry, p, settings)
            row.append(power(a, geometry))
            if n == 0:
                first.append(flux.copy())
                sums.append(flux.copy())
                square_sums.append(flux**2)
            else:
                sums[i] += flux
                square_sums[i] += flux**2
        powers.append(row)
        if (n+1) % 16 == 0:
            print('Haze-only SPM ensemble: {}/{}'.format(n+1, count), flush=True)
    powers = np.asarray(powers)
    mean = np.array(sums)/count
    variance = (np.array(square_sums)-count*mean**2)/(count-1)
    selected = abs(image_axis) <= 10
    cv = []
    for i in range(3):
        one = first[i][np.ix_(selected, selected)]
        average = mean[i][np.ix_(selected, selected)]
        cv.append(dict(single_spatial_cv=float(one.std()/one.mean()),
                       mean_spatial_cv=float(average.std()/average.mean())))
    difference = powers[:, 0]-powers[:, 1]
    stats = dict(surfaces=count, model='Whole-surface first-order SPM haze ONLY; no local FEM correction',
        illumination_x_reflection_relative_l2=beam_error,
        mean_dsc_um2=powers.mean(axis=0).tolist(), std_dsc_um2=powers.std(axis=0, ddof=1).tolist(),
        paired_difference_mean=float(difference.mean()),
        paired_difference_95ci_halfwidth=float(1.998*difference.std(ddof=1)/np.sqrt(count)),
        central_20um_spatial_cv=cv,
        scope='Independent synthetic surface intensities averaged, NOT fields. Not an actual exposure '
              'unless the measurement averages corresponding decorrelated surface/illumination states. '
              'Spatial CV includes the illumination envelope, not a detrended pure speckle contrast.')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTPUT/'haze_ensemble.npz', image_axis_um=image_axis, single=np.array(first),
                        mean=mean, variance=variance, powers_dsc_um2=powers)
    plot_haze(image_axis, first, mean, count)
    return stats


def plot_haze(axis, single, mean, count):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected = abs(axis) <= 30
    maximum = max(float(np.max(a)) for a in single)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    for row, images in enumerate((single, mean)):
        for i, values in enumerate(images):
            im = axes[row, i].pcolormesh(axis[selected], axis[selected],
                values[np.ix_(selected, selected)], vmin=0, vmax=maximum, shading='auto', cmap='inferno')
            axes[row, i].set(title='Ch{} | {}'.format(i+1, 'One surface' if row == 0 else '{}-surface mean'.format(count)),
                            xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
    fig.colorbar(im, ax=list(axes.ravel()), label='Haze W/um^2 | same absolute scale for all panels')
    fig.suptitle('266 nm | 1 W | unchanged coherent beam | whole SPM only (not averaged hybrid FEM)')
    fig.savefig(OUTPUT/'haze_single_vs_ensemble.png', dpi=170)
    plt.close(fig)


def main():
    data, metadata = load_sources(SOURCE)
    p = Parameters(**metadata['parameters'])
    values = json.loads((SOURCE/'hybrid_input.json').read_text(encoding='utf-8'))['imaging_settings']
    settings = replace(optical_settings(values, p), pupil_samples=64, padding=4)
    result = dict(parameters=metadata['parameters'], roughness=metadata['stationary_roughness'],
        flat_reference=flat_reference(p, settings), components=components_check(data, p),
        ensemble=ensemble(data, metadata, p, settings))
    (OUTPUT/'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result['components'], indent=2), flush=True)


if __name__ == '__main__':
    main()
