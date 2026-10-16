"""Smooth-wafer and mirrored hybrid-source checks, without changing the model."""
import gc
import numpy as np

from fullwave_prototype.far_field import radiation_amplitude
from fullwave_prototype.hybrid_radiation import HybridRadiation
from optical_propagation import pupil_geometry
from .channel_symmetry import channel_image, reflect_sources, relative_error, reflect_centered_grid


def smooth_check(native, p, settings, camera):
    from pySCATMECH.local import Local_BRDF_Model
    model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {
        'lambda': p.wavelength_nm/1000, 'radius': p.diameter_nm/2000,
        'sphere': complex(p.particle_n, p.particle_k),
        'substrate': complex(p.substrate_n, p.substrate_k), 'delta': p.gap_nm/1000,
        'type': 0, 'lmax': 14, 'order': -1, 'Norm_Inc_Approx': 0})
    values, images = {}, {}
    for name in ('plane_Bobbert_Vlieger', 'gaussian_FEM'):
        pairs, pictures = [], []
        for tilt in (45, -45):
            geometry = pupil_geometry(.5, tilt, settings.pupil_samples)
            directions = geometry[1][geometry[3]]
            if name == 'gaussian_FEM':
                amplitude = radiation_amplitude(native['source_points'],
                    native['source_polarization_weights'], directions, p)
            else:
                amplitude = np.array([model.JonesDSC(np.deg2rad(p.incidence_deg),
                    np.arccos(q[2]), np.arctan2(q[1], q[0])-np.pi/2,
                    coords='psps')[:, 0 if p.polarization == 's' else 1] for q in directions])
            metrics, image = channel_image(amplitude, geometry, p, settings, camera)
            pairs.append(metrics)
            pictures.append(image)
        errors = {key: abs(pairs[0][key]-pairs[1][key])/np.mean([v[key] for v in pairs])
                  for key in pairs[0]}
        errors['mirrored_image_relative_l2'] = relative_error(
            pictures[0]['flux'], reflect_centered_grid(pictures[1]['flux']))
        values[name] = dict(channels=pairs, errors=errors,
            power_symmetry_within_1pct=bool(max(errors[k] for k in pairs[0]) < .01))
        images[name] = pictures
    values['scope'] = ('Plane-wave BV and current finite-Gaussian FEM are separate symmetry tests, '
        'not an absolute-power cross-validation between different illuminations. '
        'A smooth wafer without a particle has zero diffuse source; specular reflection is excluded.')
    return values, images


def hybrid_measurements(data, p, settings, camera):
    radiation = HybridRadiation(data, p)
    measurements, pictures = [], []
    for tilt in (45, -45):
        geometry = pupil_geometry(.5, tilt, settings.pupil_samples)
        terms = radiation.components(geometry[1][geometry[3]])
        terms['particle_change'] = terms['with_particle']-terms['without']
        row, image_row = {}, {}
        for state in ('with_particle', 'without', 'particle_change'):
            row[state], image_row[state] = channel_image(terms[state], geometry, p, settings, camera)
        row['interference'] = {key: row['with_particle'][key]-row['without'][key]-row['particle_change'][key]
                               for key in row['with_particle']}
        measurements.append(row)
        pictures.append(image_row)
    del radiation
    gc.collect()
    return measurements, pictures


def compare_mirror(original, mirrored, original_images, mirrored_images):
    result = []
    for index in (0, 1):
        errors = {}
        for state in ('with_particle', 'without', 'particle_change'):
            a, b = original[index][state], mirrored[1-index][state]
            errors[state] = {key: abs(a[key]-b[key])/max(abs(a[key]), 1e-300) for key in a}
            left, right = original_images[index][state], mirrored_images[1-index][state]
            errors[state].update(
                pupil_intensity_relative_l2=relative_error(left['pupil_intensity'], right['pupil_intensity'][:, ::-1]),
                image_relative_l2=relative_error(left['flux'], reflect_centered_grid(right['flux'])),
                pixel_power_relative_l2=relative_error(left['pixels'], right['pixels'][:, ::-1]))
        result.append(errors)
    return dict(original_to_opposite_channel=result,
        power_swap_within_1pct=bool(max(row[state][key] for row in result
            for state in row for key in ('aperture_power_w', 'roi_power_w')) < .01))


def mirror_check(data, fresh, background, p, settings, camera):
    print('Integrating original hybrid sources...', flush=True)
    original, images = hybrid_measurements(data, p, settings, camera)
    print('Integrating reflected archived sources (algebraic check)...', flush=True)
    reflected = reflect_sources(data)
    algebraic, algebraic_images = hybrid_measurements(reflected, p, settings, camera)
    exact = compare_mirror(original, algebraic, images, algebraic_images)
    del algebraic_images
    # Re-illuminate the reflected height with the unchanged beam, rather than just mirroring its field.
    physical_fields = (reflect_centered_grid(data['whole_height_um']*data['whole_height_window'])
                       [..., None]*background)
    field_error = relative_error(physical_fields, reflected['whole_weighted_fields'])
    reflected['whole_weighted_fields'] = physical_fields
    for state in ('with', 'without'):
        for prefix in ('source_points_', 'source_polarization_weights_'):
            reflected[prefix+state] = fresh[prefix+state]
    print('Integrating independently solved mirrored hybrid sources...', flush=True)
    independent, independent_images = hybrid_measurements(reflected, p, settings, camera)
    remesh = compare_mirror(original, independent, images, independent_images)
    return dict(original_channels=original, reflected_archived_channels=algebraic,
        independently_solved_mirror_channels=independent, algebraic_mirror=exact,
        independent_remesh_mirror=remesh, unchanged_beam_reflection_relative_l2=field_error,
        scope='Algebraic check reflects saved volume sources. Independent check rebuilds exact '
              'mirrored triangular CAD and re-solves both local Maxwell states; the whole surface '
              'is still first-order SPM, not full-wafer FEM. No camera noise or PRNU is used.'), images, independent_images
