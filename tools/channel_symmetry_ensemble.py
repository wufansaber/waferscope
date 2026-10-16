"""Fixed-sample independent-surface haze symmetry, with paired confidence intervals."""
import gc
import numpy as np

from fullwave_prototype.rough_surface import Roughness
from fullwave_prototype.surface_scattering import SurfaceRadiation
from fullwave_prototype.whole_surface import random_surface
from optical_propagation import pupil_geometry
from .channel_symmetry import channel_image, paired_summary


def run_ensemble(data, metadata, background, p, settings, camera, count, seed):
    if count < 16:
        raise ValueError('Use at least 16 independent surfaces for this audit')
    rough = Roughness(**metadata['stationary_roughness'])
    axis = data['whole_axis_um']
    plan = dict(grid_points=len(axis), spacing_um=float(axis[1]-axis[0]))
    geometries = [pupil_geometry(.5, tilt, settings.pupil_samples) for tilt in (45, -45)]
    powers = np.empty((count, 2, 2))
    mean_images = None
    first_images = None
    # Fix sample count and seeds before seeing any result. No antithetic/mirror pairing.
    for n in range(count):
        height = random_surface(rough, plan, seed=seed+n)
        fields = height[..., None]*data['whole_height_window'][..., None]*background
        radiation = SurfaceRadiation(axis, axis, fields, p)
        frames = []
        for i, geometry in enumerate(geometries):
            amplitude = radiation.amplitude(geometry[1][geometry[3]])
            metrics, image = channel_image(amplitude, geometry, p, settings, camera)
            powers[n, i] = metrics['aperture_power_w'], metrics['roi_power_w']
            frames.append(image['flux'])
        if n == 0:
            first_images = np.array(frames)
            mean_images = np.zeros_like(first_images)
            image_axis = image['axis_um']
        mean_images += frames
        del radiation, fields, height
        gc.collect()
        if (n+1) % 8 == 0:
            print('Independent haze surfaces: {}/{}'.format(n+1, count), flush=True)
    summary = {}
    for i, name in enumerate(('aperture', 'camera_roi')):
        summary[name] = paired_summary(powers[:, :, i])
        # Two predeclared power comparisons: Bonferroni simultaneous 95% coverage.
        summary[name]['simultaneous_95pct'] = paired_summary(powers[:, :, i], confidence=.975)
    summary.update(surfaces=count, seed_start=seed, seed_stop_inclusive=seed+count-1,
        model='Whole-surface first-order SPM, haze only',
        scope='Intensities, not fields, averaged over independent PSD surfaces at identical '
              'illumination. Does not validate a full independent-FEM particle/roughness ensemble. '
              'No detector noise, no PRNU, no forced mirroring. Repeated frames of a static '
              'surface do not represent this ensemble.')
    return summary, dict(powers_w=powers, image_axis_um=image_axis,
                         first_images=first_images, mean_images=mean_images/count)
