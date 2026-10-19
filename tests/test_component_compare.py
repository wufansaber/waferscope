"""Component-comparison settings and plotting without a heavy FEM solve."""
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from fullwave_prototype.component_compare import _clean_source, _settings
from fullwave_prototype.component_compare_plotting import _component_norm, plot_comparison
from fullwave_prototype.plot_directions import mark_illumination_source
from fullwave_prototype.imaging import ImagingSettings
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness


class ComponentCompareTests(unittest.TestCase):
    def test_weak_component_uses_its_own_numeric_range(self):
        weak = _component_norm(np.array([1e-9]), logarithmic=False)
        strong = _component_norm(np.array([1e-6]), logarithmic=False)
        self.assertEqual(weak.vmax, 1e-9)
        self.assertEqual(strong.vmax, 1e-6)
        self.assertAlmostEqual(weak(1e-9), 1)
        self.assertLess(strong(1e-9), .04)
        self.assertAlmostEqual(_component_norm(np.array([1e-9]), True)(1e-9), 1)

    def test_source_direction_is_minus_ky(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        mark_illumination_source(ax, 75)
        azimuth, radius = ax.collections[-1].get_offsets()[0]
        self.assertAlmostEqual(azimuth, -np.pi/2)
        self.assertAlmostEqual(radius, np.sin(np.deg2rad(75)))
        self.assertIn('-KY', ax.texts[-1].get_text())
        plt.close(fig)

    def test_saved_optics_reconstruct_without_sampling_camera_noise(self):
        saved = asdict(ImagingSettings(numerical_apertures=(0, .4125, .96)))
        saved['detector']['enabled'] = True
        settings = _settings(saved)
        self.assertEqual(settings.numerical_apertures, (0, .4125, .96))
        self.assertFalse(settings.detector.enabled)
        self.assertEqual(settings.detector.exposure_ms, saved['detector']['exposure_ms'])

    def test_valid_smooth_source_is_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'scattering'
            source.mkdir()
            clean = root/'particle_only_scattering'
            clean.mkdir()
            p = Parameters()
            r = Roughness()
            report = dict(source_format_version=1, parameters=asdict(p),
                          roughness=asdict(Roughness(rms_nm=0)), angular_step_deg=5,
                          metrics={'with': {'source_quadrature_order': 7}})
            (clean/'report.json').write_text(json.dumps(report), encoding='utf-8')
            (clean/'paired_fields.npz').touch()
            with patch('fullwave_prototype.component_compare.run_scattering') as run:
                selected, saved = _clean_source(root, source, p, r, 5, 7)
            self.assertEqual(selected, clean)
            self.assertEqual(saved, report)
            run.assert_not_called()

    def test_component_scales_and_export_totals(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            theta = np.linspace(0, np.pi/2, 19)
            phi = np.linspace(0, 2*np.pi, 72, endpoint=False)
            far = dict(theta=theta, phi=phi,
                       particle_dsc=np.ones((19, 72))*1e-8,
                       haze_dsc=np.ones((19, 72))*1e-6)
            axis = np.linspace(-20, 20, 51)
            fields = {}
            channels = []
            for index in (1, 2, 3):
                fields[f'ch{index}_axis_um'] = axis
                fields[f'ch{index}_particle_w_um2'] = np.ones((51, 51))*index*1e-18
                fields[f'ch{index}_haze_w_um2'] = np.ones((51, 51))*index*1e-16
                channels.append(dict(channel=index, na=.5, tilt_deg=(45, -45, 0)[index-1],
                                     particle_power_w=index*1e-15, haze_power_w=index*1e-13,
                                     particle_energy_j=index*1e-18, haze_energy_j=index*1e-16))
            report = dict(parameters=dict(incidence_deg=75),
                          imaging_settings=dict(magnification=20,
                              numerical_apertures=[.5]*3, tilt_degrees=[45, -45, 0]),
                          far_hemisphere=dict(particle_dsc_um2=1e-8, haze_dsc_um2=1e-6,
                              particle_equivalent_power_w=1e-14,
                              haze_equivalent_power_w=1e-12),
                          exposure_ms=1., channels=channels)
            plot_comparison(output, far, fields, report)
            for name in ('far_components.png', 'far_components_log.png',
                         'image_components.png', 'image_components_log.png'):
                self.assertGreater((output/name).stat().st_size, 1000)
            totals = (output/'totals.txt').read_text(encoding='utf-8')
            self.assertIn('Ch3', totals)
            self.assertIn('Haze (fJ)', totals)
            self.assertIn('separate numerical color scales', totals)


if __name__ == '__main__':
    unittest.main()
