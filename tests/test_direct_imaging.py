"""Independent diffraction: rotations, Hankel theory and production FFT comparisons."""
import unittest
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import numpy as np

from fullwave_prototype.direct_imaging import direct_image, polar_cap, ray_frame, rotate_between


class DirectImagingTests(unittest.TestCase):
    def test_saved_audit_detects_corrupted_complex_field(self):
        from fullwave_prototype.imaging import ImagingSettings, calculate_images
        from fullwave_prototype.parameters import Parameters
        from fullwave_prototype.imaging_audit import check_saved
        p = Parameters(substrate_n=1, substrate_k=0)
        settings = ImagingSettings(pupil_samples=32, padding=4, aberrations_nm=(15, 25, 10))
        data = dict(source_points_with=np.zeros((1, 3)), source_points_without=np.zeros((1, 3)),
                    source_polarization_weights_with=np.array([[.3, .2j, 1]]),
                    source_polarization_weights_without=np.array([[.2j, .3, .5]]))
        arrays, report = calculate_images(data, p, settings)
        report['scattering_parameters'] = asdict(p)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            np.savez(path/'image_fields.npz', **arrays)
            actual, _ = check_saved(path)
            self.assertTrue(all(c['passed'] for c in actual['channels']))
            compact = {k: v for k, v in arrays.items() if '_image_e_' not in k and '_image_hz0_' not in k}
            np.savez(path/'image_fields.npz', **compact)
            compact_report, _ = check_saved(path)
            self.assertTrue(all(c['passed'] and not c['complex_fields_checked'] for c in compact_report['channels']))
            arrays['ch3_image_e_with'] *= 1.01
            np.savez(path/'image_fields.npz', **arrays)
            actual, _ = check_saved(path)
            self.assertFalse(actual['channels'][2]['passed'])
            arrays['ch3_image_hz0_with'][...] = np.nan
            np.savez(path/'image_fields.npz', **arrays)
            with self.assertRaisesRegex(ValueError, 'Non-finite'):
                check_saved(path)

    def test_rotations_and_disk_measure(self):
        uv, weight, _ = polar_cap(.8, 0, 24, 96)
        local, _, _, _, _ = ray_frame(uv, 0)
        z = np.broadcast_to((0., 0., 1.), local.shape)
        np.testing.assert_allclose(rotate_between(local, local, z), z, atol=1e-15)
        self.assertAlmostEqual(weight.sum(), np.pi*.8**2, places=13)
        uv, weight, _ = polar_cap(0, 0)
        e, h, flux = direct_image(np.empty((0, 2)), uv, weight, [[0, 0]], .355)
        np.testing.assert_array_equal(flux, 0)

    def test_direct_against_bessel_hankel_theory(self):
        from fullwave_prototype.theory_optics import hankel_slice
        uv, weight, _ = polar_cap(.8, 0, 64, 192)
        local, _, s, p, _ = ray_frame(uv, 0)
        z = np.broadcast_to((0., 0., 1.), local.shape)
        # Invert the ideal objective: a uniformly x-polarized collimated disk.
        far = rotate_between(np.broadcast_to((1., 0., 0.), local.shape), z, local)
        far *= np.sqrt(local[:, 2])[:, None]
        amplitude = np.column_stack((np.sum(far*s, axis=1), np.sum(far*p, axis=1)))
        for m, dz in ((1, 0), (1, .4), (20, 0)):
            x = np.linspace(-1.5, 1.5, 31)*m
            e, h, _ = direct_image(amplitude, uv, weight, np.column_stack((x, np.zeros_like(x))),
                                   .355, magnification=m, defocus_um=dz)
            expected_e, expected_h = hankel_slice(x, .8, .355, m, dz)
            np.testing.assert_allclose(e, expected_e, rtol=1e-10, atol=1e-13)
            np.testing.assert_allclose(h, expected_h, rtol=1e-10, atol=1e-13)

    def test_fft_absolute_complex_fields(self):
        from optical_propagation import pupil_geometry
        from fullwave_prototype.imaging import ImagingSettings, image_from_amplitude
        rng = np.random.default_rng(123)
        for na, tilt, m, analyzer in ((.5, 0, 20, 'none'), (.8, 45, 1, 'x'),
                                      (.5, -45, 20, '45'), (.5, 20, 2, 'y')):
            local, k, rotation, mask, du = pupil_geometry(na, tilt, 32)
            amp = rng.normal(size=mask.shape+(2,))+1j*rng.normal(size=mask.shape+(2,))
            amp[~mask] = 0
            settings = ImagingSettings(pupil_samples=32, padding=4, magnification=m,
                                       analyzer=analyzer, image_defocus_um=.6, transmission=.7)
            actual = image_from_amplitude(amp, local, k, rotation, du, .355, settings)
            indices = np.array((37, 59, 64, 69, 93))
            yy, xx = np.meshgrid(indices, indices, indexing='ij')
            xy = np.column_stack((actual['axis_um'][xx.ravel()], actual['axis_um'][yy.ravel()]))
            e, h, flux = direct_image(amp[mask], local[mask, :2], np.full(mask.sum(), du**2),
                xy, .355, tilt, m, .6, .7, analyzer)
            np.testing.assert_allclose(e, actual['e'][yy, xx].reshape(-1, 3), rtol=1e-11, atol=1e-14)
            np.testing.assert_allclose(h, actual['hz0'][yy, xx].reshape(-1, 3), rtol=1e-11, atol=1e-14)
            np.testing.assert_allclose(flux, actual['flux'][yy, xx].ravel(), rtol=1e-11, atol=1e-15)

    def test_invalid_inputs(self):
        for uv in ([[1, 0]], [[float('nan'), 0]]):
            with self.assertRaises(ValueError):
                ray_frame(uv, 0)
        with self.assertRaises(ValueError):
            direct_image([[1, 0]], [[.1, 0]], [-1], [[0, 0]], .355)


if __name__ == '__main__':
    unittest.main()
