"""Independent Parseval, covariance and experimental preprocessing checks."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from fullwave_prototype.experimental_input import prepare
from fullwave_prototype.experimental_surface import (height_periodogram, load_height,
    load_psd, measurement_report, measured_nodes, psd_covariance_kernel, sample_psd)
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, height_mesh
from fullwave_prototype.rough_psd import expected_finite_psd, height_basis_transform


class ExperimentalSurfaceTests(unittest.TestCase):
    def test_plane_removal_units_and_parseval(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'height.npz'
            y, x = np.meshgrid(np.arange(9)*.025, np.arange(13)*.025, indexing='ij')
            height = 5+3*x-2*y+np.cos(2*np.pi*x/.1)
            np.savez(source, height_nm=height, dx_um=.025, dy_um=.025)
            processed, report = load_height(source)
            self.assertAlmostEqual(processed.mean(), 0, places=13)
            np.testing.assert_allclose(report['footprint_xy_um'], [.3, .2])
            fx, fy, psd = height_periodogram(processed, .025, .025)
            variance = psd.sum()*(fx[1]-fx[0])*(fy[1]-fy[0])
            self.assertAlmostEqual(variance, np.mean(processed**2)*1e-6, places=18)
            output = Path(directory)/'audit'
            audited = prepare(source, output)
            self.assertLess(audited['parseval_relative_error'], 1e-12)
            load_psd(output/'measured_psd.npz')

    def test_missing_values_need_explicit_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'height.npz'
            height = np.ones((9, 9))
            height[4, 4] = np.nan
            np.savez(source, height_nm=height, dx_um=.025, dy_um=.025)
            with self.assertRaisesRegex(ValueError, 'missing'):
                load_height(source)
            filled, report = load_height(source, missing='nearest')
            self.assertTrue(np.isfinite(filled).all())
            self.assertEqual(report['missing_points'], 1)

    def test_exported_afm_sized_periodogram_can_be_reimported(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'psd.npz'
            height = np.random.default_rng(8).normal(size=(512, 512))
            fx, fy, density = height_periodogram(height, .001, .001)
            np.savez(source, fx_per_um=fx, fy_per_um=fy, psd_um4=density)
            _, _, _, report = load_psd(source)
            self.assertAlmostEqual(report['rms_nm'], np.sqrt(np.mean(height**2)), places=12)

    def test_measured_height_is_absolute_and_seed_independent(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'height.npz'
            height = np.arange(81).reshape(9, 9)**2*.001
            np.savez(source, height_nm=height, dx_um=.025, dy_um=.025)
            r = Roughness(model='measured_height', measurement_path=str(source), detrend='mean',
                          grid_points=9, taper_fraction=.125)
            r = replace(r, rms_nm=measurement_report(r)['rms_nm'])
            r.validate(Parameters(gap_nm=10), particle_enabled=False)
            np.testing.assert_allclose(measured_nodes(r).T,
                                       (height-height.mean())/1000, atol=1e-15)
            np.testing.assert_array_equal(height_mesh(r)[0], height_mesh(replace(r, seed=99))[0])
            with self.assertRaisesRegex(ValueError, 'without amplitude'):
                replace(r, rms_nm=1).validate(Parameters(), particle_enabled=False)
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                height_mesh(replace(r, width_um=.3))

    def test_psd_covariance_and_finite_reference_match_direct_quadratic_form(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'psd.npz'
            f = np.arange(-4, 5)*2.
            fx, fy = np.meshgrid(f, f)
            psd = 2e-8*np.exp(-(.2*fx)**2-(.4*fy)**2)
            np.savez(source, fx_per_um=f, fy_per_um=f, psd_um4=psd)
            r = Roughness(model='measured_psd', measurement_path=str(source), grid_points=5,
                          taper_fraction=.25)
            metadata = measurement_report(r)
            r = replace(r, rms_nm=metadata['rms_nm'])
            kernel = psd_covariance_kernel(r)
            self.assertAlmostEqual(kernel[4, 4], psd.sum()*4, places=18)
            self.assertNotAlmostEqual(kernel[5, 4], kernel[4, 5], places=10)
            theta, phi = np.deg2rad([0, 45, 90]), np.deg2rad([0, 90, 180, 270])
            polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
            frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
                np.sin(polar).ravel()*np.sin(azimuth).ravel()-np.sin(np.deg2rad(75))))/.266
            from fullwave_prototype.rough_surface import _taper_grid, taper_effective_area
            basis = height_basis_transform(r, frequency)*_taper_grid(r).ravel()
            x, y = np.meshgrid(np.arange(5), np.arange(5), indexing='ij')
            x, y = x.ravel(), y.ravel()
            covariance = kernel[x[:, None]-x[None, :]+4, y[:, None]-y[None, :]+4]
            direct = np.einsum('fi,ij,fj->f', basis.conj(), covariance, basis).real/taper_effective_area(r)
            np.testing.assert_allclose(expected_finite_psd(r, theta, phi, .266, 75).ravel(),
                                       direct, rtol=1e-10, atol=1e-20)
            self.assertEqual(float(sample_psd(source, np.array([[100, 0]]))[0]), 0.)
            psd[0, 0] *= 2
            np.savez(source, fx_per_um=f, fy_per_um=f, psd_um4=psd)
            with self.assertRaisesRegex(ValueError, 'real height'):
                load_psd(source)

    def test_radial_psd_variance_against_analytic_integral(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'psd.csv'
            f = np.linspace(0, 20, 201)
            a, b = 2e-8, .5
            np.savetxt(source, np.column_stack((f, a*np.exp(-(b*f)**2))), delimiter=',')
            _, _, _, report = load_psd(source)
            expected = np.pi*a/b**2
            self.assertAlmostEqual((report['rms_nm']/1000)**2/expected, 1, delta=.004)
            np.savetxt(source, np.column_stack((f[10:], a*np.exp(-(b*f[10:])**2))), delimiter=',')
            _, _, _, report = load_psd(source)
            self.assertEqual(report['measured_radial_band_per_um'][0], 1.)
            self.assertEqual(sample_psd(source, [[.1, 0]])[0], 0.)


if __name__ == '__main__':
    unittest.main()
