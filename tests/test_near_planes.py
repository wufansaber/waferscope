"""Orthogonal sampling, archive compatibility and shared-scale near plots."""
from dataclasses import asdict, replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.near_planes import (make_planes, sample_planes, finish_planes,
                                          saved_planes, electric_field, PLANE_AXES)
from fullwave_prototype.near_plotting import (surface_segments, shared_norm,
                                            draw_plane, plot_near_planes)
from fullwave_prototype.rough_surface import Roughness, height_mesh


class NearPlaneTests(unittest.TestCase):
    def setUp(self):
        self.p = Parameters(gap_nm=10)

    def data(self):
        planes = make_planes(self.p, samples=11)
        data = {}

        def evaluate(mesh, field, points, p):
            e = points.astype(complex)+field*(1+2j)
            return e, 2j*e, np.zeros_like(e)

        with patch('fullwave_prototype.solver.evaluate_fields', side_effect=evaluate) as evaluator:
            sample_planes(None, 1, self.p, planes, 'with', data)
            sample_planes(None, 0, self.p, planes, 'without', data)
        self.assertEqual(evaluator.call_count, 6)
        finish_planes(data, planes)
        data['surface_vertices_um'], data['surface_triangles_um'] = height_mesh(Roughness())
        return planes, data

    def test_coordinates_order_and_physical_domain(self):
        for p in (self.p, replace(self.p, diameter_nm=10, gap_nm=205)):
            planes = make_planes(p)
            self.assertEqual([s.name for s in planes], ['xz', 'xy', 'yz'])
            for plane in planes:
                points = plane.points.reshape(101, 101, 3)
                horizontal, vertical, fixed = PLANE_AXES[plane.name]
                np.testing.assert_array_equal(points[0, :, horizontal], plane.u)
                np.testing.assert_array_equal(points[:, 0, vertical], plane.v)
                np.testing.assert_array_equal(points[..., fixed], plane.offset)
                self.assertLess(np.linalg.norm(points, axis=-1).max(), p.domain_radius_um)
            self.assertAlmostEqual(planes[1].offset, (p.diameter_nm/2+p.gap_nm)/1000)

    def test_native_samples_complex_fields_and_legacy_aliases(self):
        planes, data = self.data()
        for plane in planes:
            prefix = 'near_'+plane.name+'_'
            expected = plane.points.reshape(11, 11, 3)+(1+2j)
            np.testing.assert_allclose(data[prefix+'e_with'], expected)
            np.testing.assert_allclose(data[prefix+'e_change'], 1+2j)
            np.testing.assert_allclose(data[prefix+'hz0_change'], 2j*(1+2j))
            self.assertEqual(plane.metadata()['shape'], [11, 11])
        np.testing.assert_array_equal(data['near_e_with'], data['near_yz_e_with'])
        np.testing.assert_array_equal(data['near_e_change'], data['near_yz_e_change'])
        self.assertEqual([p.name for p in saved_planes(data)], ['xz', 'xy', 'yz'])
        self.assertFalse(np.array_equal(data['near_xz_e_with'], data['near_yz_e_with']))

    def test_legacy_archive_does_not_invent_missing_planes(self):
        _, data = self.data()
        legacy = {k: v for k, v in data.items()
                  if not any(k.startswith('near_'+p+'_') for p in PLANE_AXES)}
        self.assertEqual([p.name for p in saved_planes(legacy)], ['yz'])
        with tempfile.TemporaryDirectory() as directory:
            plot_near_planes(directory, legacy, {'parameters': asdict(self.p)})
            self.assertTrue((Path(directory)/'near_field_yz.png').is_file())
            self.assertFalse((Path(directory)/'near_field_xy.png').exists())

    def test_triangle_intersection_includes_diagonal_kink(self):
        a, b, c, d = np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 1], [-1, 1, 0.]])
        segments = surface_segments(np.array([[a, b, c], [a, c, d]]), 0, 0)
        self.assertEqual(len(segments), 2)
        hits = np.concatenate(segments)
        self.assertTrue(any(np.allclose(point, [0, 0, .5]) for point in hits))
        np.testing.assert_array_equal(hits[:, 0], 0)

    def test_shared_scale_correct_particle_centers_and_export(self):
        from fullwave_prototype.imaging_plotting import _pyplot
        plt = _pyplot()
        planes, data = self.data()
        norm = shared_norm(data, planes)
        maximum = max(np.sum(abs(electric_field(data, p))**2, axis=-1).max() for p in planes)
        self.assertEqual(norm.vmax, maximum)
        report = {'parameters': asdict(self.p)}
        for plane in planes:
            fig, ax = plt.subplots()
            im = draw_plane(ax, data, report, plane, norm)
            self.assertIs(im.norm, norm)
            np.testing.assert_allclose(ax.patches[0].center, (0, 0 if plane.name == 'xy' else 60))
            self.assertEqual(ax.patches[0].radius, 50)
            plt.close(fig)
        with tempfile.TemporaryDirectory() as directory:
            plot_near_planes(directory, data, report)
            for name in ('near_field', 'near_field_xz', 'near_field_xy', 'near_field_yz'):
                self.assertGreater((Path(directory)/(name+'.png')).stat().st_size, 1000)


if __name__ == '__main__':
    unittest.main()
