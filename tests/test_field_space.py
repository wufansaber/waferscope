import importlib.util
from dataclasses import replace
import unittest

from fullwave_prototype.parameters import Parameters


class FieldSpaceParameters(unittest.TestCase):
    def test_old_reports_and_order_validation(self):
        self.assertEqual(Parameters(**{'order': 3}).pml_order, 0)
        for order in (-1, 4, 1.5):
            with self.assertRaises(ValueError):
                Parameters(order=3, pml_order=order).validate()
        Parameters(order=3, pml_order=2).validate()
        for value in (-.1, 1, float('nan')):
            with self.assertRaises(ValueError):
                Parameters(pml_inner_fraction=value).validate()

    def test_controls_affect_local_mesh_only(self):
        from fullwave_prototype.imaging_run import build_parser, scattering_settings
        args = build_parser().parse_args(['--accuracy', 'balanced', '--pml-order', '2', '--local-mesh-scale', '.8'])
        p, r = scattering_settings(args)
        self.assertEqual(p.pml_order, 2)
        self.assertAlmostEqual(p.particle_mesh_um, .012)
        self.assertAlmostEqual(r.mesh_um, .020)
        self.assertEqual((p.air_mesh_um, p.substrate_mesh_um, p.pml_mesh_scale), (.08, .04, 1.25))


@unittest.skipIf(importlib.util.find_spec('ngsolve') is None, 'Native NGSolve required')
class NativeFieldSpaceTests(unittest.TestCase):
    def test_physical_cells_and_shared_entities_keep_full_order(self):
        import ngsolve as ng
        from fullwave_prototype.field_space import make_space
        from fullwave_prototype.solver import build_mesh
        p = Parameters(order=3, gap_nm=10, domain_radius_um=.14, pml_thickness_um=.10,
                       air_mesh_um=.1, substrate_mesh_um=.08, particle_mesh_um=.035)
        mesh = build_mesh(p)
        uniform = make_space(mesh, p)
        reduced = make_space(mesh, replace(p, pml_order=2))
        self.assertLess(reduced.ndof, uniform.ndof)
        pml_reduced = 0
        for cell in mesh.Elements(ng.VOL):
            if cell.mat.startswith('pml_'):
                pml_reduced += reduced.GetFE(cell).ndof < uniform.GetFE(cell).ndof
            else:
                self.assertEqual(reduced.GetFE(cell).ndof, uniform.GetFE(cell).ndof)
                for edge in cell.edges:
                    self.assertEqual(reduced.GetOrder(edge), 3)
                for face in cell.faces:
                    self.assertEqual(reduced.GetOrder(face), 3)
        self.assertGreater(pml_reduced, 0)
        import numpy as np
        protected = make_space(mesh, replace(p, pml_order=2, pml_inner_fraction=.5))
        self.assertGreater(protected.ndof, reduced.ndof)
        self.assertLess(protected.ndof, uniform.ndof)
        radius = p.domain_radius_um+.5*p.pml_thickness_um
        for cell in mesh.Elements(ng.VOL):
            center = np.mean([mesh[vertex].point for vertex in cell.vertices], axis=0)
            if not cell.mat.startswith('pml_') or np.linalg.norm(center) <= radius:
                self.assertEqual(protected.GetFE(cell).ndof, uniform.GetFE(cell).ndof)


if __name__ == '__main__':
    unittest.main()
