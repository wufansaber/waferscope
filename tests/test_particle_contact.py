"""The default flat-wafer geometry is contact, not a fixed-height sphere."""
from dataclasses import asdict, replace
import os
import unittest

from fullwave_prototype.defects import center_height_um, metadata
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.parameters import Parameters

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


class ParticleContactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from PyQt5.QtWidgets import QApplication
        cls.app = QApplication.instance() or QApplication([])

    def test_imaging_entry_defaults_to_flat_wafer_contact(self):
        p, roughness = scattering_settings(build_parser().parse_args([]))
        self.assertEqual((p.gap_nm, roughness.rms_nm), (0., 0.))

    def test_diameter_changes_center_not_bottom(self):
        for diameter in range(10, 101, 10):
            p = Parameters(diameter_nm=diameter)
            self.assertEqual(metadata(p)['bottom_z_nm'], 0.)
            self.assertAlmostEqual(center_height_um(p), diameter/2000.)

    def test_explicit_nonzero_gap_remains_reproducible(self):
        p, _ = scattering_settings(build_parser().parse_args(['--gap-nm', '10']))
        self.assertEqual(p.gap_nm, 10.)

    def test_rough_surface_overlap_is_not_silently_repaired(self):
        with self.assertRaisesRegex(ValueError, 'highest roughness peak'):
            scattering_settings(build_parser().parse_args(['--rms-nm', '2']))

    def test_gui_contact_presets_and_saved_geometry(self):
        from fullwave_ui.controls import SimulationControls
        controls = SimulationControls()
        try:
            for preset in (0, 1, 3, 4, 5):
                controls.preset.setCurrentIndex(preset)
                self.assertEqual(controls.fields['gap-nm'].value(), 0.)
            controls.preset.setCurrentIndex(1)
            for diameter in (10., 50., 100.):
                controls.fields['diameter-nm'].setValue(diameter)
                self.assertEqual(controls.fields['gap-nm'].value(), 0.)
                self.assertIn('center z = {:g} nm'.format(diameter/2), controls.defect_controls.info.text())
            controls.preset.setCurrentIndex(2)
            self.assertIn('10 nm gap', controls.preset.currentText())
            self.assertEqual(controls.fields['gap-nm'].value(), 10.)
            controls.load_physical_parameters(asdict(replace(Parameters(), gap_nm=12.)))
            self.assertEqual(controls.fields['gap-nm'].value(), 12.)
        finally:
            controls.close()
