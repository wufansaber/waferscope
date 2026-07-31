import unittest

import numpy as np

from scatmech_ui.parameters import (
    CoatingLayerParameters,
    SimulationParameters,
)
from scatmech_ui.surface_model import (
    SCATMECH_IMPORT_ERROR,
    abc_psd_a_um4,
    build_scattering_model,
)


class SurfaceParameterTests(unittest.TestCase):
    def test_rectangular_footprint_area(self):
        parameters = SimulationParameters(
            roughness_enabled=True,
            illumination_length_um=20.0,
            illumination_width_um=8.0,
        )

        parameters.validate()

        self.assertEqual(parameters.illuminated_area_um2, 160.0)

    def test_abc_psd_parameter_integrates_to_requested_rms(self):
        rms_nm = 0.2
        length_um = 0.35
        exponent = 3.4

        a_um4 = abc_psd_a_um4(rms_nm, length_um, exponent)
        integrated_variance_um2 = (
            2.0
            * np.pi
            * a_um4
            / (length_um**2 * (exponent - 2.0))
        )

        self.assertAlmostEqual(
            integrated_variance_um2,
            (rms_nm / 1000.0) ** 2,
        )

    def test_invalid_layer_and_roughness_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            SimulationParameters(
                coating_layers=(
                    CoatingLayerParameters(thickness_nm=0.0),
                )
            ).validate()
        with self.assertRaises(ValueError):
            SimulationParameters(
                roughness_enabled=True,
                roughness_exponent=2.0,
            ).validate()


@unittest.skipIf(
    SCATMECH_IMPORT_ERROR is not None,
    "pySCATMECH is unavailable",
)
class ScatmechSurfaceModelTests(unittest.TestCase):
    def test_particle_only_model_preserves_local_dsc_path(self):
        bundle = build_scattering_model(SimulationParameters())

        self.assertEqual(bundle.model.getModelName(), "Bobbert_Vlieger_BRDF_Model")
        self.assertEqual(bundle.output_kind, "dsc")
        self.assertEqual(bundle.stack_layer_count, 0)

    def test_multilayer_roughness_uses_two_source_brdf(self):
        parameters = SimulationParameters(
            coating_layers=(
                CoatingLayerParameters(
                    material="Fused silica (SiO2)",
                    n=1.476,
                    thickness_nm=80.0,
                ),
                CoatingLayerParameters(
                    material="Alumina (Al2O3)",
                    n=1.79,
                    thickness_nm=35.0,
                ),
            ),
            roughness_enabled=True,
            roughness_interface_mode="top_interface",
            illumination_length_um=20.0,
            illumination_width_um=8.0,
        )

        bundle = build_scattering_model(parameters)
        model_parameters = bundle.model.getParameters()

        self.assertEqual(bundle.model.getModelName(), "Two_Source_BRDF_Model")
        self.assertEqual(bundle.output_kind, "brdf_to_dsc")
        self.assertEqual(bundle.effective_area_um2, 160.0)
        self.assertEqual(bundle.stack_layer_count, 2)
        self.assertEqual(
            model_parameters["source1"],
            "Bobbert_Vlieger_BRDF_Model",
        )
        self.assertEqual(
            model_parameters["source2"],
            "Roughness_Stack_BRDF_Model",
        )
        self.assertEqual(model_parameters["source2.this_layer"], "2")


if __name__ == "__main__":
    unittest.main()
