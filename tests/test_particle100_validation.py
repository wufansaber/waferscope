"""Independent camera audit guards; no archived FEM data required."""
from dataclasses import replace

import numpy as np
import pytest

from defect_detection import calculate_detection
from detector_model import DetectorParameters
from tests.test_defect_detection import channel
from tools.particle100_camera import independent_budget, monte_carlo


@pytest.mark.parametrize('reference_frames', (1, 4))
def test_independent_camera_budget_and_sampling(reference_frames):
    p = DetectorParameters(enabled=True, pixel_pitch_um=1, pixels_x=2, pixels_y=2,
        reference_frames=reference_frames, prnu_fraction=.01, dsnu_e=.5)
    d = calculate_detection(channel(110, p), channel(30, p), 266, p)
    budget = independent_budget(d, p, 266)
    assert budget['variance_relative_max_error'] < 1e-12
    assert budget['roi_snr'] == pytest.approx(d['detection_roi_snr'])
    assert budget['template_snr'] == pytest.approx(d['detection_template_snr'])
    assert monte_carlo(d, p, frames=12000)['passed']


def test_audit_rejects_saturation_and_unknown_haze_mismatch():
    p = DetectorParameters(enabled=True, pixel_pitch_um=1, pixels_x=2, pixels_y=2)
    for case in (replace(p, full_well_e=50), replace(p, haze_mismatch_fraction=.1)):
        d = calculate_detection(channel(110, case), channel(30, case), 266, case,
                                np.random.default_rng(0))
        with pytest.raises(ValueError, match='unclipped, same-surface'):
            monte_carlo(d, case)
