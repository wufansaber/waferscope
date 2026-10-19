from dataclasses import asdict, replace
import json
import numpy as np
import pytest
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness
from fullwave_prototype.imaging import ImagingSettings
from fullwave_prototype.finite_illumination import power_area
from fullwave_prototype.whole_haze import run
from fullwave_prototype.whole_haze_resume import load_resume


def test_extension_equals_one_fixed_count_job_and_preserves_source(tmp_path, monkeypatch):
    monkeypatch.setattr('fullwave_prototype.whole_haze_plotting.plot_results', lambda *args: None)
    p = Parameters(wavelength_nm=266, beam_profile='gaussian', beam_length_um=2, beam_width_um=2, beam_order=31)
    r = Roughness(model='gaussian', rms_nm=1, correlation_length_um=.1)
    area = power_area(p)
    settings = ImagingSettings(illumination_model='gaussian', illumination_power_w=1,
        illumination_length_um=2, illumination_width_um=2, illumination_incidence_deg=75,
        irradiance_w_um2=1/area, illumination_normalization_area_um2=area, pupil_samples=32, padding=2)
    inputs = dict(parameters=asdict(p), roughness=asdict(r), imaging_settings=asdict(settings), samples=2)
    first, extended, full = (tmp_path/name for name in ('first', 'extended', 'full'))
    run(inputs, first)
    original = (first/'whole_haze_fields.npz').read_bytes()
    inputs['samples'] = 4
    result = run(inputs, extended, first)
    direct = run(inputs, full)
    np.testing.assert_allclose(result['statistics']['sample_observables_um2'],
                               direct['statistics']['sample_observables_um2'], rtol=1e-12)
    with np.load(extended/'whole_haze_fields.npz') as a, np.load(full/'whole_haze_fields.npz') as b:
        for key in a.files:
            np.testing.assert_allclose(a[key], b[key], rtol=1e-12, atol=1e-20)
    assert (first/'whole_haze_fields.npz').read_bytes() == original
    with pytest.raises(ValueError, match='identical scattering'):
        load_resume(first, replace(p, wavelength_nm=355), r, settings, 4)
    (first/'whole_haze_fields.npz').write_bytes(original+b'changed')
    with pytest.raises(ValueError, match='fields changed'):
        load_resume(first, p, r, settings, 4)
