"""Real small SCATMECH calculations exercise the UI exporters."""
import numpy as np
from matplotlib.figure import Figure
from scatmech_ui.parameters import (SimulationParameters, ObjectiveCollectionParameters,
                                    ParticleSizeSweepParameters)
from scatmech_ui.size_sweep import calculate_particle_size_sweep, export_particle_size_sweep
from scatmech_ui.imaging import ImagingParameters, calculate_imaging, export_imaging
from result_archive import load_result_archive


def test_sweep_exports_each_size_native_field(tmp_path):
    result = calculate_particle_size_sweep(SimulationParameters(angular_step_deg=15.),
        ObjectiveCollectionParameters(grid_size=51), ParticleSizeSweepParameters(30., 40., 10.))
    export_particle_size_sweep(result, tmp_path, Figure())
    saved, source = load_result_archive(tmp_path/'particle_size_sweep_data.npz')
    assert source is None
    assert len(saved['far_field_samples']) == 2
    for original, restored in zip(result['far_field_samples'], saved['far_field_samples']):
        np.testing.assert_array_equal(original['mueller_raw'], restored['mueller_raw'])
    np.testing.assert_array_equal(saved['sampled_dsc_um2_sr'], result['sampled_dsc_um2_sr'])
    assert len(saved['collection_sampling']['channel_masks']) == 3


def test_imaging_exports_raw_phases_and_legacy_keys(tmp_path):
    result = calculate_imaging(SimulationParameters(), ObjectiveCollectionParameters(),
                               ImagingParameters(pupil_samples=32, padding=2))
    export_imaging(result, tmp_path, Figure())
    restored, _ = load_result_archive(tmp_path/'coherent_imaging.npz')
    assert len(restored['channels']) == 3
    with np.load(tmp_path/'coherent_imaging.npz', allow_pickle=False) as saved:
        for index in range(3):
            original = result['channels'][index]['image_e_modes']
            np.testing.assert_array_equal(restored['channels'][index]['image_e_modes'], original)
            np.testing.assert_array_equal(saved[f'ch{index+1}_image_e_modes'], original)
