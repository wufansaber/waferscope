from copy import deepcopy
from dataclasses import asdict

import pytest

from rcwa_prototype.geometry import Scene
from rcwa_prototype.supercell import compare_cells


def cells():
    scene = Scene(harmonics_x=2, harmonics_y=40)
    ch = dict(na=.6, tilt_deg=45., collected_orders=4,
              difference_field_efficiency=.004, actual_power_difference=.003)
    coarse = dict(info=dict(backend='S4', scene=asdict(scene), executable_sha256='test'),
                  summary=dict(channels=[ch], hemisphere_difference_field_efficiency=.008))
    fine = deepcopy(coarse)
    fine['info']['scene'].update(lx_um=1., harmonics_x=4)
    fine['summary']['hemisphere_difference_field_efficiency'] = .004
    fine['summary']['channels'][0].update(difference_field_efficiency=.002,
                                         actual_power_difference=.0015, collected_orders=8)
    return coarse, fine


def test_cell_comparison_uses_cross_section_not_raw_efficiency():
    result = compare_cells(*cells())
    assert result['channels'][0]['relative_cross_section_change'] == pytest.approx(0.)
    assert result['field_difference_computed'] is False
    assert result['physical_frequency_halfwidth_cycles_per_um'] == [4., 40.]


@pytest.mark.parametrize('change', [dict(harmonics_x=2), dict(slices=12), dict(substrate_n=4.)])
def test_cell_comparison_rejects_unmatched_cutoff_or_inputs(change):
    coarse, fine = cells()
    fine['info']['scene'].update(change)
    with pytest.raises(ValueError):
        compare_cells(coarse, fine)


def test_cell_comparison_requires_same_solver_and_apertures():
    coarse, fine = cells()
    fine['info']['backend'] = 'grcwa'
    with pytest.raises(ValueError, match='solver'):
        compare_cells(coarse, fine)
    coarse, fine = cells()
    fine['summary']['channels'][0]['na'] = .5
    with pytest.raises(ValueError, match='apertures'):
        compare_cells(coarse, fine)


def test_report_rejects_mixed_sample_reference_material_rules(monkeypatch):
    from tools import rcwa_supercell_report as report
    info = dict(backend='grcwa', scene=asdict(Scene()), particle=True, particle_method='analytic')
    reference_info = dict(info, particle=False, particle_method='fft')
    entries = iter([(info, {}), (reference_info, {})])
    monkeypatch.setattr(report, 'checked_load', lambda path: next(entries))
    with pytest.raises(ValueError, match='numerical settings'):
        report.pair(dict(sample='unused', reference='unused'))
