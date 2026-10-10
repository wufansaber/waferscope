import numpy as np
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.surface_planning import coverage_plan
from fullwave_prototype.energy_balance import sphere_points


def test_bisection_tracks_a_small_bracket_without_claiming_scattering_accuracy(monkeypatch):
    p = Parameters(beam_profile='gaussian', beam_length_um=2, beam_width_um=2)

    def evaluate(parameters, width, spacing, taper):
        return parameters, dict(width_um=width, planning_power_coverage=1-np.exp(-width))

    monkeypatch.setattr('fullwave_prototype.surface_planning.coverage_at_width', evaluate)
    _, plan = coverage_plan(p, .1, .99, .02, 2048)
    low, high = plan['coverage_bracket_um']
    assert high-low <= .01*high
    assert low < np.log(100) <= high
    assert not plan['scattering_convergence_verified']


def test_energy_surface_quadrature_is_closed_and_has_correct_area():
    normals, areas = sphere_points(.2, 16)
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1, atol=1e-14)
    assert abs(areas.sum()-4*np.pi*.2**2) < 1e-14
    np.testing.assert_allclose(np.sum(normals*areas[:, None], axis=0), 0, atol=1e-14)
