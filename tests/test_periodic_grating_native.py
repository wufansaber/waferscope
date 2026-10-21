import numpy as np
import pytest

ng = pytest.importorskip('ngsolve')

from fullwave_prototype.periodic_grating import GratingInputs, solve, extract
from fullwave_prototype.planar import fresnel


def test_periodic_flat_wafer_is_analytic_fresnel_null_scatter():
    p = GratingInputs(ridge_width_um=.5)
    mesh, field, h, dh, epsinv, info = solve(p, order=2, factor=1.4, pml=.45)
    arrays, checks = extract(p, mesh, field, h, dh, epsinv, samples=1024)
    rp = fresnel(np.cos(np.deg2rad(75)), p.epsilon)[1]
    assert np.linalg.norm(field.vec.FV().NumPy()) == 0
    assert checks['reflectance'] == pytest.approx(abs(rp)**2, rel=1e-12)
    assert abs(checks['energy_closure_error']) < 1e-6


def test_periodic_pattern_energy_and_two_plane_checks():
    p = GratingInputs()
    mesh, field, h, dh, epsinv, info = solve(p, order=3, pml=.45)
    arrays, checks = extract(p, mesh, field, h, dh, epsinv, samples=2048)
    assert .2 < checks['reflectance'] < .4
    assert abs(checks['energy_closure_error']) < .002
    assert checks['modal_plane_relative_l2'] < .001
