from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.periodic_grating import GratingInputs, rayleigh_orders, modal_power
from fullwave_prototype.reciprocity_surface import box_quadrature, amplitude
from fullwave_prototype.parameters import Parameters


def test_periodic_diffraction_angles_and_power_obliquity():
    m, tangent, kz = rayleigh_orders(GratingInputs())
    mask = kz.imag == 0
    np.testing.assert_array_equal(m[mask], [-3, -2, -1, 0])
    assert tangent[m == 0][0] == pytest.approx(np.sin(np.deg2rad(75)))
    assert modal_power(np.array([1., 1.]), np.array([1., .5]), 60).tolist() == pytest.approx([2., 1.])


@pytest.mark.parametrize('change', [dict(period_um=-1), dict(ridge_width_um=.6),
    dict(depth_um=float('nan')), dict(substrate_k=-1), dict(incidence_deg=90)])
def test_grating_rejects_invalid_inputs(change):
    with pytest.raises(ValueError):
        replace(GratingInputs(), **change).validate()


def test_closed_box_orientation_and_surface_area():
    points, normals, areas = box_quadrature([[-.4, .4], [-.3, .3], [-.2, .5]], order=8)
    assert areas.sum() == pytest.approx(2*(.8*.6+.8*.7+.6*.7))
    np.testing.assert_allclose(np.sum(normals*areas[:, None], axis=0), 0, atol=1e-15)
    assert np.sum(points*normals*areas[:, None]) == pytest.approx(3*.8*.6*.7)


def test_reciprocity_extraction_against_exact_displaced_electric_dipole():
    p = Parameters(wavelength_nm=266, substrate_n=1, substrate_k=0)
    points, normals, areas = box_quadrature([[-.3, .4], [-.35, .3], [-.25, .4]], order=48)
    origin = np.array((.03, -.04, .02))
    dipole = np.array((.5+.1j, -.3j, 1.))
    r = points-origin
    radius = np.linalg.norm(r, axis=1)
    n = r/radius[:, None]
    ndotp = n@dipole
    transverse = dipole-n*ndotp[:, None]
    longitudinal = 3*n*ndotp[:, None]-dipole
    phase = np.exp(1j*p.k0*radius)/(4*np.pi*radius)
    e = phase[:, None]*(p.k0**2*transverse+
        (1/radius**2-1j*p.k0/radius)[:, None]*longitudinal)
    h = (phase*(p.k0**2+1j*p.k0/radius))[:, None]*np.cross(n, dipole)
    _, tangent, kz = rayleigh_orders(GratingInputs())
    valid = kz.imag == 0
    directions = np.column_stack((np.zeros(valid.sum()), tangent[valid], kz[valid].real))
    actual = amplitude(points, normals, areas, e, h, directions, p)
    s = np.cross(np.array((0., 0., 1.)), directions)
    s /= np.linalg.norm(s, axis=1)[:, None]
    polarization = np.cross(directions, s)
    expected = (p.k0**2/(4*np.pi)*np.column_stack((s@dipole, polarization@dipole))*
        np.exp(-1j*p.k0*(directions@origin))[:, None])
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-9)


def test_reciprocity_closed_surface_is_zero_for_source_free_planar_background():
    from fullwave_prototype.planar import background_fields
    p = Parameters(wavelength_nm=266)
    points, normals, areas = box_quadrature([[-.4, .4], [-.4, .4], [-.14, .23]], order=48)
    e, h = background_fields(points, p)
    directions = np.array([[.3, .2, np.sqrt(.87)], [-.4, .5, np.sqrt(.59)]])
    result = amplitude(points, normals, areas, e, h, directions, p)
    np.testing.assert_allclose(result, 0, atol=1e-8)
