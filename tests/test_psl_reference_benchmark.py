from dataclasses import asdict

import numpy as np
import pytest

from tools.psl_reference_benchmark import PAPER_MATERIALS, homogeneous_parameters


def test_air_reference_is_not_a_silicon_wafer_case():
    p = homogeneous_parameters(2)
    p.validate()
    assert p.wavelength_nm == 266 and p.diameter_nm == 100
    assert p.eps_substrate == 1 and p.gap_nm == 0
    assert p.polarization == 'p' and p.incidence_deg == 75
    assert p.particle_n == 1.59 and p.particle_k == 0


def test_mesh_pair_does_not_change_geometry_materials_or_illumination():
    a, b = asdict(homogeneous_parameters(2)), asdict(homogeneous_parameters(2, .05))
    assert {key for key in a if a[key] != b[key]} == {'air_mesh_um', 'substrate_mesh_um'}


def test_literature_constants_are_not_the_legacy_uv_control():
    assert PAPER_MATERIALS[442] == (1.615, 4.76, .17, 1.466)
    assert PAPER_MATERIALS[633] == (1.588, 3.88, .02, 1.457)


def test_angular_mie_oracle_has_absolute_cross_section_normalization():
    pytest.importorskip('miepython')
    pytest.importorskip('matplotlib')
    from tools.psl_reference_report import mie_dsc
    from fullwave_prototype.theory_reference import bessel_mie_cross_sections
    p = homogeneous_parameters(2)
    mu, w = np.polynomial.legendre.leggauss(24)
    mm, phi = np.meshgrid(mu, np.arange(48)*np.pi/24, indexing='ij')
    rho = np.sqrt(1-mm**2)
    rays = np.stack((rho*np.cos(phi), rho*np.sin(phi), mm), axis=-1).reshape(-1, 3)
    integral = np.repeat(w*np.pi/24, 48) @ mie_dsc(p, rays)
    np.testing.assert_allclose(integral, bessel_mie_cross_sections(p)['csca_um2'], rtol=1e-12, atol=0)
