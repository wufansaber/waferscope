import numpy as np
import pytest

from rcwa_prototype.current_overlap import compare_source_fields


def source(fields, wavevectors):
    return dict(source_electric_fourier=np.asarray(fields, complex),
        source_kxy_per_um=np.asarray(wavevectors, float),
        source_z_um=np.array([.02, .08]), source_weight_um=np.array([.05, .05]),
        source_radius_um=np.array([.03, .04]))


def test_absolute_field_amplitude_and_phase_remain_visible():
    a = source(np.ones((2, 1, 3)), [[0, 0]])
    b = source(2*np.ones((2, 1, 3)), [[0, 0]])
    assert compare_source_fields(a, b, 1.5)['internal_e_relative_l2'] == pytest.approx(1.)
    b = source(1j*np.ones((2, 1, 3)), [[0, 0]])
    result = compare_source_fields(a, b, 1.5)
    assert result['internal_e_relative_l2'] == pytest.approx(np.sqrt(2))
    assert result['contrast_moment_relative_change'] == pytest.approx(np.sqrt(2))


def test_different_fourier_bases_can_represent_the_same_field():
    a = source(np.ones((2, 1, 3)), [[0, 0]])
    fields = np.zeros((2, 2, 3), complex)
    fields[:, 0] = 1
    b = source(fields, [[0, 0], [5, 7]])
    result = compare_source_fields(a, b, 1.5)
    assert result['internal_e_relative_l2'] < 1e-7
    assert result['contrast_moment_relative_change'] < 1e-13


def test_disk_overlap_matches_independent_polar_integration():
    rng = np.random.default_rng(31)
    a = source(rng.normal(size=(2, 3, 3))+1j*rng.normal(size=(2, 3, 3)),
               rng.normal(size=(3, 2))*20)
    b = source(rng.normal(size=(2, 4, 3))+1j*rng.normal(size=(2, 4, 3)),
               rng.normal(size=(4, 2))*20)
    norm, delta = 0., 0.
    nodes, weights = np.polynomial.legendre.leggauss(20)
    phi = np.arange(80)*2*np.pi/80
    for i, radius in enumerate(a['source_radius_um']):
        r = (nodes+1)*radius/2
        rr, pp = np.meshgrid(r, phi, indexing='ij')
        points = np.column_stack((rr.ravel()*np.cos(pp.ravel()), rr.ravel()*np.sin(pp.ravel())))
        w = np.repeat(weights*radius/2*r, 80)*2*np.pi/80*a['source_weight_um'][i]
        ea = np.exp(1j*points@a['source_kxy_per_um'].T)@a['source_electric_fourier'][i]
        eb = np.exp(1j*points@b['source_kxy_per_um'].T)@b['source_electric_fourier'][i]
        norm += w@np.sum(abs(ea)**2, axis=1)
        delta += w@np.sum(abs(eb-ea)**2, axis=1)
    result = compare_source_fields(a, b, 1.5)
    assert result['first_internal_e_volume_norm_um3'] == pytest.approx(norm, rel=1e-12)
    assert result['internal_e_relative_l2'] == pytest.approx(np.sqrt(delta/norm), rel=1e-12)


def test_different_particle_sampling_is_rejected():
    a = source(np.ones((2, 1, 3)), [[0, 0]])
    b = source(np.ones((2, 1, 3)), [[0, 0]])
    b['source_z_um'][0] += .001
    with pytest.raises(ValueError, match='Matching physical particle'):
        compare_source_fields(a, b, 1.5)


def test_nonfinite_data_cannot_produce_a_valid_report():
    a = source(np.ones((2, 1, 3)), [[0, 0]])
    a['source_electric_fourier'][0, 0, 0] = np.nan
    with pytest.raises(ValueError, match='Finite Fourier fields'):
        compare_source_fields(a, a, 1.5)
