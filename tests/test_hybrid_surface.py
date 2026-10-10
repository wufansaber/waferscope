"""Shared coordinates and coherent hybrid identities, without a native FEM solve."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import tempfile

import numpy as np
import pytest

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, height_mesh
from fullwave_prototype.shared_surface import (sample_mother, local_crop, save_mother,
    load_mother, surface_hash, make_mother)
from fullwave_prototype.hybrid_radiation import combine, triangle_source, triangle_amplitude
from fullwave_prototype.hybrid_run import optical_settings, source_signature, load_sources
from fullwave_prototype.imaging import ImagingSettings, calculate_channel, calculate_images
from fullwave_prototype.far_field import radiation_amplitude
from fullwave_prototype.planar import background_fields
from fullwave_prototype.surface_scattering import radiation_kernel


def test_crop_preserves_global_origin_orientation_and_height_offset(tmp_path):
    axis = (np.arange(32)-16)*.02
    x, y = np.meshgrid(axis, axis)
    height = .002+.001*x+.003*y
    r = Roughness(model='gaussian', width_um=.2, grid_points=9, taper_fraction=.125)
    local = local_crop(tmp_path, axis, height, r)
    assert local.detrend == 'none'
    assert local.model == 'measured_height'
    vertices, _ = height_mesh(local)
    xx, yy, zz = np.moveaxis(vertices, -1, 0)
    # The FEM mesh is X,Y indexed; the measurement file is Y,X indexed.
    np.testing.assert_allclose(zz[1:-1, 1:-1],
        (.002+.001*xx+.003*yy)[1:-1, 1:-1], atol=1e-16)
    assert zz[4, 4] == pytest.approx(.002)
    assert local.rms_nm != r.rms_nm
    with pytest.raises(ValueError):
        sample_mother(axis, height, [1], [0])


def test_mother_reproducibility_absolute_scale_and_mutation_guard(tmp_path):
    r = Roughness(model='gaussian', correlation_length_um=.05)
    plan = dict(grid_points=32, spacing_um=.02)
    axis, height = make_mother(r, plan)
    a, h = make_mother(r, plan)
    np.testing.assert_array_equal(h, height)
    np.testing.assert_array_equal(a, axis)
    np.testing.assert_allclose(make_mother(replace(r, rms_nm=4), plan)[1], 2*height)
    saved = save_mother(tmp_path, axis, height, r)
    assert saved['sha256'] == surface_hash(axis, height)
    assert load_mother(tmp_path)[2] == saved
    height[0, 0] += .001
    np.savez_compressed(tmp_path/'shared_surface.npz', axis_um=axis, height_um=height)
    with pytest.raises(ValueError, match='changed'):
        load_mother(tmp_path)


def test_complex_subtraction_keeps_local_particle_change_without_double_counting():
    rng = np.random.default_rng(7)
    fields = [rng.normal(size=(10, 2))+1j*rng.normal(size=(10, 2)) for _ in range(4)]
    whole, local, fem_with, fem_without = fields
    with_particle, without = combine(*fields)
    np.testing.assert_allclose(with_particle-without, fem_with-fem_without, atol=1e-14)
    np.testing.assert_allclose(combine(whole, local, local, local)[0], whole, atol=1e-14)
    # An explicit counterexample to summing independent intensities.
    a, _ = combine(np.array([1j]), np.array([0j]), np.array([1j]), np.array([0j]))
    assert abs(a[0])**2 == 4


def test_triangular_spm_integrates_projected_height_with_global_phase():
    p = Parameters(incidence_deg=0, polarization='s')
    w, h = .03, .002
    vertices = np.array([[-w/2, -w/2, h], [w/2, -w/2, h],
                         [w/2, w/2, h], [-w/2, w/2, h]])
    triangles = vertices[np.array([[0, 1, 2], [0, 2, 3]])]
    positions, weights = triangle_source(triangles, p)
    q = np.array([[0., 0., 1.], [.3, .4, np.sqrt(.75)]])
    fields = background_fields(np.zeros(3), p)[0]
    integral = h*w*w*np.sinc(p.k0*q[:, 0]*w/(2*np.pi))*np.sinc(p.k0*q[:, 1]*w/(2*np.pi))
    exact = (radiation_kernel(q, p)@fields)*integral[:, None]
    np.testing.assert_allclose(triangle_amplitude(positions, weights, q, p), exact, rtol=1e-12, atol=1e-20)
    shifted = triangles+np.array([.12, -.07, 0])
    sp, sw = triangle_source(shifted, p)
    expected = exact*np.exp(-1j*p.k0*(q[:, :2]@np.array([.12, -.07])))[:, None]
    np.testing.assert_allclose(triangle_amplitude(sp, sw, q, p), expected, rtol=1e-12, atol=1e-20)


def test_custom_imaging_matches_original_source_and_rejects_local_cache():
    p = Parameters(substrate_n=1, substrate_k=0)
    settings = ImagingSettings(pupil_samples=32, padding=2)
    data = {key+name: value for name, weight in [('with', 1), ('without', .2j)]
            for key, value in [('source_points_', np.zeros((1, 3))),
                              ('source_polarization_weights_', np.array([[0, weight, 0]], complex))]}
    def provider(q, state):
        return radiation_amplitude(data['source_points_'+state],
                                   data['source_polarization_weights_'+state], q, p)
    base, _ = calculate_channel(data, p, settings, 0)
    custom, _ = calculate_channel({}, p, settings, 0, amplitude_provider=provider)
    for key in base:
        np.testing.assert_allclose(base[key], custom[key], rtol=1e-14, atol=1e-20)
    with pytest.raises(ValueError, match='cache'):
        calculate_channel({}, p, settings, 0, cache=object(), amplitude_provider=provider)
    with pytest.raises(ValueError, match='cache'):
        calculate_images({}, p, settings, cache_directory='unused', amplitude_provider=provider)


def test_hybrid_reimage_rejects_changed_beam_shape_and_accepts_power_only():
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    from fullwave_prototype.imaging_options import imaging_settings
    args = build_parser().parse_args(['--beam-profile', 'gaussian'])
    p, _ = scattering_settings(args)
    settings = imaging_settings(args, parameters=p)
    result = optical_settings(asdict(settings), replace(p, beam_order=41))
    result.validate()
    doubled = replace(settings, illumination_power_w=2*settings.illumination_power_w,
                      irradiance_w_um2=2*settings.irradiance_w_um2)
    assert optical_settings(asdict(doubled), p).irradiance_w_um2 == pytest.approx(2*settings.irradiance_w_um2)
    with pytest.raises(ValueError, match='same finite coherent beam'):
        optical_settings(asdict(settings), replace(p, beam_length_um=3))


def test_source_signature_covers_beam_psd_plan_and_mother(tmp_path):
    from fullwave_prototype.hybrid_run import MODEL
    metadata = dict(model=MODEL, source_format_version=2, parameters={}, roughness={},
        stationary_roughness={}, plan={}, source_sha256='source', mother_surface=dict(sha256='mother'))
    metadata['physical_signature'] = source_signature(metadata)
    metadata['parameters']['beam_length_um'] = 3
    (tmp_path/'hybrid_source_report.json').write_text(json.dumps(metadata), encoding='utf-8')
    with pytest.raises(ValueError, match='metadata changed'):
        load_sources(tmp_path)


def test_reimage_verifies_original_mother_and_complex_source_files(tmp_path):
    import hashlib
    from fullwave_prototype.hybrid_run import MODEL
    axis = (np.arange(16)-8)*.02
    height = np.ones((16, 16))*.001
    mother = save_mother(tmp_path, axis, height, Roughness())
    path = tmp_path/'hybrid_fields.npz'
    np.savez_compressed(path, field=np.ones(2, complex))
    metadata = dict(model=MODEL, source_format_version=2, parameters={}, roughness={},
        stationary_roughness={}, plan={}, mother_surface=mother,
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    metadata['physical_signature'] = source_signature(metadata)
    (tmp_path/'hybrid_source_report.json').write_text(json.dumps(metadata), encoding='utf-8')
    assert load_sources(tmp_path)[0]['field'].shape == (2,)
    original = path.read_bytes()
    np.savez_compressed(path, field=np.zeros(2, complex))
    with pytest.raises(ValueError, match='complex source changed'):
        load_sources(tmp_path)
    path.write_bytes(original)
    # Even a self-consistent replacement of both mother files cannot reuse the old FEM solve.
    save_mother(tmp_path, axis, height*2, Roughness())
    with pytest.raises(ValueError, match='no longer matches'):
        load_sources(tmp_path)


def test_fourier_grid_refinement_keeps_the_same_realization_and_global_origin():
    from fullwave_prototype.hybrid_surface_audit import refined_mother
    axis = (np.arange(32)-16)*.02
    x, y = np.meshgrid(axis, axis)
    height = .002+.001*np.cos(2*np.pi*x/.64)+.003*np.sin(4*np.pi*y/.64)
    refined_axis, refined = refined_mother(axis, height)
    np.testing.assert_allclose(refined_axis[::2], axis, atol=1e-16)
    np.testing.assert_allclose(refined[::2, ::2], height, atol=1e-17)
    xx, yy = np.meshgrid(refined_axis, refined_axis)
    np.testing.assert_allclose(refined,
        .002+.001*np.cos(2*np.pi*xx/.64)+.003*np.sin(4*np.pi*yy/.64), atol=1e-17)
