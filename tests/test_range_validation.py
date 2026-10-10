"""Exact parent recovery and physically weighted phase-sensitive aperture checks."""
from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.rough_surface import Roughness
from fullwave_prototype.shared_surface import (make_mother, save_mother, recover_synthesis_parent,
                                              surface_hash)
from fullwave_prototype.range_observables import sample_caps, cap_powers, cap_field_changes


@pytest.mark.parametrize('model', ['gaussian', 'abc'])
def test_recovered_parent_preserves_the_original_mother_exactly(tmp_path, model):
    rough = Roughness(model=model, correlation_length_um=.1)
    axis, height = make_mother(rough, dict(grid_points=32, spacing_um=.02))
    saved = save_mother(tmp_path, axis, height, rough)
    a, parent, report = recover_synthesis_parent(tmp_path)
    n = len(axis)
    np.testing.assert_array_equal(parent[n//2:n//2+n, n//2:n//2+n], height)
    np.testing.assert_array_equal(a[n//2:n//2+n], axis)
    assert report['mother_sha256'] == saved['sha256']
    assert report['parent_sha256'] == surface_hash(a, parent)
    save_mother(tmp_path, axis, height, replace(rough, seed=8))
    with pytest.raises(ValueError, match='does not reproduce'):
        recover_synthesis_parent(tmp_path)


def test_phase_change_is_detected_even_when_collected_energy_is_identical():
    class Constant:
        def __init__(self, phase):
            self.phase = phase
        def amplitude(self, directions):
            return np.full((len(directions), 2), self.phase, complex)
    reference = sample_caps(Constant(1), (.5, 0, .4), (45, 0, -45), 16)
    moved = sample_caps(Constant(1j), (.5, 0, .4), (45, 0, -45), 16)
    np.testing.assert_allclose(cap_powers(reference), cap_powers(moved))
    np.testing.assert_allclose(cap_field_changes(moved, reference), [np.sqrt(2), 0, np.sqrt(2)])
    assert cap_powers(reference)[0] == pytest.approx(4*np.pi*(1-np.sqrt(1-.5**2)))
