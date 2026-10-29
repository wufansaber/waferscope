from dataclasses import replace
import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from s4_prototype.finite_patch import FinitePatch, volume_quadrature
from s4_prototype.lua_scene import script


def scene():
    return Scene(lx_um=1.5,periods_y=4,harmonics_x=6,harmonics_y=8,gap_um=0.)


def test_finite_sources_have_physical_volumes_and_centroids():
    s=scene(); patch=FinitePatch()
    q,w,particle=volume_quadrature(s,patch,lateral=4,vertical=4)
    assert np.all(w>0)
    assert w[~particle].sum()==pytest.approx(2*.8*.25*.1)
    assert w[particle].sum()==pytest.approx(np.pi*.1**3/6*(1+1/(2*s.slices**2)))
    np.testing.assert_allclose(np.sum(q*w[:,None],axis=0)[:2],0,atol=1e-16)
    assert np.sum(q[~particle,2]*w[~particle])/w[~particle].sum()==pytest.approx(-s.depth_um/2)
    assert np.sum(q[particle,2]*w[particle])/w[particle].sum()==pytest.approx(s.gap_um+s.diameter_um/2)


def test_native_script_has_flat_buffers_and_real_trench_ends():
    s=scene(); source,_=script(s,particle=True,finite_patch=FinitePatch())
    assert "S:AddLayer('grating',0.10000000000000001,'Wafer')" in source
    assert source.count("S:SetLayerPatternRectangle('grating','Air'")==2
    assert 'UseAnalyticStripeRule' not in source
    assert "S:SetLayerPatternCircle('particle_" in source


@pytest.mark.parametrize('change', [dict(length_um=1.5),dict(trench_count=3),dict(trench_count=True),dict(trench_count=6)])
def test_bad_or_unbuffered_patch_is_rejected(change):
    with pytest.raises(ValueError): replace(FinitePatch(),**change).validate(scene())


def test_finite_edges_cannot_use_invariant_stripe_shortcut():
    with pytest.raises(ValueError,match='full 2D'):
        script(scene(),finite_patch=FinitePatch(),stripe_rule=True)
