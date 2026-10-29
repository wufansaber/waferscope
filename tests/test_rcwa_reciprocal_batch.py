from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip('grcwa')

from rcwa_prototype.geometry import Scene
from rcwa_prototype.reciprocal_background import model_for, polarizations, spectrum


@pytest.mark.parametrize('kind', ['flat', 'grating'])
@pytest.mark.parametrize('direction', [[.3, -.4, np.sqrt(.75)], [0., 0., 1.]])
def test_batched_polarizations_match_original_separate_exterior_solves(kind, direction):
    from grcwa.rcwa import SolveExterior
    scene = Scene(kind=kind, periods_y=1, grating_method='analytic_li')
    orders = np.column_stack((np.zeros(9, int), np.arange(-4, 5)))
    direction = np.asarray(direction)
    pol = polarizations(direction)
    k0 = 2*np.pi/scene.wavelength_um
    target = -k0*direction[:2]
    model = model_for(scene, orders, target)
    index, = np.flatnonzero(np.isclose(model.kx, target[0]) & np.isclose(model.ky, target[1]))
    incoming_phase = np.exp(-1j*(target[0]*scene.lx_um/2+target[1]*scene.ly_um/2))
    outgoing_phase = np.exp(1j*(model.kx*scene.lx_um/2+model.ky*scene.ly_um/2))
    h = -np.cross(-direction, pol.T).T[:2]
    columns = []
    for i in range(2):
        a0, bN = np.zeros(2*model.nG, complex), np.zeros(2*model.nG, complex)
        a0[index], a0[index+model.nG] = h[:, i]*incoming_phase
        _, back = SolveExterior(a0, bN, model.q_list, model.phi_list, model.kp_list, model.thickness_list)
        hxy = model.phi_list[0]@back
        t = model.kp_list[0]@(model.phi_list[0]@(-back/model.omega/model.q_list[0]))
        fields = np.column_stack((t[model.nG:], -t[:model.nG],
            -(model.ky*hxy[:model.nG]-model.kx*hxy[model.nG:])/model.omega))
        columns.append(fields*outgoing_phase[:, None])
    expected = np.concatenate((pol[None], np.stack(columns, axis=-1)), axis=0)
    result = spectrum(scene, direction, orders)
    np.testing.assert_allclose(result['electric'], expected, rtol=2e-12, atol=2e-13)
    enlarged = spectrum(replace(scene, lx_um=.8333333333333334), direction, orders)
    np.testing.assert_allclose(enlarged['electric'], result['electric'], rtol=2e-12, atol=2e-13)
