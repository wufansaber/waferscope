"""Opposite-ray reference illumination in a reciprocal straight-grating background."""
import numpy as np


def polarizations(direction):
    n = np.asarray(direction, float)
    if n.shape != (3,) or not np.all(np.isfinite(n)) or abs(np.linalg.norm(n)-1) > 1e-10 or n[2] <= 0:
        raise ValueError('A finite upward unit direction is required')
    transverse = np.hypot(n[0], n[1])
    s = np.array([-n[1], n[0], 0.])/transverse if transverse > 1e-12 else np.array([0., 1., 0.])
    return np.column_stack((s, np.cross(s, n)))


def model_for(scene, orders, bloch):
    import grcwa
    from .basis import configure_orders
    grcwa.set_backend('numpy')
    k0 = 2*np.pi/scene.wavelength_um
    theta = np.arcsin(np.clip(np.linalg.norm(bloch)/k0, 0., 1.))
    phi = np.arctan2(bloch[1], bloch[0])
    model = grcwa.obj(3, [scene.lx_um, 0.], [0., scene.ly_um],
                       1/scene.wavelength_um, theta, phi, verbose=0)
    model.Add_LayerUniform(0., 1.)
    if scene.kind == 'grating':
        model.Add_LayerFourier(scene.depth_um, dict(kind='rectangular_lines'))
    else:
        model.Add_LayerUniform(scene.depth_um, scene.eps_substrate)
    model.Add_LayerUniform(0., scene.eps_substrate)
    model.Init_Setup()
    configure_orders(model, orders)
    if scene.kind == 'grating':
        from .grating_fourier import setup_layer
        setup_layer(model, scene, 1, block_diagonal=False)
    return model


def spectrum(scene, direction, orders, bloch=None):
    from grcwa.rcwa import SolveExterior
    n = np.asarray(direction, float)
    pol = polarizations(n)
    k0 = 2*np.pi/scene.wavelength_um
    target = -k0*n[:2]
    model = model_for(scene, orders, target if bloch is None else bloch)
    take = np.flatnonzero(np.isclose(model.kx, target[0], atol=1e-10, rtol=1e-10) &
                          np.isclose(model.ky, target[1], atol=1e-10, rtol=1e-10))
    if len(take) != 1:
        raise ValueError('Reciprocal incident ray must occur exactly once in the reference basis')
    count = model.nG
    h = -np.cross(-n, pol.T).T[:2]
    incoming_phase = np.exp(-1j*(target[0]*scene.lx_um/2+target[1]*scene.ly_um/2))
    outgoing_phase = np.exp(1j*(model.kx*scene.lx_um/2+model.ky*scene.ly_um/2))
    # Both polarizations share one scattering matrix; keep independent RHS columns.
    a0, bN = np.zeros((2*count, 2), complex), np.zeros((2*count, 2), complex)
    a0[take[0]] = h[0]*incoming_phase
    a0[take[0]+count] = h[1]*incoming_phase
    _, back = SolveExterior(a0, bN, model.q_list, model.phi_list, model.kp_list, model.thickness_list)
    hxy = model.phi_list[0]@back
    t = model.kp_list[0]@(model.phi_list[0]@(-back/model.omega/model.q_list[0][:, None]))
    electric = np.stack((t[count:], -t[:count],
        -(model.ky[:, None]*hxy[:count]-model.kx[:, None]*hxy[count:])/model.omega), axis=1)
    reflected = electric*outgoing_phase[:, None, None]
    return dict(polarizations=pol,
        kxy=np.vstack((target, np.column_stack((model.kx.real, model.ky.real)))),
        kz=np.r_[-k0*n[2], model.q_list[0][:count]],
        electric=np.concatenate((pol[None, :, :], reflected), axis=0))
