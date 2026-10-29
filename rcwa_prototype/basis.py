"""Explicit anisotropic Fourier truncation using grcwa's uniform modes."""
import numpy as np


def rectangular_orders(scene):
    gx, gy = np.meshgrid(np.arange(-scene.harmonics_x, scene.harmonics_x+1),
                         np.arange(-scene.harmonics_y, scene.harmonics_y+1), indexing='ij')
    orders = np.column_stack((gx.ravel(), gy.ravel()))
    radius = (orders[:, 0]/scene.lx_um)**2+(orders[:, 1]/scene.ly_um)**2
    return orders[np.argsort(radius, kind='stable')]


def configure_rectangular(model, scene):
    configure_orders(model, rectangular_orders(scene))


def configure_orders(model, orders):
    """Set an explicit basis; useful also for exact reciprocal background sectors."""
    from grcwa.kbloch import Lattice_SetKs
    from grcwa.rcwa import MakeKPMatrix, SolveLayerEigensystem_uniform
    model.G = np.asarray(orders, dtype=int)
    model.nG = len(model.G)
    model.kx, model.ky = Lattice_SetKs(model.G,
        model.omega*np.sin(model.theta)*np.cos(model.phi),
        model.omega*np.sin(model.theta)*np.sin(model.phi), model.Lk1, model.Lk2)
    model.kp_list, model.q_list, model.phi_list = [[None]*model.Layer_N for _ in range(3)]
    for index, entry in enumerate(model.id_list):
        if entry[0] != 0:
            continue
        epsilon = model.Uniform_ep_list[entry[2]]
        model.kp_list[index] = MakeKPMatrix(model.omega, 0, 1/epsilon, model.kx, model.ky)
        model.q_list[index], model.phi_list[index] = SolveLayerEigensystem_uniform(
            model.omega, model.kx, model.ky, epsilon)
