"""grcwa adapter with explicit physical +Z-up fields and absolute modal powers."""
from dataclasses import asdict
import importlib.metadata
import json
import time

import numpy as np

from .geometry import layers


def compress_uniform(stack):
    result = []
    for thickness, epsilon in stack:
        uniform = not isinstance(epsilon, dict) and (epsilon is None or np.all(epsilon == epsilon.flat[0]))
        value = (1. if epsilon is None else complex(epsilon.flat[0])) if uniform else epsilon
        if uniform and result and np.ndim(result[-1][1]) == 0 and result[-1][1] == value:
            result[-1] = (result[-1][0]+thickness, value)
        else:
            result.append((thickness, value))
    return result


def outgoing_fields(model, scene):
    from grcwa.rcwa import SolveExterior
    _, backward = SolveExterior(model.a0, model.bN, model.q_list, model.phi_list,
        model.kp_list, model.thickness_list)
    n = model.nG
    hxy = model.phi_list[0]@backward
    t = model.kp_list[0]@(model.phi_list[0]@(-backward/model.omega/model.q_list[0]))
    # grcwa z points down. E is polar and H is axial under z reflection.
    electric = np.column_stack((t[n:], -t[:n], -(model.ky*hxy[:n]-model.kx*hxy[n:])/model.omega))
    kz = np.asarray(model.q_list[0][:n])
    k0 = 2*np.pi/scene.wavelength_um
    kx0, ky0 = 0., k0*np.sin(np.deg2rad(scene.incidence_deg))
    incident_phase = np.exp(-1j*(kx0*scene.lx_um/2+ky0*scene.ly_um/2+
        k0*np.cos(np.deg2rad(scene.incidence_deg))*scene.cap_top_um))
    phase = incident_phase*np.exp(1j*(model.kx*scene.lx_um/2+
        model.ky*scene.ly_um/2-kz*scene.cap_top_um))
    electric *= phase[:, None]
    return modal_result(model, scene, electric)


def modal_result(model, scene, electric):
    n = model.nG
    kz = np.asarray(model.q_list[0][:n])
    k0 = 2*np.pi/scene.wavelength_um
    valid = (abs(kz.imag) < 1e-8*k0) & (kz.real > 1e-9*k0)
    power = np.zeros(n)
    power[valid] = kz[valid].real/k0/np.cos(np.deg2rad(scene.incidence_deg))*np.sum(abs(electric[valid])**2, axis=1)
    directions = np.full((n, 3), np.nan)
    directions[valid] = np.column_stack((model.kx[valid].real, model.ky[valid].real, kz[valid].real))/k0
    return dict(orders=np.asarray(model.G), electric_reflected_z0=electric,
        reflectance_by_order=power, directions=directions, propagating=valid)


def xy_fields(model, scene, stack, offset=.02, samples=80):
    depth = scene.cap_top_um-offset
    traversed = 0.
    for layer, (thickness, _) in enumerate(stack, 1):
        if depth < traversed+thickness-1e-10:
            break
        traversed += thickness
    else:
        raise ValueError('Near plane lies outside patterned stack')
    data = model.Solve_FieldFourier(layer, float(depth-traversed))
    # The pinned implementation returns a one-element list of z samples.
    fe = np.asarray(data[0][0]).T.copy()
    fe[:, 2] *= -1
    axis_x = np.linspace(-scene.lx_um/2, scene.lx_um/2, samples, endpoint=False)
    axis_y = np.linspace(-scene.ly_um/2, scene.ly_um/2, samples, endpoint=False)
    xx, yy = np.meshgrid(axis_x, axis_y)
    k0 = 2*np.pi/scene.wavelength_um
    common = np.exp(-1j*(k0*np.sin(np.deg2rad(scene.incidence_deg))*scene.ly_um/2+
        k0*np.cos(np.deg2rad(scene.incidence_deg))*scene.cap_top_um))
    output = np.empty((xx.size, 3), complex)
    for begin in range(0, xx.size, 128):
        phase = np.exp(1j*((xx.ravel()[begin:begin+128, None]+scene.lx_um/2)*model.kx+
            (yy.ravel()[begin:begin+128, None]+scene.ly_um/2)*model.ky))
        output[begin:begin+128] = common*phase@fe
    return dict(near_x_um=axis_x, near_y_um=axis_y, near_z_um=np.asarray(offset),
        near_electric=output.reshape(samples, samples, 3))


def solve(scene, particle=False, near=False, particle_method='fft', current_z_order=0):
    import grcwa
    scene.validate()
    if particle_method not in ('fft', 'analytic'):
        raise ValueError('Unsupported particle Fourier method')
    if particle_method == 'analytic' and scene.factorization != 'laurent':
        raise ValueError('Analytic particle check requires Laurent factorization')
    if current_z_order and (not particle or particle_method != 'analytic' or
            isinstance(current_z_order, bool) or not isinstance(current_z_order, int) or not 2 <= current_z_order <= 16):
        raise ValueError('Contrast-current export requires an analytic particle and 2..16 z nodes per slice')
    if current_z_order and scene.kind == 'grating' and scene.grating_method != 'analytic_li':
        raise ValueError('Reciprocal grating current export requires the matching analytic Li background')
    grcwa.set_backend('numpy')
    start = time.perf_counter()
    if particle and particle_method == 'analytic':
        from .sphere_fourier import sphere_stack
        stack, geometry = sphere_stack(scene)
    else:
        stack, geometry = layers(scene, particle)
    stack = compress_uniform(stack)
    model = grcwa.obj(3 if scene.harmonics_x is not None else scene.modes,
        [scene.lx_um, 0], [0, scene.ly_um],
        1/scene.wavelength_um, np.deg2rad(scene.incidence_deg), np.pi/2,
        verbose=0, fmm_method='pol' if scene.factorization == 'pol' else None)
    model.Add_LayerUniform(0., 1.)
    grids = []
    analytic_layer = None
    analytic_spheres = {}
    for stack_index, (thickness, epsilon) in enumerate(stack):
        if scene.kind == 'grating' and scene.grating_method == 'analytic_li' and stack_index == len(stack)-1:
            analytic_layer = model.Layer_N
            model.Add_LayerFourier(thickness, dict(kind='rectangular_lines'))
            continue
        if isinstance(epsilon, dict):
            analytic_spheres[model.Layer_N] = epsilon['radius_um']
            model.Add_LayerFourier(thickness, epsilon)
        elif np.ndim(epsilon) == 0:
            model.Add_LayerUniform(thickness, epsilon)
        else:
            model.Add_LayerGrid(thickness, scene.grid_x, scene.grid_y)
            grids.append(epsilon.ravel())
    model.Add_LayerUniform(0., scene.eps_substrate)
    model.Init_Setup()
    if scene.harmonics_x is not None:
        from .basis import configure_rectangular
        configure_rectangular(model, scene)
    # FFT convolution needs every retained reciprocal-vector difference represented.
    extent = 2*np.ptp(np.asarray(model.G), axis=0)+1
    if grids and np.any(extent > [scene.grid_x, scene.grid_y]):
        raise ValueError('Dielectric FFT grid too small for reciprocal-vector differences: {}'.format(extent))
    print('RCWA: requested {}; actual {} modes; {} patterned layers'.format(
        scene.requested_modes, model.nG, model.Patterned_N), flush=True)
    if grids:
        model.GridLayer_geteps(np.concatenate(grids))
    if analytic_spheres:
        from .sphere_fourier import setup_layer as setup_sphere
        for index, radius in analytic_spheres.items():
            setup_sphere(model, scene, index, radius)
    if analytic_layer is not None:
        from .grating_fourier import setup_layer
        setup_layer(model, scene, analytic_layer)
    model.MakeExcitationPlanewave(p_amp=1, p_phase=0, s_amp=0, s_phase=0)
    active, take = model, np.arange(model.nG)
    if analytic_layer is not None and not particle:
        from .grating_fourier import incident_block
        active, take = incident_block(model, scene)
    result = outgoing_fields(active, scene)
    if active is not model:
        electric = np.zeros((model.nG, 3), complex)
        electric[take] = result['electric_reflected_z0']
        result = modal_result(model, scene, electric)
    reflectance, transmission = active.RT_Solve(normalize=1)
    check = abs(float(reflectance)-result['reflectance_by_order'].sum())
    if check > 1e-8 or not np.isfinite(result['electric_reflected_z0']).all():
        raise RuntimeError('Outgoing E-field/power normalization check failed: {}'.format(check))
    if near:
        result.update(xy_fields(active, scene, stack))
    if current_z_order:
        from .contrast_current import sample_current
        result.update(sample_current(model, scene, analytic_spheres, current_z_order))
    info = dict(scene=asdict(scene), geometry=geometry, actual_modes=model.nG,
        patterned_layers=model.Patterned_N, grid_layers=len(grids),
        analytic_grating_layer=analytic_layer,
        particle_method=particle_method, analytic_particle_layers=len(analytic_spheres),
        analytic_particle_branch='positive real / positive imaginary for roundoff-real lossless eigenvalues',
        current_z_order=current_z_order,
        excited_solver_modes=active.nG,
        retained_order_extent=np.max(np.abs(model.G), axis=0).tolist(),
        grating_harmonics_on_incident_block=sorted((model.G[(model.G[:, 0] == 0) &
            (model.G[:, 1] % scene.periods_y == 0), 1]//scene.periods_y).tolist()),
        elapsed_seconds=time.perf_counter()-start,
        reflectance=float(reflectance), flux_into_substrate=float(transmission),
        field_vs_library_power_absolute_difference=float(check),
        backend='grcwa', backend_version=importlib.metadata.version('grcwa'),
        coordinates='Global z up; incidence from -KY; reflected E amplitudes referenced to z=0.',
        status='EXPERIMENTAL_NOT_CONVERGED',
        material_caveat='Legacy Si n=5.57+0.387i retained for controlled comparisons, not a 266 nm calibration.',
        power_caveat='Discrete supercell order efficiencies; not a continuous DSC or isolated-particle result.')
    print(json.dumps(info, indent=2), flush=True)
    return result, info
