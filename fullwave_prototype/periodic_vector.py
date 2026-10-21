"""Periodic Maxwell bridge using the existing 3D H(curl) solver and sources."""
import numpy as np


def build_cell(g, p, line_period=.08, top=.4, bottom=.3, corner_mesh=None):
    import ngsolve as ng
    from netgen.occ import Box, Glue, OCCGeometry, Translation
    scale = 1000.
    xs = np.array([-line_period/2, 0, line_period/2])
    ys = np.array([-g.period_um/2, -g.ridge_width_um/2, g.ridge_width_um/2, g.period_um/2])
    zs = np.array([-bottom-p.pml_thickness_um, -bottom, -g.depth_um, 0, top,
                   top+p.pml_thickness_um])
    regions = []
    for ix in range(2):
        for iy in range(3):
            for iz in range(5):
                name = ('pml_substrate' if iz == 0 else 'substrate' if iz == 1 else
                        ('substrate' if iy == 1 else 'pit') if iz == 2 else
                        'air' if iz == 3 else 'pml_air')
                box = Box(tuple(scale*np.array((xs[ix], ys[iy], zs[iz]))),
                          tuple(scale*np.array((xs[ix+1], ys[iy+1], zs[iz+1]))))
                box.mat(name)
                size = dict(substrate=p.substrate_mesh_um, pit=min(.035, p.air_mesh_um),
                    air=p.air_mesh_um, pml_substrate=p.substrate_mesh_um*p.pml_mesh_scale,
                    pml_air=p.air_mesh_um*p.pml_mesh_scale)[name]
                box.maxh = size*scale
                regions.append(box)
    shape = Glue(regions)
    if corner_mesh is not None:
        if not np.isfinite(corner_mesh) or corner_mesh <= 0:
            raise ValueError('Positive corner mesh size required')
        for edge in shape.edges:
            y, z = edge.center.y/scale, edge.center.z/scale
            if abs(abs(y)-g.ridge_width_um/2) < 1e-9 and min(abs(z), abs(z+g.depth_um)) < 1e-9:
                edge.maxh = corner_mesh*scale
    shape.faces.name = 'interface'
    tolerance = 1e-6
    for axis, width, name in ((0, line_period, 'px'), (1, g.period_um, 'py')):
        left = [f for f in shape.faces if abs(tuple(f.center)[axis]+width*scale/2) < tolerance]
        right = [f for f in shape.faces if abs(tuple(f.center)[axis]-width*scale/2) < tolerance]
        for face in left:
            target = np.array(tuple(face.center))
            target[axis] += width*scale
            matches = [f for f in right if np.linalg.norm(np.array(tuple(f.center))-target) < tolerance]
            if len(matches) != 1:
                raise ValueError('Ambiguous periodic face pairing')
            vector = np.zeros(3)
            vector[axis] = width*scale
            face.name, matches[0].name = name+'_left', name+'_right'
            face.Identify(matches[0], name, trafo=Translation(tuple(vector)))
    for face in shape.faces:
        if abs(face.center.z-zs[0]*scale) < tolerance or abs(face.center.z-zs[-1]*scale) < tolerance:
            face.name = 'outer'
    native = OCCGeometry(shape).GenerateMesh(maxh=p.air_mesh_um*scale)
    native.Scale(1/scale)
    native.SetGeometry(OCCGeometry(shape.Scale((0, 0, 0), 1/scale)))
    mesh = ng.Mesh(native)
    transform = ng.pml.Cartesian(mins=(-line_period/2, -g.period_um/2, -bottom),
        maxs=(line_period/2, g.period_um/2, top), alpha=1j*p.pml_strength)
    for name in ('pml_air', 'pml_substrate'):
        mesh.SetPML(transform, definedon=name)
    # Netgen lists master/minion point pairs; derive each identification's phase.
    shifts = {}
    beta = p.k0*np.sin(np.deg2rad(p.incidence_deg))
    for first, second, number in native.GetIdentifications():
        delta = np.array(native[second].p)-np.array(native[first].p)
        phase = complex(np.exp(1j*beta*delta[1]))
        if number in shifts and not np.isclose(phase, shifts[number], rtol=0, atol=1e-9):
            raise ValueError('Inconsistent periodic phase within an identification')
        shifts[number] = phase
    if set(shifts) != set(range(1, len(shifts)+1)):
        raise ValueError('Non-contiguous periodic identifications')
    from .field_space import make_space
    space = ng.Periodic(make_space(mesh, p), phase=[shifts[i] for i in sorted(shifts)])
    return mesh, space


def extract_modes(g, p, mesh, field, line_period=.08, heights=(.15, .25), samples=2048):
    import ngsolve as ng
    from .solver import evaluate_fields
    from .periodic_grating import rayleigh_orders, modal_power
    from .planar import fresnel
    m, tangent, kz = rayleigh_orders(g)
    y = (np.arange(samples)+.5)*g.period_um/samples-g.period_um/2
    beta = p.k0*np.sin(np.deg2rad(p.incidence_deg))
    fourier = np.exp(-1j*p.k0*tangent[:, None]*y[None, :])
    amplitudes, invariance = [], []
    for z in heights:
        # Average along two inequivalent line coordinates as an invariance check.
        rows = []
        for x in (-line_period*.19, line_period*.17):
            points = np.column_stack((np.full(samples, x), y, np.full(samples, z)))
            _, magnetic, _ = evaluate_fields(mesh, field, points, p)
            h = magnetic[:, 0]
            # Remove the incoming plane wave before modal extraction.
            h = h-np.exp(1j*beta*y-1j*p.k0*np.cos(np.deg2rad(p.incidence_deg))*z)
            rows.append((fourier@h/samples)*np.exp(-1j*p.k0*kz*z))
        amplitudes.append(np.mean(rows, axis=0))
        valid = abs(kz.imag) < 1e-12
        invariance.append(float(np.linalg.norm((rows[0]-rows[1])[valid])/
            max(np.linalg.norm(rows[1][valid]), 1e-300)))
    valid = abs(kz.imag) < 1e-12
    efficiency = modal_power(amplitudes[0], kz, p.incidence_deg)
    z = np.linspace(.08, .3, 13)
    left = np.column_stack((np.zeros(len(z)), np.full(len(z), -g.period_um/2), z))
    right = left.copy()
    right[:, 1] = g.period_um/2
    e_left = evaluate_fields(mesh, field, left, p)[0]
    e_right = evaluate_fields(mesh, field, right, p)[0]
    tangential = [0, 2]
    trace_error = float(np.linalg.norm((e_right-e_left*np.exp(1j*beta*g.period_um))[:, tangential])/
        np.linalg.norm(e_right[:, tangential]))
    return dict(orders=m, tangent=tangent, kz=kz, amplitude=np.asarray(amplitudes),
        efficiency=efficiency, propagating=valid), dict(bloch_trace_relative_l2=trace_error,
        line_invariance_relative_l2=max(invariance),
        modal_plane_relative_l2=float(np.linalg.norm((amplitudes[0]-amplitudes[1])[valid])/
            max(np.linalg.norm(amplitudes[1][valid]), 1e-300)))


def volume_modes(g, p, mesh, field, line_period):
    from .volume_sources import polarization_sources
    from .radiation import radiation_amplitude
    from .periodic_grating import rayleigh_orders, modal_power
    from .planar import fresnel
    m, tangent, kz = rayleigh_orders(g)
    valid = abs(kz.imag) < 1e-12
    directions = np.column_stack((np.zeros(valid.sum()), tangent[valid], kz[valid].real))
    _, rp, _ = fresnel(np.cos(np.deg2rad(p.incidence_deg)), p.eps_substrate)
    arrays = {}
    for order in (9, 11):
        points, weights, volumes = polarization_sources(mesh, field, p, False, order=order)
        expected = (g.period_um-g.ridge_width_um)*g.depth_um*line_period
        if not np.isclose(volumes['pit'], expected, rtol=1e-8):
            raise ValueError('Periodic trench contrast volume mismatch')
        far = radiation_amplitude(points, weights, directions, p)
        # Poisson-summed outgoing Green function: i/(2*A*k_z) exp(i*k.r).
        modal = np.zeros(len(m), complex)
        modal[valid] = np.sign(tangent[valid])*2j*np.pi/(line_period*g.period_um*p.k0*kz[valid])*far[:, 1]
        modal[m == 0] += rp
        arrays['volume_amplitude_q'+str(order)] = modal
        arrays['volume_efficiency_q'+str(order)] = modal_power(modal, kz, p.incidence_deg)
    return arrays


def energy_check(g, p, mesh, field, arrays, line_period, bottom=-.3):
    import ngsolve as ng
    from .solver import background_coefficient, evaluate_fields
    c = np.cos(np.deg2rad(p.incidence_deg))
    electric = field+background_coefficient(p)
    absorption = float(ng.Integrate(p.k0*p.eps_substrate.imag*ng.Norm(electric)**2,
        mesh, definedon=mesh.Materials('substrate'), order=12))/(line_period*g.period_um*c)
    xx, yy = np.meshgrid((np.arange(8)+.5)*line_period/8-line_period/2,
                        (np.arange(1024)+.5)*g.period_um/1024-g.period_um/2)
    points = np.column_stack((xx.ravel(), yy.ravel(), np.full(xx.size, bottom+1e-7)))
    e, h, _ = evaluate_fields(mesh, field, points, p)
    downward = float(-np.cross(e, h.conj())[:, 2].real.mean()/c)
    reflectance = float(arrays['efficiency'].sum())
    v = arrays['propagating']
    plane, volume = arrays['amplitude'][0, v], arrays['volume_amplitude_q11'][v]
    q9 = arrays['volume_amplitude_q9'][v]
    return dict(reflectance=reflectance, absorption_to_cut=absorption,
        downward_flux_at_cut=downward, energy_closure_error=1-reflectance-absorption-downward,
        volume_vs_plane_relative_l2=float(np.linalg.norm(volume-plane)/np.linalg.norm(plane)),
        volume_quadrature_relative_l2=float(np.linalg.norm(volume-q9)/np.linalg.norm(volume)))
