"""Vector H(curl) scattered-field solve with a curved sphere and radial PML."""
import time
import numpy as np

from .planar import background_fields
from .absorbing_boundary import set_boundary


def build_mesh(p):
    import ngsolve as ng
    from netgen.occ import Sphere, Box, Glue, OCCGeometry
    outer_r = p.domain_radius_um + p.pml_thickness_um
    outer = Sphere((0, 0, 0), outer_r)
    outer.faces.name = 'outer'
    inner = Sphere((0, 0, 0), p.domain_radius_um)
    lower = Box((-2*outer_r, -2*outer_r, -2*outer_r), (2*outer_r, 2*outer_r, 0))
    ball = Sphere((0, 0, (p.diameter_nm/2+p.gap_nm)/1000), p.diameter_nm/2000)
    ball.mat('particle')
    ball.maxh = p.particle_mesh_um
    ball.faces.maxh = p.particle_mesh_um
    regions = [(inner-lower-ball, 'air', p.air_mesh_um),
               (inner*lower, 'substrate', p.substrate_mesh_um),
               ((outer-inner)-lower, 'pml_air', p.air_mesh_um*p.pml_mesh_scale),
               ((outer-inner)*lower, 'pml_substrate', p.substrate_mesh_um*p.pml_mesh_scale)]
    for shape, name, h in regions:
        shape.mat(name)
        shape.maxh = h
    mesh = ng.Mesh(OCCGeometry(Glue([s for s, _, _ in regions]+[ball])).GenerateMesh(
        maxh=p.air_mesh_um))
    mesh.Curve(max(3, p.order+1))
    set_boundary(mesh, p)
    return mesh


def contrast_materials(p, particle_enabled=True):
    return {'particle': p.eps_particle-1 if particle_enabled else 0,
            'bump': p.eps_substrate-1, 'pit': 1-p.eps_substrate}


def background_coefficient(p):
    import ngsolve as ng
    from .illumination import illumination_rays
    pos = ng.CF((ng.x, ng.y, ng.z))
    above, below = ng.CF((0j, 0j, 0j)), ng.CF((0j, 0j, 0j))
    for down, up, kt, ein, er, et in illumination_rays(p):
        above += ng.CF(tuple(ein))*ng.exp(1j*p.k0*(ng.CF(tuple(down))*pos))
        above += ng.CF(tuple(er))*ng.exp(1j*p.k0*(ng.CF(tuple(up))*pos))
        below += ng.CF(tuple(et))*ng.exp(1j*p.k0*(ng.CF(tuple(kt))*pos))
    return ng.IfPos(ng.z, above, below)


def solve(p, mesh=None, particle_enabled=True):
    import ngsolve as ng
    p.validate()
    ng.SetNumThreads(p.threads)
    print('NGSolve worker threads: {}'.format(p.threads), flush=True)
    start = time.perf_counter()
    if mesh is None:
        print('Generating curved particle/substrate mesh...', flush=True)
        mesh = build_mesh(p)
    eps = mesh.MaterialCF({'particle': p.eps_particle if particle_enabled else 1,
                           'substrate': p.eps_substrate, 'bump': p.eps_substrate,
                           'pml_substrate': p.eps_substrate}, default=1)
    from .field_space import make_space
    space = make_space(mesh, p)
    print('Elements: {}; degrees of freedom: {}'.format(mesh.ne, space.ndof), flush=True)
    if space.ndof > p.max_dofs:
        raise ValueError('FEM size exceeds max_dofs; coarsen the mesh or explicitly raise the limit.')
    e_bg = background_coefficient(p)
    u, v = space.TnT()
    a = ng.BilinearForm(space, symmetric=True, condense=p.static_condensation)
    a += (ng.curl(u)*ng.curl(v)-p.k0**2*eps*u*v)*ng.dx(bonus_intorder=p.integration_bonus)
    rhs = ng.LinearForm(space)
    for material, contrast in contrast_materials(p, particle_enabled).items():
        if material in mesh.GetMaterials() and contrast != 0:
            rhs += p.k0**2*contrast*e_bg*v*ng.dx(material, bonus_intorder=p.integration_bonus)
    field = ng.GridFunction(space)
    with ng.TaskManager():
        assembly_start = time.perf_counter()
        a.Assemble()
        rhs.Assemble()
        assembly_seconds = time.perf_counter()-assembly_start
        print('Assembly complete in {:.1f}s; factoring with {}...'.format(
            assembly_seconds, p.direct_solver), flush=True)
        solve_rhs = rhs.vec.CreateVector()
        solve_rhs.data = rhs.vec
        if p.static_condensation:
            solve_rhs.data += a.harmonic_extension_trans*rhs.vec
        factor_start = time.perf_counter()
        inverse = a.mat.Inverse(space.FreeDofs(coupling=p.static_condensation), inverse=p.direct_solver)
        factor_seconds = time.perf_counter()-factor_start
        print('Factorization complete in {:.1f}s; recovering the full field...'.format(factor_seconds), flush=True)
        field.vec.data = inverse*solve_rhs
        residual = rhs.vec.CreateVector()
        residual.data = a.mat*field.vec-solve_rhs
        if p.static_condensation:
            field.vec.data += a.harmonic_extension*field.vec
            field.vec.data += a.inner_solve*rhs.vec
    free = np.array(list(space.FreeDofs(coupling=p.static_condensation)), dtype=bool)
    rel = np.linalg.norm(residual.FV().NumPy()[free])/max(np.linalg.norm(solve_rhs.FV().NumPy()[free]), 1e-300)
    if not np.isfinite(rel) or rel > 1e-8 or not np.isfinite(field.vec.FV().NumPy()).all():
        raise RuntimeError('Linear solve failed the true algebraic residual/finite-field check: {}'.format(rel))
    volume = ng.Integrate(1, mesh, definedon=mesh.Materials('particle'))
    metadata = dict(elements=mesh.ne, dofs=space.ndof, algebraic_residual=float(rel),
                    threads=p.threads,
                    assembly_seconds=assembly_seconds, factor_seconds=factor_seconds,
                    factored_dofs=sum(space.FreeDofs(coupling=p.static_condensation)),
                    static_condensation=p.static_condensation,
                    direct_solver=p.direct_solver,
                    integration_bonus=p.integration_bonus,
                    pml_order=p.pml_order or p.order,
                    residual_system='Schur complement' if p.static_condensation else 'full Maxwell matrix',
                    solve_seconds=time.perf_counter()-start, particle_volume_um3=float(volume),
                    ngsolve_version=ng.__version__)
    print('Maxwell solve complete in {:.1f}s; relative residual {:.3g}'.format(
        metadata['solve_seconds'], rel), flush=True)
    return mesh, field, metadata


def evaluate_fields(mesh, field, points, p):
    import ngsolve as ng
    points = np.asarray(points)
    locations = mesh(points[:, 0], points[:, 1], points[:, 2])
    esc = np.asarray(field(locations)).reshape(-1, 3)
    hsc = np.asarray(ng.curl(field)(locations)).reshape(-1, 3)/(1j*p.k0)
    ebg, hbg = background_fields(points, p)
    return esc+ebg, hsc+hbg, esc
