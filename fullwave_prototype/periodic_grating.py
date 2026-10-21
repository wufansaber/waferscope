"""Independent 2D TM Maxwell benchmark for x-invariant periodic ridges.

NGSolve x is the physical transverse y coordinate; NGSolve y is height z.
H_x = exp(i*beta*y) h(y,z). This scalar reduction is exact only for
in-plane p illumination and an infinitely long, one-dimensional grating.
It is not the finite 3D particle/pattern solver.
"""
from dataclasses import dataclass, asdict
import time

import numpy as np


@dataclass(frozen=True)
class GratingInputs:
    wavelength_um: float = .266
    period_um: float = .5
    ridge_width_um: float = .25
    depth_um: float = .1
    substrate_n: float = 5.57
    substrate_k: float = .387
    incidence_deg: float = 75.

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite grating inputs required')
        if min(self.wavelength_um, self.period_um, self.ridge_width_um,
               self.depth_um, self.substrate_n) <= 0:
            raise ValueError('Positive grating dimensions/index required')
        if self.ridge_width_um > self.period_um or self.substrate_k < 0 or not 0 <= self.incidence_deg < 90:
            raise ValueError('Invalid ridge width, loss or incidence')

    @property
    def epsilon(self):
        return complex(self.substrate_n, self.substrate_k)**2


def rayleigh_orders(p, span=12):
    p.validate()
    m = np.arange(-span, span+1)
    tangent = np.sin(np.deg2rad(p.incidence_deg))+m*p.wavelength_um/p.period_um
    kz = np.sqrt(1-tangent**2+0j)
    return m, tangent, kz


def modal_power(amplitude, kz, incidence_deg):
    return np.maximum(np.real(kz), 0)/np.cos(np.deg2rad(incidence_deg))*abs(amplitude)**2


def build_mesh(p, factor=1., pml=.3, top=.4, bottom=.3):
    import ngsolve as ng
    from netgen.geom2d import SplineGeometry
    p.validate()
    if min(factor, pml, top, bottom-p.depth_um) <= 0:
        raise ValueError('Positive mesh/PML size and clearance required')
    # Preserve a nondegenerate grid in the flat-interface null test.
    width = min(p.ridge_width_um, .999*p.period_um)
    xs = np.array([-p.period_um/2, -width/2, width/2, p.period_um/2])
    zs = np.array([-bottom-pml, -bottom, -p.depth_um, 0, top, top+pml])
    geometry = SplineGeometry()
    points = [[geometry.AppendPoint(float(x), float(z)) for x in xs] for z in zs]
    domain = lambda ix, iz: 1+ix+3*iz if 0 <= ix < 3 and 0 <= iz < 5 else 0
    for iz in range(5):
        for ix in range(3):
            name = ('pml_substrate' if iz == 0 else 'substrate' if iz == 1 else
                    ('substrate' if ix == 1 or p.ridge_width_um == p.period_um else 'pit')
                    if iz == 2 else 'air' if iz == 3 else 'pml_air')
            number = domain(ix, iz)
            geometry.SetMaterial(number, name)
            geometry.SetDomainMaxH(number, factor*dict(substrate=.012, pit=.02,
                air=.035, pml_air=.045, pml_substrate=.018)[name])
    for iz in range(6):
        for ix in range(3):
            geometry.Append(['line', points[iz][ix], points[iz][ix+1]],
                leftdomain=domain(ix, iz), rightdomain=domain(ix, iz-1),
                bc='outer' if iz in (0, 5) else 'interface')
    for iz in range(5):
        left = geometry.Append(['line', points[iz][0], points[iz+1][0]],
            leftdomain=0, rightdomain=domain(0, iz), bc='left')
        for ix in (1, 2):
            geometry.Append(['line', points[iz][ix], points[iz+1][ix]],
                leftdomain=domain(ix-1, iz), rightdomain=domain(ix, iz), bc='interface')
        geometry.Append(['line', points[iz][3], points[iz+1][3]],
            leftdomain=domain(2, iz), rightdomain=0, bc='right', copy=left)
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=.045*factor))
    return mesh


def solve(p, order=3, factor=1., pml=.3, top=.4, bottom=.3, bonus=8):
    import ngsolve as ng
    from .planar import fresnel
    start = time.perf_counter()
    mesh = build_mesh(p, factor, pml, top, bottom)
    ng.SetNumThreads(2)
    space = ng.Periodic(ng.H1(mesh, order=order, complex=True, dirichlet='outer'))
    if space.ndof > 600000:
        raise ValueError('Periodic benchmark exceeds its 600000 DOF safety limit')
    k0 = 2*np.pi/p.wavelength_um
    c, s = np.cos(np.deg2rad(p.incidence_deg)), np.sin(np.deg2rad(p.incidence_deg))
    _, rp, q = fresnel(c, p.epsilon)
    beta = k0*s
    h_air = ng.exp(-1j*k0*c*ng.y)+complex(rp)*ng.exp(1j*k0*c*ng.y)
    dh_air = -1j*k0*c*ng.exp(-1j*k0*c*ng.y)+1j*k0*c*complex(rp)*ng.exp(1j*k0*c*ng.y)
    h_sub = (1+complex(rp))*ng.exp(-1j*k0*complex(q)*ng.y)
    dh_sub = -1j*k0*complex(q)*h_sub
    hbg = ng.IfPos(ng.y, h_air, h_sub)
    dhbg = ng.IfPos(ng.y, dh_air, dh_sub)
    inv_eps = mesh.MaterialCF({'substrate': 1/p.epsilon, 'pml_substrate': 1/p.epsilon}, default=1)
    u, v = space.TnT()
    du = ng.CF((ng.grad(u)[0]+1j*beta*u, ng.grad(u)[1]))
    dv = ng.CF((ng.grad(v)[0]-1j*beta*v, ng.grad(v)[1]))
    # Explicit diagonal coordinate stretch includes the Bloch carrier terms.
    stretch = mesh.MaterialCF({'pml_air': 1+2j, 'pml_substrate': 1+2j}, default=1)
    a = ng.BilinearForm(space)
    a += (inv_eps*(stretch*du[0]*dv[0]+du[1]*dv[1]/stretch)-
          k0*k0*stretch*u*v)*ng.dx(bonus_intorder=bonus)
    rhs = ng.LinearForm(space)
    if 'pit' in mesh.GetMaterials():
        rhs += (1/p.epsilon-1)*ng.CF((1j*beta*h_sub, dh_sub))*dv*ng.dx('pit', bonus_intorder=bonus)
    print('TM periodic H1 P{}: {} elements, {} DOFs'.format(order, mesh.ne, space.ndof), flush=True)
    field = ng.GridFunction(space)
    with ng.TaskManager():
        a.Assemble()
        rhs.Assemble()
        field.vec.data = a.mat.Inverse(space.FreeDofs(), inverse='umfpack')*rhs.vec
    residual = rhs.vec.CreateVector()
    residual.data = a.mat*field.vec-rhs.vec
    free = np.array(list(space.FreeDofs()), bool)
    rel = np.linalg.norm(residual.FV().NumPy()[free])/max(np.linalg.norm(rhs.vec.FV().NumPy()[free]), 1e-300)
    if not np.isfinite(rel) or rel > 1e-8:
        raise RuntimeError('Periodic linear solve failed residual check')
    return mesh, field, hbg, dhbg, inv_eps, dict(order=order, factor=factor, pml_um=pml,
        top_um=top, bottom_um=bottom, bonus=bonus, dofs=space.ndof, elements=mesh.ne,
        algebraic_residual=float(rel), seconds=time.perf_counter()-start)


def extract(p, mesh, field, hbg, dhbg, inv_eps, samples=4096, heights=(.15, .25), bottom=-.3+1e-7):
    import ngsolve as ng
    from .planar import fresnel
    k0 = 2*np.pi/p.wavelength_um
    c = np.cos(np.deg2rad(p.incidence_deg))
    _, rp, _ = fresnel(c, p.epsilon)
    m, tangent, kz = rayleigh_orders(p)
    x = (np.arange(samples)+.5)*p.period_um/samples-p.period_um/2
    amplitudes = []
    for z in heights:
        values = np.asarray(field(mesh(x, np.full(samples, z)))).ravel()
        coeff = np.exp(-2j*np.pi*m[:, None]*x[None, :]/p.period_um)@values/samples
        coeff *= np.exp(-1j*k0*kz*z)
        coeff[m == 0] += rp
        amplitudes.append(coeff)
    beta = k0*np.sin(np.deg2rad(p.incidence_deg))
    h = field+hbg
    ey = 1j*inv_eps*(ng.grad(field)[1]+dhbg)/k0
    ez = -1j*inv_eps*(ng.grad(field)[0]+1j*beta*h)/k0
    absorption_density = k0*p.epsilon.imag*(ng.Norm(ey)**2+ng.Norm(ez)**2)
    region = ng.IfPos(ng.y-bottom, 1, 0)
    absorption = float(ng.Integrate(absorption_density*region, mesh,
        definedon=mesh.Materials('substrate'), order=12))/(p.period_um*c)
    location = mesh(x, np.full(samples, bottom))
    downward = float(np.mean(np.real(np.asarray(ey(location)).ravel()*np.conj(np.asarray(h(location)).ravel())))/c)
    reflectance = modal_power(amplitudes[0], kz, p.incidence_deg)
    valid = abs(kz.imag) < 1e-12
    return dict(orders=m, tangent=tangent, kz=kz, amplitude=np.asarray(amplitudes),
        efficiency=reflectance, propagating=valid), dict(reflectance=float(reflectance.sum()),
        absorption_to_cut=absorption, downward_flux_at_cut=downward,
        energy_closure_error=float(1-reflectance.sum()-absorption-downward),
        modal_plane_relative_l2=float(np.linalg.norm((amplitudes[0]-amplitudes[1])[valid])/
            max(np.linalg.norm(amplitudes[1][valid]), 1e-300)),
        absorption_cut_z_um=bottom)
