"""Independent Mie near-field, optical theorem, PEC and vector-diffraction audit."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .parameters import Parameters
from .theory_reference import mie_fields, bessel_mie_cross_sections


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-300))


def sphere_check(p, output):
    import miepython
    import ngsolve as ng
    from .physical_audit import sphere_rule
    from .solver import solve, evaluate_fields, background_coefficient
    from .volume_sources import polarization_sources
    from .far_field import radiation_amplitude
    from .planar import background_fields, background_vectors
    if (1.7*p.diameter_nm+p.gap_nm)/1000 >= p.domain_radius_um:
        raise ValueError('Mie probe shells must fit entirely inside the physical domain')
    output.mkdir(parents=True, exist_ok=True)
    mesh, field, info = solve(p)
    directions, weights = sphere_rule(20)
    sources, moments, _ = polarization_sources(mesh, field, p, order=9)
    amplitude = radiation_amplitude(sources, moments, directions, p)
    csca = float(np.sum(weights*np.sum(abs(amplitude)**2, axis=-1)))
    ein, _ = background_fields(sources, p)
    cext = float(p.k0*np.imag(np.sum(ein.conj()*moments)))
    total = field+background_coefficient(p)
    cabs = float(p.k0*p.eps_particle.imag*ng.Integrate(ng.InnerProduct(total, total).real,
                  mesh, definedon=mesh.Materials('particle'), order=10))
    down, _, _, polarization, _, _ = background_vectors(p)
    aforward = radiation_amplitude(sources, moments, down[None], p)[0]
    s = np.cross((0, 0, 1), down)
    s = s/np.linalg.norm(s) if np.linalg.norm(s) > 1e-12 else np.array((-1., 0, 0))
    forward_cart = aforward[0]*s+aforward[1]*np.cross(down, s)
    optical = float(4*np.pi/p.k0*np.imag(polarization.conj()@forward_cart))
    exact = bessel_mie_cross_sections(p)
    qext, qsca, _, _ = miepython.efficiencies(complex(p.particle_n, -p.particle_k),
                                           p.diameter_nm/1000, p.wavelength_nm/1000)
    area = np.pi*(p.diameter_nm/2000)**2
    center = np.array((0, 0, (p.diameter_nm/2+p.gap_nm)/1000))
    radius = p.diameter_nm/2000
    radii = np.array((.6, 1.05, 1.5, 2.4))*radius
    points = np.concatenate([center+r*directions for r in radii])
    e, h, esc = evaluate_fields(mesh, field, points, p)
    exact_e, exact_h = mie_fields(points, p)
    ebg, _ = background_fields(points, p)
    near = []
    for i, r in enumerate(radii):
        sl = slice(i*len(directions), (i+1)*len(directions))
        near.append(dict(radius_nm=float(r*1000), region='inside' if r < radius else 'outside',
            electric_relative_l2=relative(e[sl], exact_e[sl]),
            magnetic_relative_l2=relative(h[sl], exact_h[sl]),
            scattered_e_relative_l2=relative(esc[sl], exact_e[sl]-ebg[sl]) if r > radius else None))
    result = dict(parameters=asdict(p), solver=info, csca_um2=csca, cext_um2=cext, cabs_um2=cabs,
        exact=exact, scattering_relative_error=abs(csca/exact['csca_um2']-1),
        extinction_relative_error=abs(cext/exact['cext_um2']-1),
        energy_balance_relative=abs(cext-csca-cabs)/max(abs(cext), 1e-300),
        optical_theorem_vs_volume_relative=abs(optical-cext)/max(abs(cext), 1e-300),
        miepython_vs_bessel_scattering_relative=abs(qsca*area/exact['csca_um2']-1),
        miepython_vs_bessel_extinction_relative=abs(qext*area/exact['cext_um2']-1),
        near_shells=near, miepython_version=miepython.__version__,
        note='No fitted phase, magnitude or position. Optical theorem and volume extinction share the same solved sources; near E/H oracle is independent.')
    np.savez_compressed(output/'fields.npz', points_um=points, e=e, hz0=h,
                        exact_e=exact_e, exact_hz0=exact_h, directions=directions, amplitude=amplitude)
    (output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('sphere', 'optics'), required=True)
    parser.add_argument('--output', default='build/theory_audit')
    parser.add_argument('--order', type=int, default=2)
    parser.add_argument('--bonus', type=int, default=0)
    parser.add_argument('--strength', type=float, default=2.)
    parser.add_argument('--pml-mesh-scale', type=float, default=1.)
    parser.add_argument('--pml-thickness', type=float, default=.16)
    parser.add_argument('--pml-order', type=int, default=0)
    parser.add_argument('--pml-inner-fraction', type=float, default=0)
    parser.add_argument('--loss', type=float, default=0.)
    parser.add_argument('--diameter-nm', type=float, default=100.)
    parser.add_argument('--wavelength-nm', type=float, default=355.)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    if args.case == 'sphere':
        p = Parameters(substrate_n=1, substrate_k=0, gap_nm=10,
            order=args.order, integration_bonus=args.bonus, static_condensation=True,
            direct_solver='sparsecholesky', pml_strength=args.strength, particle_k=args.loss,
            pml_mesh_scale=args.pml_mesh_scale, pml_thickness_um=args.pml_thickness,
            diameter_nm=args.diameter_nm, wavelength_nm=args.wavelength_nm, pml_order=args.pml_order,
            pml_inner_fraction=args.pml_inner_fraction)
        sphere_check(p, output)
    else:
        from .theory_optics import pec_image_check, vector_focus_check
        optics, arrays = vector_focus_check()
        result = dict(pec=pec_image_check(), vector_focus=optics)
        np.savez_compressed(output/'optics.npz', **arrays)
        (output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
