"""Reproducible physical checks; no 'validated' flag inferred from solver residuals."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss

from .parameters import Parameters
from .solver import solve, evaluate_fields, background_coefficient
from .volume_sources import polarization_sources
from .far_field import radiation_amplitude
from .planar import background_fields


def sphere_rule(n=24, hemisphere=False):
    mu, weights = leggauss(n)
    if hemisphere:
        mu, weights = (mu+1)/2, weights/2
    phi = np.arange(2*n)*np.pi/n
    mm, pp = np.meshgrid(mu, phi, indexing='ij')
    rho = np.sqrt(1-mm*mm)
    directions = np.stack((rho*np.cos(pp), rho*np.sin(pp), mm), axis=-1).reshape(-1, 3)
    return directions, np.repeat(weights*np.pi/n, 2*n)


def reference(p, output, python, kind, directions, weights):
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'directions.npz', directions=directions, weights=weights)
    payload = dict(parameters=asdict(p), kind=kind, directions_file=str((output/'directions.npz').resolve()))
    (output/'request.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    subprocess.run([python, '-m', 'fullwave_prototype.audit_reference', '--request',
                    str(output/'request.json'), '--output', str(output/'reference.npz')], check=True)
    with np.load(output/'reference.npz') as data:
        return dict(data)


def free_sphere(p, output, python):
    import ngsolve as ng
    mesh, field, metrics = solve(p)
    points, sources, _ = polarization_sources(mesh, field, p, order=7)
    directions, weights = sphere_rule()
    amplitude = radiation_amplitude(points, sources, directions, p)
    csca = float(np.sum(weights*np.sum(abs(amplitude)**2, axis=-1)))
    ebg, _ = background_fields(points, p)
    cext = float(p.k0*np.imag(np.sum(np.conj(ebg)*sources)))
    total = field+background_coefficient(p)
    cabs = float(p.k0*p.eps_particle.imag*ng.Integrate(
        ng.InnerProduct(total, total).real, mesh, definedon=mesh.Materials('particle'), order=8))
    ref = reference(p, output, python, 'mie', directions, weights)
    shells = []
    # Surface flux is evaluated from native curl(E), independently of the far-field extractor.
    sd, sw = sphere_rule(48)
    for radius in (.15, .19):
        e, h, esc = evaluate_fields(mesh, field, radius*sd, p)
        _, hb = background_fields(radius*sd, p)
        sc = float(np.sum(sw*np.real(np.sum(np.cross(esc, np.conj(h-hb))*sd, axis=1)))*radius**2)
        net = float(np.sum(sw*np.real(np.sum(np.cross(e, np.conj(h))*sd, axis=1)))*radius**2)
        shells.append(dict(radius_um=radius, scattered_flux_um2=sc, total_outward_flux_um2=net,
                           scattered_vs_far_relative=abs(sc/csca-1),
                           absorption_balance_relative=abs(net+cabs)/max(csca+cabs, 1e-300)))
    report = dict(parameters=asdict(p), solver=metrics, csca_um2=csca, cext_um2=cext, cabs_um2=cabs,
                  mie_csca_um2=float(ref['csca_um2']), mie_cext_um2=float(ref['cext_um2']),
                  scattering_relative_error=abs(csca/float(ref['csca_um2'])-1),
                  extinction_relative_error=abs(cext/float(ref['cext_um2'])-1),
                  energy_balance_relative=abs(cext-csca-cabs)/max(cext, 1e-300), flux_shells=shells)
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def flat_pupil(p, output, python):
    from optical_propagation import pupil_geometry
    from .imaging import ImagingSettings, image_from_amplitude
    mesh, field, info = solve(p)
    print('Extracting complex volume sources for the {} pupil benchmark...'.format(p.polarization), flush=True)
    points, sources, _ = polarization_sources(mesh, field, p, order=7)
    settings = ImagingSettings(pupil_samples=32, padding=4)
    geometries = [pupil_geometry(.5, tilt, 32) for tilt in settings.tilt_degrees]
    directions = np.concatenate([k[mask] for _, k, _, mask, _ in geometries])
    ref = reference(p, output, python, 'bv', directions, np.ones(len(directions)))['amplitude']
    print('Independent BV reference complete; evaluating FEM pupils...', flush=True)
    # BV uses a sphere-centred outgoing expansion; FEM uses the interface origin.
    height = (p.diameter_nm/2+p.gap_nm)/1000
    corrected = ref*np.exp(-1j*p.k0*height*directions[:, 2])[:, None]
    actual = radiation_amplitude(points, sources, directions, p)
    phase = np.vdot(corrected, actual)
    phase /= abs(phase)
    results, arrays, start = [], {}, 0
    for index, (local, k, rotation, mask, du) in enumerate(geometries):
        count = int(mask.sum())
        aa, bb = np.zeros(mask.shape+(2,), complex), np.zeros(mask.shape+(2,), complex)
        aa[mask], bb[mask] = actual[start:start+count], corrected[start:start+count]*phase
        start += count
        a = image_from_amplitude(aa, local, k, rotation, du, p.wavelength_nm/1000, settings)
        b = image_from_amplitude(bb, local, k, rotation, du, p.wavelength_nm/1000, settings)
        results.append(dict(tilt_deg=settings.tilt_degrees[index],
            pupil_complex_relative_l2=float(np.linalg.norm(aa-bb)/np.linalg.norm(bb)),
            image_flux_relative_l2=float(np.linalg.norm(a['flux']-b['flux'])/np.linalg.norm(b['flux'])),
            collected_power_relative_error=abs(a['collected_dsc_um2']/b['collected_dsc_um2']-1)))
        arrays['ch{}_fem'.format(index+1)], arrays['ch{}_bv'.format(index+1)] = a['flux'], b['flux']
        arrays['axis_um'] = a['axis_um']
    np.savez_compressed(output/'images.npz', **arrays)
    report = dict(parameters=asdict(p), solver=info, channels=results,
                  fitted_global_phase_rad=float(np.angle(phase)),
                  note='Known origin translation plus ONE constant unit phase for all channels; no amplitude scaling or spatial shift fitted. Imaging propagator is shared, so this checks the FEM pupil against BV, not independent lens theory.')
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def rough_pml(output, order=2, pml_power=1, pml_mesh_scale=1., static_condensation=False,
              direct_solver='umfpack', pml_strength=2., integration_bonus=0):
    import ngsolve as ng
    from .rough_surface import Roughness, build_rough_mesh
    from .absorbing_boundary import set_boundary
    p = Parameters(gap_nm=10, order=order, pml_power=pml_power, pml_mesh_scale=pml_mesh_scale,
                   static_condensation=static_condensation, direct_solver=direct_solver,
                   integration_bonus=integration_bonus)
    p.validate()
    roughness = Roughness()
    ng.SetNumThreads(p.threads)
    mesh = build_rough_mesh(p, roughness)
    directions, weights = sphere_rule(20, hemisphere=True)
    probes = np.array(((0, 0, .007), (0, .02, .012), (0, .08, .025), (0, .06, -.02)))
    results, arrays = [], {}
    for strength in (1., 2., 3.):
        case = replace(p, pml_strength=strength)
        set_boundary(mesh, case)
        _, field, info = solve(case, mesh=mesh)
        points, sources, _ = polarization_sources(mesh, field, case)
        a = radiation_amplitude(points, sources, directions, case)
        e, _, _ = evaluate_fields(mesh, field, probes, case)
        results.append(dict(strength=strength, solver=info,
                            upper_dsc_um2=float(np.sum(weights*np.sum(abs(a)**2, axis=-1)))))
        arrays['amplitude_'+str(int(strength))], arrays['near_e_'+str(int(strength))] = a, e
        del field
    for item in results:
        key = str(int(item['strength']))
        item['power_change_vs_strength2'] = abs(item['upper_dsc_um2']/results[1]['upper_dsc_um2']-1)
        item['near_probe_relative_l2_vs_strength2'] = float(np.linalg.norm(arrays['near_e_'+key]-arrays['near_e_2'])/np.linalg.norm(arrays['near_e_2']))
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'fields.npz', **arrays, probes_um=probes)
    report = dict(parameters=asdict(p), roughness=asdict(roughness), cases=results,
                  note='Same geometry and mesh; only PML strength changes. Not domain-size or thickness convergence. Particle-present only.')
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='build/physical_audit')
    parser.add_argument('--reference-python', default='D:/anaconda/envs/NIST2/python.exe')
    parser.add_argument('--case', choices=('mie', 'flat', 'pml'), required=True)
    parser.add_argument('--order', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--pml-power', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--pml-mesh-scale', type=float, default=1.)
    parser.add_argument('--pml-strength', type=float, default=2.)
    parser.add_argument('--integration-bonus', type=int, default=0)
    parser.add_argument('--static-condensation', action='store_true')
    parser.add_argument('--direct-solver', choices=('umfpack', 'sparsecholesky'), default='umfpack')
    parser.add_argument('--loss', type=float, nargs='+', default=(0., .1), help='Particle k values for Mie checks')
    args = parser.parse_args()
    if args.case == 'pml' and args.pml_strength != 2:
        parser.error('The PML audit uses fixed strengths 1, 2, 3; omit --pml-strength for this case.')
    root = Path(args.output)
    controls = dict(order=args.order, pml_power=args.pml_power, pml_mesh_scale=args.pml_mesh_scale,
                    pml_strength=args.pml_strength, integration_bonus=args.integration_bonus,
                    static_condensation=args.static_condensation, direct_solver=args.direct_solver)
    if args.case == 'mie':
        for loss in args.loss:
            p = Parameters(gap_nm=10, substrate_n=1, substrate_k=0, particle_k=loss, **controls)
            print(json.dumps(free_sphere(p, root/('mie_'+str(loss)), args.reference_python), indent=2), flush=True)
    elif args.case == 'flat':
        for pol in ('p', 's'):
            print(json.dumps(flat_pupil(Parameters(gap_nm=10, polarization=pol, **controls),
                                       root/('flat_'+pol), args.reference_python), indent=2), flush=True)
    else:
        print(json.dumps(rough_pml(root/'rough_pml', **controls), indent=2), flush=True)


if __name__ == '__main__':
    main()
