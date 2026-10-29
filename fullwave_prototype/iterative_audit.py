"""Validation-only BDDC/GMRES audit; never selected by the production UI."""
import argparse
from dataclasses import asdict, replace
import gc
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .parameters import Parameters


def assemble_problem(p):
    import ngsolve as ng
    from .solver import build_mesh, background_coefficient
    from .field_space import make_space
    mesh = build_mesh(p)
    space = make_space(mesh, p)
    u, v = space.TnT()
    eps = mesh.MaterialCF({'particle': p.eps_particle, 'substrate': p.eps_substrate,
                           'pml_substrate': p.eps_substrate}, default=1)
    a = ng.BilinearForm(space, symmetric=True, condense=True)
    a += (ng.curl(u)*ng.curl(v)-p.k0**2*eps*u*v)*ng.dx(bonus_intorder=p.integration_bonus)
    rhs = ng.LinearForm(space)
    rhs += p.k0**2*(p.eps_particle-1)*background_coefficient(p)*v*ng.dx(
        'particle', bonus_intorder=p.integration_bonus)
    return mesh, space, a, rhs


def reduced_rhs(a, rhs):
    result = rhs.vec.CreateVector()
    result.data = rhs.vec+a.harmonic_extension_trans*rhs.vec
    return result


def scipy_operator(matrix, template, free):
    from scipy.sparse.linalg import LinearOperator
    source, target = template.CreateVector(), template.CreateVector()
    def multiply(values):
        source.FV().NumPy()[:] = 0
        source.FV().NumPy()[free] = values
        target.data = matrix*source
        return target.FV().NumPy()[free].copy()
    return LinearOperator((int(free.sum()), int(free.sum())), matvec=multiply, dtype=complex)


def solve_iterative(p, output, reference, maxsteps, tolerance):
    import ngsolve as ng
    ng.SetNumThreads(p.threads)
    start = time.perf_counter()
    mesh, space, a, rhs = assemble_problem(p)
    free = np.array(list(space.FreeDofs(coupling=True)), dtype=bool)
    if space.ndof > p.max_dofs:
        raise ValueError('Audit exceeds the configured DOF limit')
    report = dict(parameters=asdict(p), dofs=space.ndof, free_dofs=int(free.sum()),
                  mode='validation-only; production UI unchanged', solver='NGSolve BDDC + SciPy restarted GMRES',
                  ngsolve_version=ng.__version__,
                  requested_solver_tolerance=tolerance, acceptance_true_residual=1e-8,
                  acceptance_direct_coefficient_l2=1e-5,
                  audit_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    direct = None
    if reference:
        with ng.TaskManager():
            a.Assemble()
            rhs.Assemble()
            reduced = reduced_rhs(a, rhs)
            inverse = a.mat.Inverse(space.FreeDofs(coupling=True), inverse='sparsecholesky')
            vector = rhs.vec.CreateVector()
            vector.data = inverse*reduced
            direct = vector.FV().NumPy().copy()
            np.savez_compressed(output/'direct_reference.npz', schur_coefficients=direct,
                                reduced_rhs=reduced.FV().NumPy().copy(), free_dofs=free)
        del inverse, vector, reduced
        gc.collect()
        print('Direct comparison coefficients saved in memory', flush=True)
    report['status'] = 'assembling_bddc'
    (output/'linear_progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    mesh.ngmesh.Save(str(output/'mesh.vol'))
    pre = ng.preconditioners.BDDC(a, coarsetype='direct', inverse='sparsecholesky')
    field = ng.GridFunction(space)
    with ng.TaskManager():
        stage = time.perf_counter()
        a.Assemble()
        rhs.Assemble()
        report['assembly_and_preconditioner_seconds'] = time.perf_counter()-stage
        print('BDDC assembled in {:.2f}s; starting GMRES'.format(
            report['assembly_and_preconditioner_seconds']), flush=True)
        def progress(iteration, residual):
            if iteration % 10 == 0:
                report.update(status='iterating', last_iteration=iteration,
                              last_preconditioned_residual=float(residual))
                (output/'linear_progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
                print('GMRES step {}: {:.6g}'.format(iteration, residual), flush=True)
        reduced = reduced_rhs(a, rhs)
        stage = time.perf_counter()
        from .iterative_linear import bounded_gmres
        solution, stats = bounded_gmres(scipy_operator(a.mat, rhs.vec, free),
            reduced.FV().NumPy()[free], scipy_operator(pre.mat, rhs.vec, free),
            maxsteps=maxsteps, tolerance=tolerance, progress=progress)
        field.vec.FV().NumPy()[free] = solution
        report.update(stats, iteration_seconds=time.perf_counter()-stage)
        residual = rhs.vec.CreateVector()
        residual.data = a.mat*field.vec-reduced
        report['true_free_residual'] = float(np.linalg.norm(residual.FV().NumPy()[free])/
                                            np.linalg.norm(reduced.FV().NumPy()[free]))
        if direct is not None:
            report['direct_coefficient_relative_l2'] = float(
                np.linalg.norm(field.vec.FV().NumPy()[free]-direct[free])/np.linalg.norm(direct[free]))
        valid = (report['scipy_info'] == 0 and report['true_free_residual'] < 1e-8 and
                 report.get('direct_coefficient_relative_l2', 0.) < 1e-5)
        report.update(status='linear_check_passed' if valid else 'linear_check_failed',
                      elapsed_seconds=time.perf_counter()-start)
        (output/'linear_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        if not valid:
            return None, None, report
        field.vec.data += a.harmonic_extension*field.vec
        field.vec.data += a.inner_solve*rhs.vec
    if not np.isfinite(field.vec.FV().NumPy()).all():
        raise RuntimeError('Non-finite recovered field; export rejected')
    return mesh, field, report


def export_sources(mesh, field, p, output, report):
    from .volume_sources import polarization_sources
    from .far_field import radiation_amplitude
    from .aperture_integrals import cap_quadrature
    points, weights, volumes = polarization_sources(mesh, field, p, True, 7)
    arrays = dict(source_points_um=points, polarization_weights_um3=weights,
                  field_coefficients=field.vec.FV().NumPy().copy())
    caps = []
    for index, tilt in enumerate((45., -45., 0.), 1):
        directions, measure = cap_quadrature(.5, tilt, 24)
        amplitude = radiation_amplitude(points, weights, directions, p)
        caps.append(float(measure@np.sum(abs(amplitude)**2, axis=1)))
        arrays.update({'ch{}_directions'.format(index): directions,
                       'ch{}_weights_sr'.format(index): measure,
                       'ch{}_amplitude_sp_um'.format(index): amplitude})
    if not np.isfinite(caps).all() or not all(np.isfinite(value).all() for value in arrays.values()):
        raise RuntimeError('Non-finite raw sources or radiation; export rejected')
    np.savez_compressed(output/'raw_sources.npz', **arrays)
    mesh.ngmesh.Save(str(output/'mesh.vol'))
    field.Save(str(output/'field.ngs'))
    report.update(cap_dsc_um2=caps, contrast_volumes_um3=volumes,
                  physical_validation='Not certified by a small algebraic residual')
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--direct-reference', action='store_true')
    parser.add_argument('--order', type=int, choices=(1, 2, 3))
    parser.add_argument('--maxsteps', type=int, default=300)
    parser.add_argument('--solver-tolerance', type=float, default=1e-10)
    args = parser.parse_args()
    saved = json.loads(args.config.read_text(encoding='utf-8'))
    p = Parameters(**saved.get('scattering_parameters', saved.get('parameters')))
    if args.order is not None:
        p = replace(p, order=args.order)
    p = replace(p, static_condensation=True)
    p.validate()
    if (p.defect_shape != 'sphere' or p.beam_profile != 'plane'
            or saved.get('roughness', {}).get('rms_nm', 0) != 0):
        raise ValueError('This first iterative audit is restricted to plane-wave spheres on a flat wafer')
    if not 1 <= args.maxsteps <= 1000:
        raise ValueError('Use 1..1000 GMRES steps')
    if not np.isfinite(args.solver_tolerance) or not 0 < args.solver_tolerance < 1:
        raise ValueError('Solver tolerance must be finite and in (0,1)')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Use a new output folder')
    args.output.mkdir(parents=True, exist_ok=True)
    mesh, field, report = solve_iterative(p, args.output, args.direct_reference, args.maxsteps,
                                        args.solver_tolerance)
    if field is None:
        print(json.dumps(report, indent=2))
        raise SystemExit(1)
    export_sources(mesh, field, p, args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
