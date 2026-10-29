"""Validation-only SuperLU comparison; no production solver or UI changes."""
import argparse
from dataclasses import asdict, replace
import gc
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .parameters import Parameters
from .iterative_audit import assemble_problem, reduced_rhs, scipy_operator, export_sources


def solve(p, output, reference, ordering):
    import ngsolve as ng
    import scipy
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import splu
    from .superlu_linear import validated_matrix
    ng.SetNumThreads(p.threads)
    start = time.perf_counter()
    mesh, space, a, rhs = assemble_problem(p)
    if space.ndof > p.max_dofs:
        raise ValueError('Audit exceeds the configured DOF limit')
    free = np.array(list(space.FreeDofs(coupling=True)), dtype=bool)
    report = dict(parameters=asdict(p), dofs=space.ndof, free_dofs=int(free.sum()),
        mode='validation-only; production UI unchanged', solver='SciPy SuperLU',
        ordering=ordering, symmetric_mode=ordering == 'MMD_AT_PLUS_A', scipy_version=scipy.__version__,
        ngsolve_version=ng.__version__, acceptance_true_residual=1e-8,
        acceptance_direct_coefficient_l2=1e-5, audit_source_sha256=hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest())
    (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with ng.TaskManager():
        a.Assemble()
        rhs.Assemble()
        reduced = reduced_rhs(a, rhs)
        report['mesh_and_assembly_seconds'] = time.perf_counter()-start
        direct = None
        if reference:
            stage = time.perf_counter()
            inverse = a.mat.Inverse(space.FreeDofs(coupling=True), inverse='sparsecholesky')
            vector = rhs.vec.CreateVector()
            vector.data = inverse*reduced
            direct = vector.FV().NumPy()[free].copy()
            report['reference_factor_solve_seconds'] = time.perf_counter()-stage
            np.savez_compressed(output/'direct_reference.npz', free_coefficients=direct, free_dofs=free)
            del inverse, vector
            gc.collect()
        entries = [np.array(v.NumPy()) for v in a.mat.CSR()]
        csr = csr_matrix(tuple(entries), shape=(space.ndof, space.ndof))[free][:, free]
        del entries
        operator = scipy_operator(a.mat, rhs.vec, free)
        matrix, check = validated_matrix(csr, operator.matvec)
        del csr
        report.update(check, matrix_nonzeros=int(matrix.nnz), status='factorizing')
        (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print('SuperLU factorization: {} free DOFs; {} nonzeros'.format(free.sum(), matrix.nnz), flush=True)
        stage = time.perf_counter()
        lu = splu(matrix, permc_spec=ordering, diag_pivot_thresh=1.,
                  options={'SymmetricMode': ordering == 'MMD_AT_PLUS_A'})
        solution = lu.solve(reduced.FV().NumPy()[free])
        report.update(factor_solve_seconds=time.perf_counter()-stage,
                      lu_nonzeros=int(lu.nnz))
        field = ng.GridFunction(space)
        field.vec.FV().NumPy()[free] = solution
        residual = rhs.vec.CreateVector()
        residual.data = a.mat*field.vec-reduced
        report['true_free_residual'] = float(np.linalg.norm(residual.FV().NumPy()[free])/
            max(np.linalg.norm(reduced.FV().NumPy()[free]), 1e-300))
        if direct is not None:
            report['direct_coefficient_relative_l2'] = float(np.linalg.norm(solution-direct)/np.linalg.norm(direct))
        valid = (np.isfinite(solution).all() and report['true_free_residual'] < 1e-8
                 and report.get('direct_coefficient_relative_l2', 0.) < 1e-5)
        report.update(status='linear_check_passed' if valid else 'linear_check_failed',
                      elapsed_seconds=time.perf_counter()-start)
        (output/'linear_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        if not valid:
            raise RuntimeError('SuperLU failed the native-operator residual/reference check')
        del lu, matrix, direct
        field.vec.data += a.harmonic_extension*field.vec
        field.vec.data += a.inner_solve*rhs.vec
    if not np.isfinite(field.vec.FV().NumPy()).all():
        raise RuntimeError('Non-finite recovered field; export rejected')
    export_sources(mesh, field, p, output, report)
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--direct-reference', action='store_true')
    parser.add_argument('--pml-mesh-scale', type=float)
    parser.add_argument('--ordering', choices=('COLAMD', 'MMD_AT_PLUS_A'), default='COLAMD')
    args = parser.parse_args()
    saved = json.loads(args.config.read_text(encoding='utf-8'))
    p = Parameters(**saved.get('scattering_parameters', saved.get('parameters')))
    p = replace(p, static_condensation=True, direct_solver='sparsecholesky',
                pml_mesh_scale=args.pml_mesh_scale if args.pml_mesh_scale is not None else p.pml_mesh_scale)
    p.validate()
    if p.defect_shape != 'sphere' or p.beam_profile != 'plane' or saved.get('roughness', {}).get('rms_nm', 0) != 0:
        raise ValueError('This audit only covers plane-wave spheres on a flat wafer')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Use a new output folder')
    args.output.mkdir(parents=True, exist_ok=True)
    solve(p, args.output, args.direct_reference, args.ordering)


if __name__ == '__main__':
    main()
